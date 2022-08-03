include("circularity.jl")
include("orfs.jl")
include("tRNAs.jl")
include("rRNAs.jl")
include("gff.jl")
include("visuals.jl")

using XGBoost

const minORF = 150

function get_overlapped_trns(trns::Vector{CMAlignment_trn}, glength::Integer)
    overlapped = Vector{Tuple{CMAlignment_trn, Int}}(undef, 0)
    for i in 1:length(trns)-1
        if trns[i].tto >= trns[i+1].tfrom
            push!(overlapped, (trns[i], trns[i+1].tfrom - 1))
        end
    end
    if last(trns).tto >= glength + first(trns).tfrom
        push!(overlapped, (last(trns), glength + first(trns).tfrom - 1))
    end
    return overlapped
end

function main(infile::String, outfile::String, svgfile::String)
    target = FASTA.Record()
    reader = open(FASTA.Reader, infile)
    read!(reader, target)
    id = FASTA.identifier(target)
    genome = CircularSequence(FASTA.sequence(LongDNA{4}, target))
    rev_genome = reverse_complement(genome)
    println(id, "\t", length(genome))

    #extend genome
    extended_genome = genome[1:length(genome)+100]
    writer = open(FASTA.Writer, "tmp.extended.fa")
    write(writer, FASTA.Record(id, extended_genome))
    close(writer)

    #find tRNAs
    trn_matches = parse_trn_alignments(cmsearch("trn", "all_trn.cm"), length(genome))
    filter!(x -> isequal(x.query, get(anticodon2trn, x.anticodon, "")), trn_matches)
    filter!(x -> x.Evalue < 1e-4, trn_matches)
    filter!(x -> x.tfrom <= length(genome), trn_matches)
    #check for overlapping tRNAs
    overlapped = get_overlapped_trns(sort(filter(x -> x.tstrand == '+', trn_matches), by=x->x.tfrom), length(genome))
    append!(overlapped, get_overlapped_trns(sort(filter(x -> x.tstrand == '-', trn_matches), by=x->x.tfrom), length(genome)))
    println(overlapped)
    #for overlapped tRNAs, generate polyadenylated version
    for (cma, trunc_end) in overlapped
        trnseq = cma.tstrand == '+' ? genome.sequence[cma.tfrom:trunc_end] : rev_genome.sequence[cma.tfrom:trunc_end]
        trnseq_polyA = trnseq * dna"AAAAAAAAAA"
        polyA_matches = parse_trn_alignments(cmsearch(cma.query, "trn", trnseq_polyA), 0)
        isempty(polyA_matches) && continue
        trn_match = polyA_matches[1]
        if trn_match.Evalue < cma.Evalue
            #construct modified CMA with new Evalue, new end point, polyA flag
            newcma = CMAlignment_trn(cma.query, cma.target, trn_match.Evalue,
            trn_match.qfrom, trn_match.qto, cma.tfrom, trunc_end, cma.tstrand,
            trn_match.tseq, trn_match.anticodon, (trn_match.tto - trn_match.tfrom) - (trunc_end - cma.tfrom))
            #delete old match from trn_matches
            deleteat!(trn_matches, findfirst(x -> x == cma, trn_matches))
            #add new CMA
            push!(trn_matches, newcma)
        end
    end
    sort!(trn_matches, by=x->x.tfrom)
    println(trn_matches)

    #find rRNAs
    rRNAs = (rrnL = rRNA(Vector{nHMMmatch}(undef, 0), Vector{CMAlignment_rrn}(undef, 0)),
            rrnS = rRNA(Vector{nHMMmatch}(undef, 0), Vector{CMAlignment_rrn}(undef, 0)))
    #search for rrn ends using cmsearch
    rrn_matches = parse_rrn_alignments(cmsearch("rrn", "all_rrn.cm"), length(genome))
    filter!(x -> x.tfrom <= length(genome), rrn_matches)
    for m in rrn_matches
        if m.query == "rrnL"
            push!(rRNAs.rrnL.stop, m)
        elseif m.query == "rrnS"
            push!(rRNAs.rrnS.stop, m)
        end
    end
    #search for rrn starts using hmmsearch
    rrn_starts = parse_tbl(rrnsearch())
    filter!(x -> x.ali_from <= length(genome), rrn_starts)
    for m in rrn_starts
        if m.query == "rrnL"
            push!(rRNAs.rrnL.start, m)
        elseif m.query == "rrnS"
            push!(rRNAs.rrnS.start, m)
        end
    end
    #fix ends using flanking trn genes

    println(rRNAs)

    #find CDSs
    fstarts, fstartcodons = getcodons(genome, startcodon)
    fstops = codonmatches(genome, stopcodon)
    rstarts, rstartcodons = getcodons(rev_genome, startcodon)
    rstops = codonmatches(rev_genome, stopcodon)
    cds_matches = parse_domt(orfsearch(id, genome, fstarts, fstops, rstarts, rstops, minORF), length(genome))
    #rationalise HMM matches to leave one per ORF
    #sort by ORF, HMM, genome position
    sort!(cds_matches, by = x -> (x.orf, x.query, x.ali_from))
    #merge adjacent matches to same HMM
    #drop HMM matches that score lower than other overlapping matches
    rationalised_cds_matches = HMMmatch[]
    current_cds_match = first(cds_matches)
    for i in 2:length(cds_matches)
        cds_match = cds_matches[i]
        if cds_match.orf == current_cds_match.orf
            if cds_match.query == current_cds_match.query
                current_cds_match = merge_hmm_matches(current_cds_match, cds_match)
            elseif cds_match.Evalue < current_cds_match.Evalue
                current_cds_match = cds_match
            end
        else
            push!(rationalised_cds_matches, current_cds_match)
            current_cds_match = cds_match
        end
    end
    push!(rationalised_cds_matches, current_cds_match)
    #fix start & stop codons
    #load XGBoost model
    startcodon_model = Booster(model_file=joinpath(emmamodels, "xgb.model"))
    ftrns = filter(x->x.tstrand == '+', trn_matches)
    rtrns = filter(x->x.tstrand == '-', trn_matches)
    fhmms = filter(x->x.strand == '+', rationalised_cds_matches)
    rhmms = filter(x->x.strand == '-', rationalised_cds_matches)
    fix_start_and_stop_codons!(fhmms, ftrns, fstarts, fstartcodons, fstops, startcodon_model, length(genome))
    fix_start_and_stop_codons!(rhmms, rtrns, rstarts, rstartcodons, rstops, startcodon_model, length(genome))
    gffs = writeGFF(outfile, id, length(genome), append!(fhmms,rhmms), trn_matches, rRNAs)
    drawgenome(svgfile, id, length(genome), gffs)
end

main(ARGS[1], ARGS[2], ARGS[3])

#ARGS[1] = fasta input
#ARGS[2] = gff output
#ARGS[3] = svg output