include("circularity.jl")
include("orfs.jl")
include("tRNAs.jl")
include("rRNAs.jl")
include("gff.jl")

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

function main(infile::String, outfile::String)
    target = FASTA.Record()
    reader = open(FASTA.Reader, infile)
    read!(reader, target)
    id = FASTA.identifier(target)
    genome = CircularSequence(FASTA.sequence(LongDNA{4}, target))
    rev_genome_seq = BioSequences.reverse_complement(genome.sequence)
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
        trnseq = cma.tstrand == '+' ? genome.sequence[cma.tfrom:trunc_end] : rev_genome_seq[cma.tfrom:trunc_end]
        trnseq_polyA = trnseq * dna"AAAAAAAAAA"
        trn_match = parse_trn_alignments(cmsearch(cma.query, "trn", trnseq_polyA), 0)[1]
        if trn_match.Evalue < cma.Evalue
            #construct modified CMA with new Evalue, new end point, polyA flag
            newcma = CMAlignment_trn(cma.query, cma.target, trn_match.Evalue,
            trn_match.qfrom, trn_match.qto, cma.tfrom, trunc_end, cma.tstrand,
            trn_match.tseq, trn_match.anticodon, true)
            #delete old match from trn_matches
            deleteat!(trn_matches, findfirst(x -> x == cma, trn_matches))
            #add new CMA
            push!(trn_matches, newcma)
        end
    end
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
    cds_matches = parse_domt(orfsearch(id, genome, minORF))
    filter!(x -> x.Evalue < 1e-10, cds_matches)

    #fix start codons

    #fix stop codons

    println(cds_matches)

    writeGFF(outfile, id, length(genome), cds_matches, trn_matches, rRNAs)
end

main(ARGS[1], ARGS[2])

#ARGS[1] = fasta input