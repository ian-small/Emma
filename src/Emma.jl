include("circularity.jl")
include("orfs.jl")
include("tRNAs.jl")
include("gff.jl")

const minORF = 150

function get_overlapped_trns(trns::Vector{CMAlignment}, glength::Integer)
    overlapped = Vector{Tuple{CMAlignment, Int}}(undef, 0)
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
    #find tRNAs
    trn_matches = parse_alignments(cmsearch(id, genome), length(genome))
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
        trn_match = parse_alignments(cmsearch(cma.query, trnseq_polyA), 0)[1]
        if trn_match.Evalue < cma.Evalue
            #construct modified CMA with new Evalue, new end point, polyA flag
            newcma = CMAlignment(cma.query, cma.target, trn_match.Evalue,
            trn_match.qfrom, trn_match.qto, cma.tfrom, trunc_end, cma.tstrand,
            trn_match.tseq, trn_match.anticodon, true)
            #delete old match from trn_matches
            deleteat!(trn_matches, findfirst(x -> x == cma, trn_matches))
            #add new CMA
            push!(trn_matches, newcma)
        end
    end
    println(trn_matches)
    cds_matches = parse_domt(hmmsearch(id, genome, minORF))
    filter!(x -> x.Evalue < 1e-10, cds_matches)
    println(cds_matches)

    writeGFF(outfile, id, length(genome), cds_matches, trn_matches)
end

main(ARGS[1], ARGS[2])

#ARGS[1] = fasta input