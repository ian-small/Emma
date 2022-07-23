# ARGS[1] = genome file(fasta)
# ARGS[2] = annotation file(gff3)

include("circularity.jl")
include("orfs.jl")

function getcodons(seq::CircularSequence, pattern)
    frames = [NamedTuple{(:pos,:codon), Tuple{Int32, Codon}}[] for f in 1:3]
    for m in eachmatch(pattern, seq.sequence)
        i::Int32 = m.captured[1]
        push!(frames[mod1(i, 3)], (pos=i, codon=matched(m)))
    end
    return frames
end

struct GFF
    ftype::String
    fstart::Int32
    fend::Int32
    strand::Char
    name::String
end

function readGFF(file::String, glength::Int32)
    features = []
    open(file, "r") do infile
        for line in eachline(infile)
            bits = split(line, "\t")
            strand = bits[7][1]
            fstart = parse(Int32, bits[4])
            fend = parse(Int32, bits[5])
            if strand == '-'
                tmp = glength - fstart + 1
                fstart = glength - fend + 1
                fend = tmp
            end
            attributes = split(bits[9], ";")
            name = first(attributes)[6:end]
            push!(features, GFF(bits[3], fstart, fend, strand, name))
        end
    end
    return features
end


function main()
    genome = FASTA.Record()
    reader = open(FASTA.Reader, ARGS[1])
    read!(reader, genome)
    close(reader)

    fseq = CircularSequence(FASTA.sequence(genome))
    rseq = reverse_complement(fseq)

    fstarts = getcodons(fseq, startcodon)
    rstarts = getcodons(rseq, startcodon)

    fstops = codonmatches(fseq, stopcodon)
    rstops = codonmatches(rseq, stopcodon)

    genes = readGFF(ARGS[2], length(fseq))
    filter!(x->x.ftype≠"ORF", genes)

    fgenes = filter(x->x.strand=='+', genes)
    rgenes = filter(x->x.strand=='-', genes)

    fcds = filter(x->x.ftype=="HMM alignment", fgenes)
    rcds = filter(x->x.ftype=="HMM alignment", rgenes)

    for cds in fcds

    end

end

main()