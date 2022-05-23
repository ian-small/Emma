using Artifacts
using BioSequences
using FASTX

const emmamodels = joinpath(artifact"Emma-models", "emma-models-0.1.1-alpha")

#NCBI translation table 5, minus UUG start
const startcodon = biore"(ATT)|(ATC)|(ATA)|(ATG)|(GTG)"d
const stopcodon = biore"(TAG)|(TAA)"d

function codonmatches(seq::CircularSequence, pattern)
    frames = [Int32[] for f in 1:3]
    for m in eachmatch(pattern, seq.sequence)
        i::Int32 = m.captured[1]
        push!(frames[mod1(i, 3)], i)
    end
    return frames
end

function getorfs!(writer::FASTA.Writer, id::String, genome::CircularSequence, strand::Char, minORF::Int)
    glength = length(genome)
    stops = codonmatches(genome, stopcodon)
    starts = codonmatches(genome, startcodon)
    for (f, frame) in enumerate(stops)
        for (s, stop) in enumerate(frame)
            stop > glength && break
            s + 1 > length(frame) && break
            nextstop = frame[s + 1]
            nextstop - stop - 1 < minORF && continue
            nextstartidx = searchsortedfirst(starts[f], stop)
            nextstartidx > length(starts[f]) && break
            nextstart = starts[f][nextstartidx]
            nextstop - nextstart + 1 < minORF && continue
            translation = translate(genome.sequence[nextstart:(nextstop-1)], code = ncbi_trans_table[5])
            translation[1] = AA_M
            write(writer, FASTA.Record(id * "_" * strand * "_" * string(nextstart) * "-" * string(nextstop), translation))
        end
    end
end

function hmmsearch(id::String, genome::CircularSequence, minORF::Int)
    writer = open(FASTA.Writer, "tmp.orfs.fa")
    getorfs!(writer, id, genome, '+', minORF)
    getorfs!(writer, id, reverse_complement(genome), '-', minORF)
    close(writer)
    hmmpath = joinpath(emmamodels, "hmms", "all.hmm")
    cmd = `hmmsearch --domtblout tmp.domt $hmmpath tmp.orfs.fa`
    outfile = "tmp.hmmsearch.out"
    run(pipeline(cmd, stdout=outfile))
    return "tmp.domt"
end

struct HMMmatch
    orf::String
    #a1::String
    tlen::Int
    query::String
    #a2::String
    qlen::Int
    Evalue::Float64
    score::Float64
    #bias::Float64
    #num::Int
    #of::Int
    #cEvalue::Float64
    #iEvalue::Float64
    #dscore::Float64
    #dbias::Float64
    hmm_from::Int
    hmm_to::Int
    ali_from::Int
    ali_to::Int
    env_from::Int
    env_to::Int
    #acc::Float64
    #description::String
end

function parse_domt(file::String)
    matches = HMMmatch[]
    open(file, "r") do infile
        while !eof(infile)
            line = readline(infile)
            startswith(line, "#") && continue
            bits = split(line, " ", keepempty=false)
            push!(matches, HMMmatch(bits[1], parse(Int, bits[3]), bits[4],
                parse(Int, bits[6]), parse(Float64, bits[7]), parse(Float64, bits[8]),
                parse(Int, bits[16]), parse(Int, bits[17]), parse(Int, bits[18]),
                parse(Int, bits[19]), parse(Int, bits[20]), parse(Int, bits[21])))
        end
    end
    return matches
end
