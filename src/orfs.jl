using Artifacts
using BioSequences
using FASTX

const emmamodels = joinpath(artifact"Emma-models", "emma-models-0.1.5-alpha")

#NCBI translation table 5, minus UUG start
const startcodon = biore"(ATT)|(ATC)|(ATA)|(ATG)|(GTG)|(TTG)"d
const stopcodon = biore"(TAG)|(TAA)"d

Codon = LongSubSeq{DNAAlphabet{2}}

function codonmatches(seq::CircularSequence, pattern)::Vector{Vector{Int32}}
    frames = [Int32[] for f in 1:3]
    for m in eachmatch(pattern, seq.sequence)
        i::Int32 = m.captured[1]
        push!(frames[mod1(i, 3)], i)
    end
    return frames
end

function getcodons(seq::CircularSequence, pattern)
    positions = [Int32[] for f in 1:3]
    codons = [Codon[] for f in 1:3]
    for m in eachmatch(pattern, seq.sequence)
        i::Int32 = m.captured[1]
        push!(positions[mod1(i, 3)], i)
        push!(codons[mod1(i, 3)], matched(m))
    end
    return positions, codons
end

function getorfs!(writer::FASTA.Writer, id::AbstractString, genome::CircularSequence, strand::Char, starts::Vector{Vector{Int32}}, stops::Vector{Vector{Int32}}, minORF::Int)
    glength = length(genome)
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
            translation = BioSequences.translate(genome.sequence[nextstart:(nextstop-1)], code = ncbi_trans_table[5])
            translation[1] = AA_M
            write(writer, FASTA.Record(id * "*" * strand * "*" * string(nextstart) * "-" * string(nextstop), translation))
        end
    end
end

function orfsearch(id::AbstractString, genome::CircularSequence, fstarts::Vector{Vector{Int32}}, fstops::Vector{Vector{Int32}},
     rstarts::Vector{Vector{Int32}}, rstops::Vector{Vector{Int32}}, minORF::Int)
    writer = open(FASTA.Writer, "tmp.orfs.fa")
    getorfs!(writer, id, genome, '+', fstarts, fstops, minORF)
    getorfs!(writer, id, reverse_complement(genome), '-', rstarts, rstops, minORF)
    close(writer)
    hmmpath = joinpath(emmamodels, "cds", "all_cds.hmm")
    cmd = `hmmsearch --domtblout tmp.domt $hmmpath tmp.orfs.fa`
    outfile = "tmp.hmmsearch.out"
    run(pipeline(cmd, stdout=outfile))
    return "tmp.domt"
end

struct HMMmatch
    orf::String
    strand::Char
    #a1::String
    #tlen::Int32
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
    hmm_from::Int32
    hmm_to::Int32
    ali_from::Int32
    ali_to::Int32
    #env_from::Int
    #env_to::Int
    #acc::Float64
    #description::String
end

function parse_domt(file::String, glength::Integer)
    matches = HMMmatch[]
    open(file, "r") do infile
        while !eof(infile)
            line = readline(infile)
            startswith(line, "#") && continue
            #convert from ORF coordinates to genome coordinates (5'-3')
            bits = split(line, " ", keepempty=false)
            orf = bits[1]
            orfbits = split(orf, "*")
            ends = split(orfbits[3], "-")
            strand = orfbits[2][1]
            orfstart = parse(Int32, ends[1])
            orfalifrom = parse(Int32, bits[18])
            orfalito = parse(Int32, bits[19])
            ali_from::Int32 = orfstart + 3*(orfalifrom-1)
            ali_to::Int32 = orfstart + 3*(orfalito - orfalifrom + 1) - 1
            if ali_from > glength
                ali_from = mod1(ali_from, glength)
                ali_to = mod1(ali_to, glength)
            end
            push!(matches, HMMmatch(orf, strand, bits[4],
                parse(Int32, bits[6]), parse(Float64, bits[7]), parse(Float64, bits[8]),
                parse(Int32, bits[16]), parse(Int32, bits[17]), ali_from, ali_to))
        end
    end
    return matches
end

function merge_hmm_matches(match1::HMMmatch, match2::HMMmatch)
    @assert match1.orf == match2.orf && match1.query == match2.query
    return HMMmatch(match1.orf, match1.strand, match1.query, match1.qlen, match1.Evalue, match1.score,
    min(match1.hmm_from, match2.hmm_from), max(match1.hmm_to, match2.hmm_to), min(match1.ali_from, match2.ali_from),
    max(match1.ali_to, match2.ali_to))
end

const target_encoding = Dict(LongSequence{DNAAlphabet{2}}("ATT") => 0.0195429, LongSequence{DNAAlphabet{2}}("TTG") => 0.0280438, LongSequence{DNAAlphabet{2}}("ATG") => 0.180693,
                            LongSequence{DNAAlphabet{2}}("GTG") => 0.0286169, LongSequence{DNAAlphabet{2}}("ATA") => 0.024437, LongSequence{DNAAlphabet{2}}("ATC") => 0.0189099)

function fix_start_and_stop_codons!(hmm_matches, trns, starts, startcodons, stops, startcodon_model, glength)
    for (i,hmm_match) in enumerate(hmm_matches)
        println(hmm_match)
        hmmstart = hmm_match.ali_from
        inframe_starts = starts[mod1(hmmstart, 3)]
        inframe_stops = stops[mod1(hmmstart, 3)]
        upstream_cds = hmm_matches[mod1(i-1, length(hmm_matches))]
        println("upstream_cds: ", upstream_cds.query)
        trnidx = searchsortedfirst(trns, hmmstart, lt=(t,x)->t.tfrom < x)
        upstream_tRNA = trns[mod1(trnidx-1, length(trns))]
        println("upstream_tRNA: ", upstream_tRNA.query)
        stop_idx = searchsortedfirst(inframe_stops, hmmstart)
        upstream_stop = inframe_stops[stop_idx] == hmmstart ? inframe_stops[stop_idx] : inframe_stops[mod1(stop_idx-1, length(inframe_stops))]
        #select relevant starts (same strand and frame, not before in-frame stop, not before tRNA end, not beyond -100 or +50 wrt to hmm start)
        possible_starts = filter(x -> hmmstart-100 <= x <= hmmstart+50, inframe_starts) #can miss starts after hmm start and after end of genome...
        filter!(x -> x > upstream_stop, possible_starts)
        filter!(x -> x > upstream_tRNA.tto, possible_starts)
        if isempty(possible_starts)
            println("no starts for ", hmm_match)
            continue
        end
        model_inputs = zeros(Float64, length(possible_starts), 5)
        for (i,ps) in enumerate(possible_starts)
            #calculate model inputs :target_encoding, :relative_to_hmm, :phase_to_hmm,:relative_to_upstream_gene, :relative_to_upstream_stop
            #println(ps, ' ', findfirst(ps .== inframe_starts), ' ', startcodons[mod1(hmmstart, 3)][findfirst(ps .== inframe_starts)])
            #, ' ', target_encoding[startcodons[findfirst(ps .== inframe_starts)]])
            model_inputs[i, 1] = target_encoding[startcodons[mod1(hmmstart, 3)][findfirst(ps .== inframe_starts)]]
            relative_to_hmm = circulardistance(ps, hmmstart, glength)
            if relative_to_hmm > glength/2; relative_to_hmm -= glength; end
            model_inputs[i, 2] = relative_to_hmm
            model_inputs[i, 3] = mod(relative_to_hmm, 3)
            upstream_gene_end = max(upstream_cds.ali_to, upstream_tRNA.tto)
            relative_to_upstream_gene = circulardistance(upstream_gene_end, ps, glength)
            if relative_to_upstream_gene > glength/2; relative_to_upstream_gene -= glength; end
            model_inputs[i, 4] = relative_to_upstream_gene
            relative_to_upstream_stop = circulardistance(upstream_stop, ps, glength)
            if relative_to_upstream_stop > glength/2; relative_to_upstream_stop -= glength; end
            model_inputs[i, 5] = relative_to_upstream_stop
            #println(model_inputs[i,:])
        end
        #predict with XGBoost model
        scores = XGBoost.predict(startcodon_model, model_inputs)
        #println(scores)
        #pick top scoring start codon
        maxscore, maxidx = findmax(scores)
        beststart = possible_starts[maxidx]
        println("best start: ", beststart, ' ', maxscore)
        #pick first stop, or if none before tRNA, see if can construct stop by polyadenylation
        stop_idx = searchsortedfirst(stops[mod1(beststart, 3)], beststart) #index of first in-frame stop following best start codon
        next_stop = stops[mod1(beststart, 3)][mod1(stop_idx, length(stops[mod1(beststart, 3)]))]
        trn_idx = searchsortedfirst(trns, beststart, lt=(t,x)->t.tfrom < x)
        next_trn = trns[mod1(trn_idx, length(trns))]
        #next_stop points to first nucleotide of stop codon
        if circulardistance(beststart, next_trn.tfrom-1, glength) < circulardistance(beststart, next_stop, glength)
            next_stop = next_trn.tfrom-1 #we can't actually check here if next_trn-1 can form a valid stop by polyadenylation; need to verify subsequently
        end
        println("first stop: ", next_stop)
        #modify HMM match
        cds =  HMMmatch(hmm_match.orf, hmm_match.strand, hmm_match.query, hmm_match.qlen, hmm_match.Evalue, hmm_match.score,
                        hmm_match.hmm_from, hmm_match.hmm_to, beststart, next_stop)
        replace!(hmm_matches, hmm_match=>cds)
    end
    return hmm_matches
end
