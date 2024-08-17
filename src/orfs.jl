#from corresponding NCBI translation tables
const ncbi_start_codons = Dict(2 => biore"(ATT)|(ATC)|(ATA)|(ATG)|(GTG)|(TTG)|(CTG)"d, 5 => biore"(ATT)|(ATC)|(ATA)|(ATG)|(GTG)|(TTG)"d)
const ncbi_stop_codons = Dict(2 => biore"(TAG)|(TAA)|(AGA)|(AGG)"d, 5 => biore"(TAG)|(TAA)"d)

const START_CODON_PENALTY_WEIGHTING = 30.0
const INTRUSION_PENALTY_WEIGHTING = 100.0
const EXPANSION_PENALTY_WEIGHTING = 20.0
const SCANNING_DISTANCE_THRESHOLD = 30
const SCANNING_PENALTY_WEIGHTING = 50.0

Codon = LongSubSeq{DNAAlphabet{4}}

#HGNC-approved symbols https://www.genenames.org/data/genegroup/#!/group/1974
const cds2symbol = Dict("ATP6" => "MT-ATP6", "ATP8" => "MT-ATP8", "COX1" => "MT-CO1", "COX2" => "MT-CO2", "COX3" => "MT-CO3", "CYTB" => "MT-CYTB",
    "ND1" => "MT-ND1", "ND2" => "MT-ND2", "ND3" => "MT-ND3", "ND4" => "MT-ND4", "ND4L" => "MT-ND4L", "ND5" => "MT-ND5", "ND6" => "MT-ND6")

const cds2product = Dict("ND1" => "NADH dehydrogenase subunit 1", "nad1" => "NADH dehydrogenase subunit 1", "ND2" => "NADH dehydrogenase subunit 2", "nad2" => "NADH dehydrogenase subunit 2",
    "ND3" => "NADH dehydrogenase subunit 3", "nad3" => "NADH dehydrogenase subunit 3", "ND4" => "NADH dehydrogenase subunit 4", "nad4" => "NADH dehydrogenase subunit 4", "ND4L" => "NADH dehydrogenase subunit 4L",
    "nad4L" => "NADH dehydrogenase subunit 4L", "ND5" => "NADH dehydrogenase subunit 5", "nad5" => "NADH dehydrogenase subunit 5", "ND6" => "NADH dehydrogenase subunit 6", "nad6" => "NADH dehydrogenase subunit 6",
    "COX1" => "cytochrome c oxidase subunit I", "cox1" => "cytochrome c oxidase subunit I", "COX2" => "cytochrome c oxidase subunit II", "cox2" => "cytochrome c oxidase subunit II",
    "COX3" => "cytochrome c oxidase subunit III", "cox3" => "cytochrome c oxidase subunit III", "ATP6" => "ATP synthase F0 subunit 6", "atp6" => "ATP synthase F0 subunit 6", "ATP8" => "ATP synthase F0 subunit 8",
    "atp8" => "ATP synthase F0 subunit 8", "CYTB" => "cytochrome b", "cob" => "cytochrome b")


function codonmatches(seq::CircularSequence, pattern)::Vector{Vector{Int32}}
    frames = [Int32[] for f in 1:3]
    for m in eachmatch(pattern, seq.sequence[1:seq.length+2])
        n_certain(matched(m)) < 3 && continue
        i::Int32 = m.captured[1]
        push!(frames[mod1(i, 3)], i)
    end
    return frames
end

function getcodons(seq::CircularSequence, pattern)
    positions = [Int32[] for f in 1:3]
    codons = Dict{Int32,Codon}()
    for m in eachmatch(pattern, seq.sequence[1:seq.length+2])
        n_certain(matched(m)) < 3 && continue
        i::Int32 = m.captured[1]
        push!(positions[mod1(i, 3)], i)
        codons[i] = matched(m)
    end
    return positions, codons
end

function getorfs!(writer::FASTA.Writer, id::AbstractString, genome::CircularSequence, translation_table::Int, strand::Char, starts::Vector{Vector{Int32}}, stops::Vector{Vector{Int32}}, minORF::Int)
    glength = length(genome)
    norfs = 0
    for (f, frame) in enumerate(starts)
        nextstop = 0
        for (s, start) in enumerate(frame)
            start < nextstop && continue
            nextstopidx = searchsortedfirst(stops[f], start + 1)
            nextstop = first(stops[mod1(f - mod1(glength, 3), 3)]) #frame of next stop when wrapping depends on genome length
            if nextstopidx <= length(stops[f])
                nextstop = stops[f][nextstopidx]
            end
            circulardistance(start, nextstop, glength) < minORF && continue
            if nextstop < start
                nextstop += glength
            end
            translation = BioSequences.translate(genome.sequence[start:(nextstop-1)], code=ncbi_trans_table[translation_table])
            translation[1] = AA_M
            write(writer, FASTA.Record(id * "*" * strand * "*" * string(start) * "-" * string(nextstop), translation))
            norfs += 1
        end
    end
    norfs
end

function orfsearch(tempfile::TempFile, id::AbstractString, genome::CircularSequence, translation_table::Int,
    fstarts::Vector{Vector{Int32}}, fstops::Vector{Vector{Int32}},
    rstarts::Vector{Vector{Int32}}, rstops::Vector{Vector{Int32}}, minORF::Int)
    out = tempfilename(tempfile, "orfs.fa")
    norfs = open(FASTA.Writer, out) do writer
        norfs = getorfs!(writer, id, genome, translation_table, '+', fstarts, fstops, minORF)
        norfs += getorfs!(writer, id, reverse_complement(genome), translation_table, '-', rstarts, rstops, minORF)
        norfs
    end
    ret = tempfilename(tempfile, "tmp.domt")
    if norfs == 0
        error("no orfs found!")
    end
    hmmpath = joinpath(emmamodels, "cds", "all_cds.hmm")
    hmmsearch = which("hmmsearch")
    cmd = `$hmmsearch --domtblout $ret $hmmpath $out`
    outfile = tempfilename(tempfile, "hmmsearch.out")
    run(pipeline(Cmd(cmd, windows_hide=true), stdout=outfile))
    return ret
end

struct HMMmatch
    orf::String
    strand::Char
    #a1::String
    #tlen::Int32
    query::String
    #a2::String
    qlen::Int
    evalue::Float64
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
    matches = FeatureMatch[]
    open(file, "r") do infile
        while !eof(infile)
            line = readline(infile)
            startswith(line, "#") && continue
            #convert from ORF coordinates to genome coordinates (5'-3')
            bits = split(line, " ", keepempty=false)
            evalue = parse(Float64, bits[13]) # i-Evalue (indepent Evalue)
            evalue > 1 && continue # filter out poor matches
            orf = bits[1]
            orfbits = split(orf, "*")
            ends = split(orfbits[3], "-")
            strand = orfbits[2][1]
            orfstart = parse(Int32, ends[1])
            orfalifrom = parse(Int32, bits[18])
            orfalito = parse(Int32, bits[19])
            ali_from = orfstart + 3 * (orfalifrom - 1)
            ali_length = 3 * (orfalito - orfalifrom + 1)
            if ali_from > glength
                ali_from = mod1(ali_from, glength)
            end

            # note that model coordinates are converted to nucleotide coordinates
            push!(matches, FeatureMatch(orf, bits[4], strand, 3 * parse(Int, bits[16]) - 2, 3 * parse(Int, bits[17]), ali_from, ali_length, evalue))
        end
    end
    @debug matches
    rationalise_matches!(matches, glength)
end

const start_codon_penalty = Dict(LongSequence{DNAAlphabet{4}}("ATT") => -9.477758266443889, LongSequence{DNAAlphabet{4}}("TTG") => -7.892795765722732,
    LongSequence{DNAAlphabet{4}}("ATA") => -8.155830171556525, LongSequence{DNAAlphabet{4}}("CTG") => -7.155830171556525, LongSequence{DNAAlphabet{4}}("GTG") => -3.3692338096657193,
    LongSequence{DNAAlphabet{4}}("ATC") => -10.0, LongSequence{DNAAlphabet{4}}("ATG") => 0.0)

function fix_start_and_stop_codons!(hmm_matches, trns, starts, startcodons, stops, glength)

    function is_possible_start(start, model_start, leftwindow, rightwindow, glength)
        d = closestdistance(start, model_start, glength)
        #must be in frame
        mod(d, 3) ≠ 0 && return false
        if d > glength / 2
            d -= glength
        end
        d > leftwindow && return false
        d < -rightwindow && return false
        return true
    end

    sort!(hmm_matches; by=x -> x.target_from)
    for (i, hmm_match) in enumerate(hmm_matches)
        #println(hmm_match)
        @debug hmm_match
        hmmstart = hmm_match.target_from

        upstream_cds = hmm_matches[mod1(i - 1, length(hmm_matches))]
        @debug "upstream_cds: $(upstream_cds.query)"

        #this works because if no trn gene is upstream of the CDS start, this will select the last trn gene in the genome, which is upstream of the first CDS
        trnidx = searchsortedfirst(trns, hmmstart, lt=(t, x) -> t.fm.target_from < x)
        upstream_tRNA = trns[mod1(trnidx - 1, length(trns))]
        distance_to_upstream_tRNA = circulardistance(upstream_tRNA.fm.target_from + upstream_tRNA.fm.target_length - 1, hmmstart, glength)
        @debug "upstream_tRNA: $(upstream_tRNA.fm.query)"

        inframe_stops = Int32[]
        for s in stops
            append!(inframe_stops, filter(x -> mod(circulardistance(x, hmmstart, glength), 3) == 0, s))
        end
        sort!(inframe_stops)
        stop_idx = searchsortedfirst(inframe_stops, hmmstart)
        upstream_stop = inframe_stops[mod1(stop_idx - 1, length(inframe_stops))]
        distance_to_upstream_stop = circulardistance(upstream_stop, hmmstart, glength)
        @debug "$hmmstart, $distance_to_upstream_stop"

        possible_starts = Int32[]
        for frame in starts, ps in frame
            if is_possible_start(ps, hmmstart, min(distance_to_upstream_tRNA, distance_to_upstream_stop), 50, glength)
                push!(possible_starts, ps)
            end
        end
        if isempty(possible_starts)
            @warn "no starts for $hmm_match"
            continue
        end
        @debug possible_starts

        #apply penalties
        penalties = zeros(Float64, length(possible_starts))
        # start codon penalty
        for (i, ps) in enumerate(possible_starts)
            penalties[i] += START_CODON_PENALTY_WEIGHTING * start_codon_penalty[startcodons[ps]]
        end
        # scanning penalty
        if distance_to_upstream_tRNA < SCANNING_DISTANCE_THRESHOLD
            idxs = sortperm(possible_starts, by=x -> closestdistance(upstream_tRNA.fm.target_from + upstream_tRNA.fm.target_length - 1, x, glength))
            for (i, idx) in enumerate(idxs)
                penalties[idx] -= SCANNING_PENALTY_WEIGHTING * (i - 1)
            end
        end
        # intrusion/expansion penalty
        for (i, ps) in enumerate(possible_starts)
            relative_to_hmm = closestdistance(hmmstart, ps, glength)
            if relative_to_hmm > 0
                penalties[i] -= INTRUSION_PENALTY_WEIGHTING * relative_to_hmm
            elseif relative_to_hmm < 0
                penalties[i] += EXPANSION_PENALTY_WEIGHTING * relative_to_hmm
            end
        end
        @debug penalties
        beststart = possible_starts[argmax(penalties)]
        @debug "best start: $beststart)"
        @info "$(hmm_match.query)"
        #pick first stop, or if none before tRNA, see if can construct stop by polyadenylation
        frame = mod1(beststart, 3)
        stop_idx = searchsortedfirst(stops[frame], beststart) #index of first in-frame stop following best start codon
        if stop_idx == length(stops[frame]) + 1 #no stop before end of genome
            frame = 3 - mod1(glength - beststart, 3)
            stop_idx = searchsortedfirst(stops[frame], 1) #restart search from start of genome in the correct frame
        end
        next_stop = stops[frame][stop_idx]
        @debug "first stop: $next_stop"
        #modify HMM match
        cds = FeatureMatch(hmm_match.id, hmm_match.query, hmm_match.strand,
            hmm_match.model_from, hmm_match.model_to, beststart, circulardistance(beststart, next_stop, glength), hmm_match.evalue)
        #println(cds)
        replace!(hmm_matches, hmm_match => cds)
    end
    return hmm_matches
end

