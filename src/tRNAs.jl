struct tRNA
    fm::FeatureMatch
    anticodon::Union{String,Missing}
    polyA::Int #if sequence completed by polyadenylation
end

const anticodon2trn = Dict("UGC" => "trnA", "GCA" => "trnC", "GUC" => "trnD", "UUC" => "trnE", "GAA" => "trnF", "UCC" => "trnG",
    "GUG" => "trnH", "GAU" => "trnI", "UUU" => "trnK", "UAA" => "trnL2", "UAG" => "trnL1", "CAU" => "trnM", "GUU" => "trnN",
    "UGG" => "trnP", "UUG" => "trnQ", "UCG" => "trnR", "GCU" => "trnS1", "UGA" => "trnS2", "UGU" => "trnT", "UAC" => "trnV",
    "UCA" => "trnW", "GUA" => "trnY")

#HGNC-approved symbols https://www.genenames.org/data/genegroup/#!/group/843
const trn2symbol = Dict("trnA" => "MT-TA", "trnC" => "MT-TC", "trnD" => "MT-TD", "trnE" => "MT-TE", "trnF" => "MT-TF", "trnG" => "MT-TG",
    "trnH" => "MT-TH", "trnI" => "MT-TI", "trnK" => "MT-TK", "trnL2" => "MT-TL2", "trnL1" => "MT-TL1", "trnM" => "MT-TM", "trnN" => "MT-TN",
    "trnP" => "MT-TP", "trnQ" => "MT-TQ", "trnR" => "MT-TR", "trnS1" => "MT-TS1", "trnS2" => "MT-TS2", "trnT" => "MT-TT", "trnV" => "MT-TV",
    "trnW" => "MT-TW", "trnY" => "MT-TY")

const trn2product = Dict("trnA" => "tRNA-Ala", "trnC" => "tRNA-Cys", "trnD" => "tRNA-Asp", "trnE" => "tRNA-Glu", "trnF" => "tRNA-Phe", "trnG" => "tRNA-Gly",
    "trnH" => "tRNA-His", "trnI" => "tRNA-Ile", "trnK" => "tRNA-Lys", "trnL2" => "tRNA-Leu", "trnL1" => "tRNA-Leu", "trnM" => "tRNA-Met", "trnN" => "tRNA-Asn",
    "trnP" => "tRNA-Pro", "trnQ" => "tRNA-Gln", "trnR" => "tRNA-Arg", "trnS1" => "tRNA-Ser", "trnS2" => "tRNA-Ser", "trnT" => "tRNA-Thr", "trnV" => "tRNA-Val",
    "trnW" => "tRNA-Trp", "trnY" => "tRNA-Tyr")

#generic vertebrate models
#= const model2trn = Dict("A.seed25-1"=>"trnA","C.seed25-1"=>"trnC","D.seed25-1"=>"trnD","E.seed25-1"=>"trnE",
    "F.seed25-1"=>"trnF","G.seed25-1"=>"trnG","H.seed25-1"=>"trnH","I.seed25-1"=>"trnI","K.seed25-1"=>"trnK",
    "L_infernalcluster1.seed25-1"=>"trnL1","L_infernalcluster2.seed25-1"=>"trnL2","M.seed25-1"=>"trnM",
    "N.seed25-1"=>"trnN","P.seed25-1"=>"trnP","Q.seed25-1"=>"trnQ","R.seed25-1"=>"trnR","S1.seed25-1"=>"trnS1",
    "S2.seed25-1"=>"trnS2","T.seed25-1"=>"trnT","V.seed25-1"=>"trnV","W.seed25-1"=>"trnW","Y.seed25-1"=>"trnY") =#

#= const trn2anticodonpos = Dict("A.seed25-1"=>31,"C.seed25-1"=>31,"D.seed25-1"=>32,"E.seed25-1"=>31,"F.seed25-1"=>32,"G.seed25-1"=>31,
    "H.seed25-1"=>33,"I.seed25-1"=>30,"K.seed25-1"=>32,"L_infernalcluster1.seed25-1"=>34,"L_infernalcluster2.seed25-1"=>31,"M.seed25-1"=>32,"N.seed25-1"=>32,
    "P.seed25-1"=>31,"Q.seed25-1"=>31,"R.seed25-1"=>32,"S1.seed25-1"=>27,"S2.seed25-1"=>31,"T.seed25-1"=>32,"V.seed25-1"=>32,
    "W.seed25-1"=>32,"Y.seed25-1"=>32) =#

#Pavels' fish-specific models
const model2trn = Dict("trnA_fish" => "trnA", "trnC_fish" => "trnC", "trnD_fish" => "trnD", "trnE_fish" => "trnE",
    "trnF_fish" => "trnF", "trnG_fish" => "trnG", "trnH_fish" => "trnH", "trnI_fish" => "trnI", "trnK_fish" => "trnK",
    "trnL1_fish" => "trnL1", "trnL2_fish" => "trnL2", "trnM_fish" => "trnM",
    "trnN_fish" => "trnN", "trnP_fish" => "trnP", "trnQ_fish" => "trnQ", "trnR_fish" => "trnR", "trnS1_fish" => "trnS1",
    "trnS2_fish" => "trnS2", "trnT_fish" => "trnT", "trnV_fish" => "trnV", "trnW_fish" => "trnW", "trnY_fish" => "trnY")

const model2anticodonpos = deserialize(joinpath(emmamodels, "trn", "model2anticodonpos.dict"))

const trn2model = Dict(value => key for (key, value) in model2trn)

const trn2anticodon = Dict(value => key for (key, value) in anticodon2trn)

function get_best_trns(alltrns::Vector{tRNA}, glength::Integer)
    besttrns = tRNA[]
    isempty(alltrns) && return besttrns
    sort!(alltrns, by=t -> t.fm.target_from)
    current_trn = first(alltrns)
    for t in alltrns[2:end]
        if current_trn.fm.target_from + current_trn.fm.target_length >= t.fm.target_from + 5
            if t.fm.evalue < current_trn.fm.evalue
                current_trn = t
            end
        else
            push!(besttrns, current_trn)
            current_trn = t
        end
    end
    push!(besttrns, current_trn)
    if last(besttrns).fm.target_from + last(besttrns).fm.target_length - glength >= first(besttrns).fm.target_from + 5
        deleteat!(besttrns, last(besttrns).fm.evalue < first(besttrns).fm.evalue ? 1 : length(besttrns))
    end
    return besttrns
end

function anticodon(qfrom::Int, qseq::AbstractString, tseq::AbstractString, pos::Int)
    @assert length(qseq) == length(tseq)
    pos < qfrom && return missing
    pointer = qfrom - 1
    trunc = false
    for (i, c) in enumerate(qseq)
        c == '.' && continue
        c == '*' && continue
        c == ' ' && continue
        isdigit(c) && continue
        if c == ']'
            trunc = false
            continue
        end
        if trunc == true
            continue
        elseif c == '['
            pointer += parse(Int, qseq[i+1:(findnext("]", qseq, i + 2)[1]-1)])
        elseif c == '<'
            trunc = true
            continue
        else
            pointer += 1
        end
        pos < pointer && return missing
        pointer == pos && return tseq[i:i+2]
    end
    return missing
end

function rationalise_trn_alignments(trns::Vector{tRNA})
    #filter out duplicate hits
    for i in trns, j in trns
        i == j && continue
        if i.fm.query == j.fm.query
            # if they are the same gene delete the poorer match
            todelete = i.fm.evalue > j.fm.evalue ? i : j
            deleteat!(trns, findfirst(x -> x == todelete, trns))
            return (rationalise_trn_alignments(trns))
        end
    end
    trns
end

function parse_trn_alignments(file::String, glength::Integer)
    trns = tRNA[]
    open(file, "r") do infile
        while !eof(infile)
            line = readline(infile)
            !startswith(line, ">>") && continue
            readline(infile) #header line
            readline(infile) #dashes
            bits = split(readline(infile), " ", keepempty=false)
            readline(infile) #blank
            readline(infile) #NC
            readline(infile) #CS
            qseqline = strip(readline(infile))
            query = qseqline[1:(findfirst(" ", qseqline)[1]-1)]
            qseqline = lstrip(qseqline[length(query)+1:end])
            qseq = qseqline[(findfirst(" ", qseqline)[1]+1):(findlast(" ", qseqline)[1]-1)]
            readline(infile) #matches
            tseqline = strip(readline(infile))
            target = tseqline[1:(findfirst(" ", tseqline)[1]-1)]
            tseqline = lstrip(tseqline[length(target)+1:end])
            tseq = tseqline[(findfirst(" ", tseqline)[1]+1):(findlast(" ", tseqline)[1]-1)]
            tstrand = bits[12][1]
            target_from = parse(Int, bits[10])
            tto = parse(Int, bits[11])
            if tstrand == '-'
                target_from = reverse_complement(target_from, glength)
                tto = reverse_complement(tto, glength)
            end
            #find anticodon position
            qfrom = parse(Int, bits[7])
            trn = model2trn[query]
            expected = trn2anticodon[trn]
            anticod = anticodon(qfrom, qseq, tseq, model2anticodonpos[query])
            if haskey(anticodon2trn, anticod) == true
                push!(trns, tRNA(FeatureMatch(target, anticodon2trn[anticod], tstrand, qfrom, parse(Int, bits[8]), target_from, tto - target_from + 1, parse(Float64, bits[3])), anticod, false))
            else
                push!(trns, tRNA(FeatureMatch(target, trn, tstrand, qfrom, parse(Int, bits[8]), target_from, tto - target_from + 1, parse(Float64, bits[3])), anticod, false))
            end
        end
    end
    filter!(x -> x.fm.evalue < 1e-5, trns)
end

function cmsearch(tempfile::TempFile, modeldir::String, modelfile::String)
    cmpath = joinpath(emmamodels, modeldir, modelfile)
    name = tempfilename(tempfile, "extended.fa")
    cmd = `cmsearch $cmpath $name`
    outfile = tempfilename(tempfile, "cmsearch.out")
    run(pipeline(cmd, stdout=outfile))
    return outfile
end

function cmsearch(tempfile::TempFile, id::String, modeldir::String, tRNA::LongDNA{4})
    name = tempfilename(tempfile, "tmp.fa")
    open(FASTA.Writer, name) do writer
        write(writer, FASTA.Record(id, tRNA))
    end
    cmpath = joinpath(emmamodels, modeldir, trn2model[id] * ".cm")
    cmd = `cmsearch $cmpath $name`
    outfile = tempfilename(tempfile, "cmsearch.out")
    run(pipeline(cmd, stdout=outfile))
    return outfile
end

function get_overlapped_trns(trns::Vector{tRNA}, glength::Integer)
    #assumes trns is sorted by feature start
    overlapped = Vector{Tuple{tRNA,Int}}(undef, 0)
    isempty(trns) && return overlapped
    for i in 1:length(trns)-1
        if circularoverlap(trns[i].fm, trns[i+1].fm, glength) > 0
            push!(overlapped, (trns[i], mod1(trns[i+1].fm.target_from - 1, glength)))
        end
    end
    if circularoverlap(last(trns).fm, first(trns).fm, glength) > 0
        push!(overlapped, (last(trns), glength + first(trns).fm.target_from - 1))
    end
    return overlapped
end
