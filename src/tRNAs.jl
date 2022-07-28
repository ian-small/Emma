struct CMAlignment_trn
    query::String
    target::String
    Evalue::Float64
    #score::Float64
    #bias::Float64
    #mdl::String
    qfrom::Int
    qto::Int
    tfrom::Int
    tto::Int
    tstrand::Char
    #acc::Float64
    #trunc::Bool
    #gc::Float64
    #qseq::String
    tseq::String
    anticodon::Union{String,Missing}
    polyA::Int #> 0 if sequence completed by polyadenylation
end

const anticodon2trn = Dict("UGC"=>"trnA","GCA"=>"trnC","GUC"=>"trnD","UUC"=>"trnE","GAA"=>"trnF","UCC"=>"trnG",
    "GUG"=>"trnH","GAU"=>"trnI","UUU"=>"trnK","UAA"=>"trnL","UAG"=>"trnL","CAU"=>"trnM","GUU"=>"trnN",
    "UGG"=>"trnP","UUG"=>"trnQ","UCG"=>"trnR","UCU"=>"trnS1","UGA"=>"trnS2","UGU"=>"trnT","UAC"=>"trnV",
    "UCA"=>"trnW","GUA"=>"trnY")

const trn2anticodonpos = Dict("trnA"=>31,"trnC"=>31,"trnD"=>32,"trnE"=>31,"trnF"=>32,"trnG"=>31,
    "trnH"=>33,"trnI"=>30,"trnK"=>32,"trnL"=>31,"trnM"=>32,"trnN"=>32,
    "trnP"=>31,"trnQ"=>31,"trnR"=>27,"trnS1"=>19,"trnS2"=>31,"trnT"=>32,"trnV"=>29,
    "trnW"=>32,"trnY"=>32)

const model2trn = Dict("A.seed25-1"=>"trnA","C.seed25-1"=>"trnC","D.seed25-1"=>"trnD","E.seed25-1"=>"trnE",
    "F.seed25-1"=>"trnF","G.seed25-1"=>"trnG","H.seed25-1"=>"trnH","I.seed25-1"=>"trnI","K.seed25-1"=>"trnK",
    "L_infernalcluster2.seed25-1"=>"trnL","M.seed25-1"=>"trnM",
    "N.seed25-1"=>"trnN","P.seed25-1"=>"trnP","Q.seed25-1"=>"trnQ","trnR"=>"trnR","trnS1"=>"trnS1",
    "S2.seed25-1"=>"trnS2","T.seed25-1"=>"trnT","trnV"=>"trnV","W.seed25-1"=>"trnW","Y.seed25-1"=>"trnY")

const trn2model = Dict(value => key for (key, value) in model2trn)

function anticodon(qfrom::Int, qseq::AbstractString, tseq::AbstractString, pos::Int)
    @assert length(qseq) == length(tseq)
    pos < qfrom && return missing
    pointer = qfrom - 1
    for (i,c) in enumerate(qseq)
        c == '.' && continue
        c == '*' && continue
        isdigit(c) && continue
        c == ']' && continue
        if c == '['
            pointer += parse(Int, qseq[i+1:(findnext("]",qseq,i+2)[1]-1)])
        else
            pointer += 1
        end
        pos < pointer && return missing
        pointer == pos && return tseq[i:i+2]
    end
    return missing
end

function parse_trn_alignments(file::String, glength::Integer)
    alignments = CMAlignment_trn[]
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
            query = qseqline[1:(findfirst(" ",qseqline)[1]-1)]
            qseqline = lstrip(qseqline[length(query)+1:end])
            qseq = qseqline[(findfirst(" ",qseqline)[1]+1):(findlast(" ",qseqline)[1]-1)]
            readline(infile) #matches
            tseqline = strip(readline(infile))
            target = tseqline[1:(findfirst(" ",tseqline)[1]-1)]
            tseqline = lstrip(tseqline[length(target)+1:end])
            tseq = tseqline[(findfirst(" ",tseqline)[1]+1):(findlast(" ",tseqline)[1]-1)]
            tstrand = bits[12][1]
            tfrom = parse(Int, bits[10])
            tto = parse(Int, bits[11])
            if tstrand == '-'
                tfrom = reverse_complement(tfrom, glength)
                tto = reverse_complement(tto, glength)
            end
            #find anticodon position
            trn = model2trn[query]
            anticodonpos = trn2anticodonpos[trn]
            qfrom = parse(Int, bits[7])
            push!(alignments, CMAlignment_trn(trn, target, parse(Float64, bits[3]),
                qfrom, parse(Int, bits[8]), tfrom, tto, tstrand,
                tseq, anticodon(qfrom, qseq, tseq, anticodonpos), 0))
        end
    end
    return alignments
end

function cmsearch(modeldir::String, modelfile::String)
    cmpath = joinpath(emmamodels, modeldir, modelfile)
    cmd = `cmsearch $cmpath tmp.extended.fa`
    outfile = "tmp.cmsearch.out"
    run(pipeline(cmd, stdout=outfile))
    return outfile
end

function cmsearch(id::String, modeldir::String, tRNA::LongDNA{2})
    writer = open(FASTA.Writer, "tmp.fa")
    write(writer, FASTA.Record(id, tRNA))
    close(writer)
    cmpath = joinpath(emmamodels, modeldir, trn2model[id] * ".cm")
    cmd = `cmsearch $cmpath tmp.fa`
    outfile = "tmp.cmsearch.out"
    run(pipeline(cmd, stdout=outfile))
    return outfile
end
    
