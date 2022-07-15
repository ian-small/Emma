struct CMAlignment_rrn
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
    #tseq::String
end

const model2rrn = Dict("endrrnl"=>"rrnL","endrrns"=>"rrnS","1-167"=>"rrnS","1-70"=>"rrnL")

function parse_rrn_alignments(file::String, glength::Integer)
    alignments = CMAlignment_rrn[]
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
            tstrand = bits[12][1]
            tfrom = parse(Int, bits[10])
            tto = parse(Int, bits[11])
            if tstrand == '-'
                tfrom = reverse_complement(tfrom, glength)
                tto = reverse_complement(tto, glength)
            end
            rrn = model2rrn[query]
            qfrom = parse(Int, bits[7])
            push!(alignments, CMAlignment_rrn(rrn, target, parse(Float64, bits[3]),
                qfrom, parse(Int, bits[8]), tfrom, tto, tstrand))
        end
    end
    return alignments
end

function rrnsearch()
    hmmpath = joinpath(emmamodels, "rrn", "all_rrn.hmm")
    cmd = `nhmmer --tblout tmp.tbl $hmmpath tmp.extended.fa`
    outfile = "tmp.nhmmer.out"
    run(pipeline(cmd, stdout=outfile))
    return "tmp.tbl"
end

struct nHMMmatch
    genome::String
    #a1::String
    query::String
    #a2::String
    hmm_from::Int
    hmm_to::Int
    ali_from::Int
    ali_to::Int
    env_from::Int
    env_to::Int
    #seqlen::Int
    strand::Char
    Evalue::Float64
    score::Float64
    #bias::Float64
    #description::String
end

function parse_tbl(file::String)
    matches = nHMMmatch[]
    open(file, "r") do infile
        while !eof(infile)
            line = readline(infile)
            startswith(line, "#") && continue
            bits = split(line, " ", keepempty=false)
            push!(matches, nHMMmatch(bits[1], model2rrn[bits[3]],
                parse(Int, bits[5]), parse(Float64, bits[6]), parse(Int, bits[7]),
                parse(Int, bits[8]), parse(Int, bits[9]), parse(Int, bits[10]),
                bits[12][1], parse(Float64, bits[13]), parse(Float64, bits[14])))
        end
    end
    return matches
end

struct rRNA
    start::Vector{nHMMmatch}
    stop::Vector{CMAlignment_rrn}
end