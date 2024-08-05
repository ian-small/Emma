const model2rrn = Dict("16srna"=>"rrnL","12srna"=>"rrnS")

const rrn2symbol = Dict("rrnS"=>"MT-RNR1", "rrnL"=>"MT-RNR2")

function rrnsearch(uid::UUID)
    hmmpath = joinpath(emmamodels, "rrn", "all_rrn.hmm")
    cmd = `nhmmer --tblout $uid.tbl $hmmpath $uid.extended.fa`
    outfile = "$uid.nhmmer.out"
    run(pipeline(cmd, stdout=outfile))
    return "$uid.tbl"
end

function parse_tbl(file::String, glength::Integer)
    matches = FeatureMatch[]
    open(file, "r") do infile
        while !eof(infile)
            line = readline(infile)
            startswith(line, "#") && continue
            bits = split(line, " ", keepempty=false)
            rrnstrand = bits[12][1]
            rrnstart = parse(Int, bits[7])
            rrnstop = parse(Int, bits[8])
            fmlength = max(rrnstart, rrnstop) - min(rrnstart, rrnstop) + 1
            fmstart = rrnstrand == '+' ? mod1(rrnstart, glength) : mod1(reverse_complement(rrnstart, glength), glength)
            push!(matches, FeatureMatch(model2rrn[bits[3]], bits[3], rrnstrand, parse(Int, bits[5]), parse(Int, bits[6]),
                fmstart, fmlength, parse(Float64, bits[13])))
        end
    end
    @debug matches
    rationalise_matches!(matches, glength)
end

function find_closest_downstream_trn(target, trns, glength)
    idx = argmin(abs.([closestdistance(target, trn.fm.target_from+trn.fm.target_length, glength) for trn in trns]))
    closest_trn = trns[idx]
    return closest_trn
end

function fix_rrn_ends!(rRNAs, ftrns, rtrns, glength)
    for rrn in rRNAs
        trns = rrn.strand == '+' ? ftrns : rtrns
        rrnstart = rrn.target_from
        rrnstop = rrn.target_from + rrn.target_length - 1
        @debug "$rrnstart $rrnstop"
        changed = false
        trnidx = searchsortedfirst(trns, rrnstart, lt=(t,x)->t.fm.target_from < x)
        upstream_tRNA = trns[mod1(trnidx-1, length(trns))]
        @debug upstream_tRNA
        clockwise_dist = circulardistance(upstream_tRNA.fm.target_from + upstream_tRNA.fm.target_length, rrnstart, glength)
        anticlockwise_dist = circulardistance(rrnstart, upstream_tRNA.fm.target_from + upstream_tRNA.fm.target_length, glength)
        if clockwise_dist < 100 || anticlockwise_dist < rrnstart + rrn.target_length #overlap, or near enough to assume contiguity
            rrnstart = mod1(upstream_tRNA.fm.target_from + upstream_tRNA.fm.target_length, glength)
            changed = true
        end
        @debug "$rrnstart $rrnstop $changed"
        trnidx = searchsortedfirst(trns, rrnstop, lt=(t,x)->t.fm.target_from+t.fm.target_length < x)
        #downstream_tRNA = trns[mod1(trnidx, length(trns))]
        downstream_tRNA = find_closest_downstream_trn(rrnstop, trns, glength)
        @debug downstream_tRNA
        if abs(closestdistance(downstream_tRNA.fm.target_from, rrnstop, glength)) < 50 #overlap, or near enough to assume contiguity
            rrnstop = rrnstart + circulardistance(rrnstart, downstream_tRNA.fm.target_from, glength) - 1
            changed = true
        end
        @debug "$rrnstart $rrnstop $changed"
        if changed
            replace!(rRNAs, rrn => FeatureMatch(rrn.id, rrn.query, rrn.strand, rrn.model_from, rrn.model_to, rrnstart,
                circulardistance(rrnstart, rrnstop, glength) + 1, rrn.evalue))
        end
    end
end
