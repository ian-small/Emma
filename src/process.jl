const minORF = 150


function remove_starts_in_tRNAs!(starts, codons, tRNAs, glength)
    fms = getproperty.(tRNAs, :fm)
    for (startvector, codonvector) in zip(starts, codons)
        startsintRNA = []
        for (i, s) in enumerate(startvector)
            if any(circularin.(s, fms, glength))
                push!(startsintRNA, i)
            end
        end
        deleteat!(startvector, startsintRNA)
    end
end

function add_stops_at_tRNAs!(stops, tRNAs, glength)
    for t in tRNAs
        push!(stops[mod1(t.fm.target_from, 3)], t.fm.target_from)
        minus1 = mod1(t.fm.target_from - 1, glength)
        push!(stops[mod1(minus1, 3)], minus1)
        minus2 = mod1(t.fm.target_from - 2, glength)
        push!(stops[mod1(minus2, 3)], minus2)
    end
    sort!(stops[1])
    sort!(stops[2])
    sort!(stops[3])
    return stops
end

#If duplicate features exist, remove the ones with higher e-values
function remove_duplicate_features(matches::Vector)
    features = FeatureMatch[]
    sorted_array = sort(matches, by=x -> x.evalue)
    filtered_dict = Dict()
    for item in sorted_array
        if !haskey(filtered_dict, item.query)
            filtered_dict[item.query] = item
        end
    end
    filtered_array = collect(values(filtered_dict))
    for feature in filtered_array
        push!(features, feature)
    end
    return features
end

function rotate(rotate_to::String, GFFs, genome::CircularSequence, glength)
    feature_idx = findfirst(x -> occursin(rotate_to, x.attributes), GFFs)
    if feature_idx ≠ nothing
        first_feature = GFFs[feature_idx]
        if first_feature.strand == '-'
            genome = reverse_complement(genome)
            for gff in GFFs
                reverse_complement!(gff, glength)
            end
        end
        offset = parse(Int32, first_feature.fstart) - 1
        for gff in GFFs
            relocate!(gff, offset, glength)
        end
        firstchunk = LongDNA{4}(genome[offset+1:glength])
        secondchunk = dna""
        if offset > 0
            secondchunk = LongDNA{4}(genome[1:offset])
        end
        genome = CircularSequence(append!(firstchunk, secondchunk))
        return GFFs, genome
    else
        @warn "Positional translation not possible due to missing feature: $rotate_to"
        return GFFs, genome
    end
end

function emma(infile::String; translation_table=2, rotate_to=nothing, outfile_gff=nothing, outfile_gb=nothing, outfile_fa=nothing, outfile_svg=nothing, loglevel="Info")

    global_logger(ConsoleLogger(loglevel == "debug" ? Logging.Debug : Logging.Info))

    @info "$infile"
    target = FASTA.Record()
    reader = open(FASTA.Reader, infile)
    read!(reader, target)
    id = FASTA.identifier(target)
    genome = CircularSequence(FASTA.sequence(LongDNA{4}, target))
    glength = length(genome)
    rev_genome = reverse_complement(genome)
    @info "$id\t$glength bp"

    #extend genome
    extended_genome = genome[1:glength+100]
    tempfile = TempFile(".")
    name = tempfilename(tempfile, "tmp.extended.fa")
    writer = open(FASTA.Writer, name)
    write(writer, FASTA.Record(id, extended_genome))
    close(writer)

    #find tRNAs
    trn_matches = parse_trn_alignments(cmsearch(tempfile, "trn", "all_trn.cm"), glength)
    filter!(x -> x.fm.target_from <= glength, trn_matches)
    rationalise_trn_alignments(trn_matches)
    @debug trn_matches
    ftrns = get_best_trns(filter(x -> x.fm.strand == '+', trn_matches), glength)
    rtrns = get_best_trns(filter(x -> x.fm.strand == '-', trn_matches), glength)
    #check for overlapping tRNAs
    overlapped = get_overlapped_trns(sort(ftrns, by=x -> x.fm.target_from), glength)
    append!(overlapped, get_overlapped_trns(sort(rtrns, by=x -> x.fm.target_from), glength))
    trn_matches = append!(ftrns, rtrns)
    #for overlapped tRNAs, generate polyadenylated version
    for (cma, trunc_end) in overlapped
        trnseq = cma.fm.strand == '+' ? genome.sequence[cma.fm.target_from:trunc_end] : rev_genome.sequence[cma.fm.target_from:trunc_end]
        trnseq_polyA = trnseq * dna"AAAAAAAAAA"
        polyA_matches = parse_trn_alignments(cmsearch(tempfile, cma.fm.query, "trn", trnseq_polyA), 0)
        isempty(polyA_matches) && continue
        trn_match = polyA_matches[1]
        if trn_match.fm.evalue < cma.fm.evalue
            #construct modified CMA with new evalue, new end point, polyA length
            newcma = tRNA(FeatureMatch(cma.fm.id, cma.fm.query, cma.fm.strand, trn_match.fm.model_from, trn_match.fm.model_to, cma.fm.target_from,
                    circulardistance(cma.fm.target_from, trunc_end, glength) + 1, trn_match.fm.evalue),
                trn_match.anticodon, (trn_match.fm.target_length) - (trunc_end - cma.fm.target_from))
            #delete old match from trn_matches
            deleteat!(trn_matches, findfirst(x -> x == cma, trn_matches))
            #add new CMA
            push!(trn_matches, newcma)
        end
    end

    @debug trn_matches
    @info "found $(length(trn_matches)) tRNA genes"

    #find rRNAs
    #search for rrns using hmmsearch
    rrns = parse_tbl(rrnsearch(tempfile), glength)
    filter!(x -> x.evalue < 1e-10, rrns)
    @debug rrns
    #fix ends using flanking trn genes
    ftRNAs = filter(x -> x.fm.strand == '+', trn_matches)
    sort!(ftRNAs, by=x -> x.fm.target_from)
    rtRNAs = filter(x -> x.fm.strand == '-', trn_matches)
    sort!(rtRNAs, by=x -> x.fm.target_from)
    fix_rrn_ends!(rrns, ftRNAs, rtRNAs, glength)
    @debug rrns

    #find CDSs
    startcodon = ncbi_start_codons[translation_table]
    stopcodon = ncbi_stop_codons[translation_table]
    fstarts, fstartcodons = getcodons(genome, startcodon)
    remove_starts_in_tRNAs!(fstarts, fstartcodons, ftRNAs, glength)
    fstops = codonmatches(genome, stopcodon)
    add_stops_at_tRNAs!(fstops, ftRNAs, glength)
    @debug fstops

    rstarts, rstartcodons = getcodons(rev_genome, startcodon)
    remove_starts_in_tRNAs!(rstarts, rstartcodons, rtRNAs, glength)
    rstops = codonmatches(rev_genome, stopcodon)
    add_stops_at_tRNAs!(rstops, rtRNAs, glength)

    cds_matches = parse_domt(orfsearch(tempfile, id, genome, translation_table, fstarts, fstops, rstarts, rstops, minORF), glength)
    @debug cds_matches
    #fix start & stop codons
    fhmms = filter(x -> x.strand == '+', cds_matches)
    rhmms = filter(x -> x.strand == '-', cds_matches)
    fix_start_and_stop_codons!(fhmms, ftRNAs, fstarts, fstartcodons, fstops, glength)
    fix_start_and_stop_codons!(rhmms, rtRNAs, rstarts, rstartcodons, rstops, glength)

    cds_matches = append!(fhmms, rhmms)
    frameshift_merge!(cds_matches, glength, genome)
    @info "found $(length(cds_matches)) protein-coding genes"
    gffs = getGFF(tempfile.uuid, genome, rev_genome, cds_matches, trn_matches, rrns, glength)

    if ~isnothing(rotate_to)
        gffs, genome = rotate(rotate_to, gffs, genome, glength)
    end

    if ~isnothing(outfile_fa)
        open(FASTA.Writer, outfile_fa) do w
            write(w, FASTA.Record(id, LongDNA{4}(genome[1:glength])))
        end
    end

    if ~isnothing(outfile_gff)
        writeGFF(id, gffs, outfile_gff, glength)
    end

    if ~isnothing(outfile_gb)
        writeGB(tempfile.uuid, id, translation_table, gffs, outfile_gb, glength)
    end

    if ~isnothing(outfile_svg)
        drawgenome(outfile_svg, id, glength, gffs)
    end

    ##cleanup
    cleanfiles(tempfile)
end
