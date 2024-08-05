mutable struct GFF
    source::String
    ftype::String
    fstart::String
    fend::String
    score::String
    strand::Char
    phase::String
    attributes::String
end


function FMcoords2GFF(strand::Char, start::Integer, length::Integer, glength::Integer)
    gffstart = strand == '+' ? start : mod1(reverse_complement(start + length - 1, glength), glength)
    gffend = strand == '+' ? start + length - 1 : reverse_complement(start, glength)
    if gffend < gffstart; gffend += glength; end
    gffstart, gffend
end

function FMcoords2GFF(fm::FeatureMatch, glength::Integer)
    FMcoords2GFF(fm.strand, fm.target_from, fm.target_length, glength::Integer)
end

function CDS2GFF(cds::FeatureMatch, genome::CircularSequence, rev_genome::CircularSequence, trns::Vector{tRNA})
    glength = length(genome)
    attributes = ""
    cdsstart = cds.target_from
    cdsstop = cds.target_from + cds.target_length - 1 + 3
    startcodon = cds.strand == '+' ? genome[cdsstart:cdsstart+2] : rev_genome[cdsstop-2:cdsstop]
    if startcodon ∈ [dna"CTG", dna"TTG"]
        attributes *= ";Note=putative non-standard start codon $startcodon"
    end
    trnidx = findfirst(t->circularin(t.fm.target_from, cdsstop-2, 3, glength), trns)
    if !isnothing(trnidx)
        stopcodon = cds.strand == '+' ? copy(genome[cdsstop-2:cdsstop-2+circulardistance(cdsstop-2, trns[trnidx].fm.target_from, glength)-1]) : copy(rev_genome[cdsstop-2:cdsstop-2+circulardistance(cdsstop-2, trns[trnidx].fm.target_from, glength)-1])
        cdsstop = trns[trnidx].fm.target_from - 1
        while length(stopcodon) < 3
            push!(stopcodon, DNA_A)
        end
        attributes *= ";Note=putative $stopcodon stop codon is completed by the addition of 3' A residues to the mRNA"
    end
    gffstart, gffend = FMcoords2GFF(cds.strand, cdsstart, circulardistance(cdsstart, cdsstop+1, glength), glength)
    return GFF("Emma", "CDS", string(gffstart), string(gffend), string(cds.evalue), cds.strand, "0", attributes)
end

function tRNA2GFF(trn::tRNA, glength::Integer)
    gffstart, gffend = FMcoords2GFF(trn.fm, glength)
    attributes = ""
    if trn.polyA > 0
        attributes *= ";Note=tRNA completed by post-transcriptional addition of " * string(trn.polyA)
        attributes *= trn.polyA > 1 ? " As" : " A"
    end
    if trn.fm.query == "trnD" && trn.anticodon == "GCC"
        attributes *= ";Note=GUC anticodon completed by RNA editing"
        @warn "trnD $(trn.anticodon) edited to GUC by RNA editing"
    elseif !haskey(anticodon2trn, trn.anticodon)
        attributes *= ";Note=tRNA has no valid anticodon"
        @warn "no valid anticodon found for $(name)"
    end
    if typeof(attributes) != Missing
        return GFF("Emma", "tRNA", string(gffstart), string(gffend), string(trn.fm.evalue), trn.fm.strand, ".", attributes)
    else 
        return nothing
    end
end

function rRNA2GFF(rrn::FeatureMatch, glength::Integer)
    gffstart, gffend = FMcoords2GFF(rrn, glength)
    return GFF("Emma", "rRNA", string(gffstart), string(gffend), string(rrn.evalue), rrn.strand, ".", "")
end

#= function add_geneGFFs(gffs::Vector{GFF}, genome::CircularSequence, glength::Integer)
    #Adds gene features, that are essentially identical to CDS features
    cdsGFFs = filter(x->x.ftype == "CDS", gffs)
    for gff in cdsGFFs
        name = cds2symbol[match(r"Name=(\w+)", gff.attributes).captures[1]]
        if name == "MT-ND3"
            group = filter(x -> x != gff && x in gffs && (x.fstart == gff.fstart || x.fend == gff.fend), cdsGFFs)
            if length(group) == 2
                deleteat!(gffs, findfirst(x->x==gff,gffs))
                g1 = first(group)
                g2 = last(group)
                deleteat!(gffs, findfirst(x->x==g1,gffs))
                deleteat!(gffs, findfirst(x->x==g2,gffs))
                push!(gffs, GFF("Emma", "gene", gff.fstart, gff.fend, gff.score, gff.strand, ".", "Name=$name;ID=$name.gene"))
                push!(gffs, GFF("Emma", "CDS", g1.fstart, g1.fend, g1.score, g1.strand, ".", "Name=$name;ID=$name.CDS.1;Note=CDS contains a frameshift where one nucleotide is skipped"))
                push!(gffs, GFF("Emma", "CDS", g2.fstart, g2.fend, g2.score, g2.strand, ".", "Name=$name;ID=$name.CDS.2;Note=CDS contains a frameshift where one nucleotide is skipped"))
            else
                push!(gffs, GFF("Emma", "gene", gff.fstart, gff.fend, gff.score, gff.strand, ".", "Name=$name;ID=$name.gene"))
            end
        else
            push!(gffs, GFF("Emma", "gene", gff.fstart, gff.fend, gff.score, gff.strand, ".", "Name=$name;ID=$name.gene"))
        end
    end
    return gffs
end =#

function getGFF(uid::UUID, genome::CircularSequence, rev_genome::CircularSequence, cds_matches::Vector{FeatureMatch},
             trn_matches::Vector{tRNA}, rRNAs::Vector{FeatureMatch}, glength::Integer)

    gffs = GFF[]
    genome_length = length(genome)    

    function writeone(gff::Union{Nothing, GFF})
        if ~isnothing(gff)
            push!(gffs, gff)
        end
    end

    for trn in trn_matches
        #gene
        name = trn2symbol[trn.fm.query]
        gene_id = uuid5(uid, name)
        #tRNA
        gff = tRNA2GFF(trn, genome_length)
        writeone(GFF("Emma", "gene", gff.fstart, gff.fend, ".", trn.fm.strand, ".", "ID=$gene_id;Name=$name"))
        trn_id = uuid5(gene_id, name)
        gff = tRNA2GFF(trn, glength)
        gff.attributes = "ID=$trn_id;Parent=$gene_id;Name=$name" * gff.attributes
        writeone(gff)
    end

    for cds in cds_matches
        #gene
        name = cds2symbol[first(split(cds.query, '.'))]
        gene_id = uuid5(uid, name)
        #mRNA
        mrna_id = uuid5(gene_id, name)
        #CDS
        gff = CDS2GFF(cds, genome, rev_genome, trn_matches)
        writeone(GFF("Emma", "gene", gff.fstart, gff.fend, ".", cds.strand, ".", "ID=$gene_id;Name=$name"))
        writeone(GFF("Emma", "mRNA", gff.fstart, gff.fend, ".", cds.strand, ".", "ID=$mrna_id;Parent=$gene_id;Name=$name"))
        cds_id = uuid5(mrna_id, name)
        gff = CDS2GFF(cds, genome, rev_genome, trn_matches)
        gff.attributes = "ID=$cds_id;Parent=$mrna_id;Name=$name" * gff.attributes
        writeone(gff)
    end

    for rrn in rRNAs
        #gene
        name = rrn2symbol[rrn.id]
        gene_id = uuid5(uid, name)
        gff = rRNA2GFF(rrn, genome_length)
        writeone(GFF("Emma", "gene", gff.fstart, gff.fend, ".", rrn.strand, ".", "ID=$gene_id;Name=$name"))
        rrn_id = uuid5(gene_id, name)
        gff = rRNA2GFF(rrn, glength)
        gff.attributes = "ID=$rrn_id;Parent=$gene_id;Name=$name" * gff.attributes
        writeone(gff)
    end
    
    return sort!(gffs; by = x ->x.fstart)
end

function writeGFF(id::AbstractString, gffs::Vector{GFF}, outfile::String, glength::Integer)
    open(outfile, "w") do out
        write(out, "##gff-version 3\n")
        write(out, "##source-version Emma beta\n")
        write(out, "##sequence-region	$id	1	$glength\n")
        write(out, "$id	Emma	region	1	$glength	.	+	0	Is_circular=true\n")
        for gff in gffs
            write(out, join([id, gff.source, gff.ftype,gff.fstart,gff.fend,gff.score,gff.strand,gff.phase,gff.attributes], "\t"), "\n")
        end
    end
end

#= function is_encompassed(outer_gene, inner_gene)
    return outer_gene.target_from <= inner_gene.target_from && (outer_gene.target_from+outer_gene.target_length) >= (inner_gene.target_from+inner_gene.target_length)
end

function add_mRNAGFFs(cds_matches::Vector{FeatureMatch}, trns::Vector{tRNA}, glength::Integer)
    mRNAs = GFF[]
    for (i,cds_match) in enumerate(cds_matches)
        name = cds2symbol[first(split(cds_match.query, '.'))]
        attributes = "Name=" * name * ";ID=$name.mRNA"
        cdsstart = cds_match.target_from
        cdsstop = cdsstart + cds_match.target_length - 1 + 3
        trnidx = searchsortedfirst(trns, cdsstart, lt=(t,x)->t.fm.target_from < x)

        upstream_tRNA = trns[mod1(trnidx-1, length(trns))]
        upstream_trn_end = upstream_tRNA.fm.target_from + upstream_tRNA.fm.target_length - 1
        dist_to_upstream_trn = abs(closestdistance(upstream_trn_end, cdsstart, glength))

        downstream_tRNA = trns[mod1(trnidx, length(trns))]
        downstream_trn_start = downstream_tRNA.fm.target_from
        dist_to_downstream_trn = abs(closestdistance(downstream_trn_start, cdsstop, glength))
        if dist_to_upstream_trn<10
            cdsstart = upstream_trn_end + 1
        else
            attributes *= ";Note=5' incomplete"
        end
        if dist_to_downstream_trn<10
            cdsstop = downstream_trn_start - 1
        else
            attributes *= ";Note=3' incomplete"
        end
        gffstart, gffend = FMcoords2GFF(cds_match.strand, cdsstart, circulardistance(cdsstart, cdsstop+1, glength), glength)
        push!(mRNAs, GFF("Emma", "mRNA", string(gffstart), string(gffend), string(cds_match.evalue), cds_match.strand, "0", attributes))
    end
    return mRNAs      
end =#

reverse_strand = Dict('+' => '-', '-' => '+', '.' => '.')

function reverse_complement!(gff::GFF, glength::Integer)
    gff.strand = reverse_strand[gff.strand]
    fstart = parse(Int32, gff.fstart)
    fend = parse(Int32, gff.fend)
    rcstart = reverse_complement(fend, glength)
    gff.fstart = string(rcstart)
    flength = fend - fstart + 1
    gff.fend = string(rcstart + flength - 1)
end

function relocate!(gff::GFF, offset::Integer, glength::Integer)
    fstart = parse(Int32, gff.fstart)
    fend = parse(Int32, gff.fend)
    flength = fend - fstart + 1
    fstart = mod1(fstart - offset, glength)
    fend = fstart + flength - 1
    gff.fstart = string(fstart)
    gff.fend = string(fend)
end
            
