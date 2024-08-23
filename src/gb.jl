
function feature_compare(x::GFF, y::GFF)
    feature_hierarchy = Dict("gene" => 1, "CDS" => 2, "mRNA" => 3, "tRNA" => 4, "rRNA" => 5)
    return feature_hierarchy[x.ftype] < feature_hierarchy[y.ftype]
end

function writeGB(outfile_gb::String, uid::UUID, id::AbstractString, translation_table::Integer, gffs::Vector{GFF})
    open(outfile_gb, "w") do out
        writeGB(out, uid, id, translation_table, gffs)
    end
end
function writeGB(out::IO, uid::UUID, id::AbstractString, translation_table::Integer, gffs::Vector{GFF})
    function triplet(a, b, c)
        write(out, join([a, b, c], '\t'), '\n')
    end
    function double(a, b)
        write(out, join([a, b], '\t'), '\n')
    end
    function skip3(a, b)
        write(out, '\t'^3, "$(a)\t$(b)", '\n')
    end
    write(out, ">Feature $id\n")
    #gffs = group_features(gffs) #Ensures that features of common locus are written together
    #for (name, vector) in gffs
    #sort!(vector, lt = feature_compare) #Ensures correct hierarchy in file
    for gff in gffs
        gff.strand == '+' ? fstart = gff.fstart : fstart = gff.fend
        gff.strand == '+' ? fend = gff.fend : fend = gff.fstart
        attributes = get_attributes(gff)
        if gff.ftype == "CDS"
            if "CDS contains a frameshift where one nucleotide is skipped" in values(attributes) #For special frameshift cases
                triplet(fstart, fend, gff.ftype)
                cdsend = parse(Int64, fend)
                second_cds = vector[findfirst(x -> x.fstart == string(cdsend + 2), vector)]
                double(second_cds.fstart, second_cds.fend)
                skip3("exception", "ribosomal slippage")
                deleteat!(vector, findfirst(x -> x == second_cds, vector))
            else
                triplet(fstart, fend, gff.ftype)

            end
            skip3("product", attributes["Product"])
            skip3("transl_table", translation_table)

            if "putative non-standard start codon" in values(attributes)
                start_position = "$fstart..$(parse(Int,fstart)+2)"
                if gff.strand == '-'
                    start_position = "$(parse(Int,fstart-2))..$fstart"
                end
                skip3("/transl_except", "(pos:$start_position,aa:Met)")
            end
            #generate protein_id; should be identical to gene_id in the .gff output
            skip3("protein_id", "gnl|Emma|$(uuid5(uid, attributes["Name"]))")
            notes = sort!(filter(x -> startswith(first(x), "Note"), collect(attributes)))
            if length(notes) > 0
                skip3("note", "$(join(last.(notes), ", "))")
            end

        elseif gff.ftype == "tRNA"
            triplet(fstart, fend, gff.ftype)
            skip3("product", attributes["Product"])

        elseif gff.ftype == "rRNA"
            triplet(fstart, fend, gff.ftype)
            skip3("product", attributes["Product"])

        elseif gff.ftype == "gene"
            triplet(fstart, fend, gff.ftype)
            skip3("gene", attributes["Name"])

        elseif gff.ftype == "mRNA"
            #= note = nothing
            if !isempty(attributes)
                if "Note=5' incomplete" in [m.match for m in attributes]
                    fstart = '<' * fstart
                end
                if "Note=3' incomplete" in [m.match for m in attributes]
                    fend = '>' * fend
                end
                note = join([last(split(m.match, '=')) for m in attributes], ", ")
            end =#
            triplet(fstart, fend, gff.ftype)
            skip3("gene", attributes["Name"])
            #= if note != nothing
                write(out, '\t'^3, "note\t$note", '\n')
            end =#
        end
    end
    #end
end

#function write

function group_features(gffs::Vector{GFF})
    gff_dict = Dict{String,Vector{GFF}}()
    for gff in gffs
        name = match(r"Name=(\w+)", gff.attributes).captures[1]
        if haskey(gff_dict, name)
            push!(gff_dict[name], gff)
        else
            gff_dict[name] = [gff]
        end
    end
    sort!(collect(gff_dict), by=x -> minimum(parse.(Int, getproperty.(last(x), :fstart))))
end
