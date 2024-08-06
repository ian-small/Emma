
function feature_compare(x::GFF, y::GFF)
    feature_hierarchy = Dict("gene" => 1, "CDS" => 2, "mRNA" => 3, "tRNA" => 4, "rRNA" => 5)
    return feature_hierarchy[x.ftype] < feature_hierarchy[y.ftype]
end

function writeGB(uid::UUID, id::AbstractString, translation_table::Int16, gffs::Vector{GFF}, outfile_gb::String, glength::Integer)
    open(outfile_gb, "w") do out
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
                        write(out, join([fstart, fend, gff.ftype], '\t'), '\n')
                        cdsend = parse(Int64, fend)
                        second_cds = vector[findfirst(x->x.fstart==string(cdsend+2),vector)]
                        write(out, join([second_cds.fstart, second_cds.fend], '\t'), '\n')
                        write(out, '\t'^3, "exception\tribosomal slippage\n")
                        deleteat!(vector, findfirst(x->x==second_cds,vector))
                    else
                        write(out, join([fstart, fend, gff.ftype], '\t'), '\n')
                    end
                    write(out, '\t'^3, "product\t$(attributes["Product"])", '\n')
                    write(out, '\t'^3, "transl_table\t$translation_table", '\n')
                    if "putative non-standard start codon" in values(attributes)
                        start_position = "$fstart..$(parse(Int,fstart)+2)"
                        if gff.strand == '-'
                            start_position = "$(parse(Int,fstart-2))..$fstart"
                        end
                        write(out, '\t'^3, "/transl_except\t(pos:$start_position,aa:Met)", '\n')
                    end
                    #generate protein_id; should be identical to gene_id in the .gff output
                    write(out, '\t'^3, "protein_id\tgnl|Emma|$(uuid5(uid, attributes["Name"]))", '\n')
                    notes = sort!(filter(x -> startswith(first(x), "Note"), collect(attributes)))
                    if length(notes) > 0
                        write(out, '\t'^3, "note\t$(join(last.(notes), ", "))", '\n')                        
                    end
                    
                elseif gff.ftype == "tRNA"
                    write(out, join([fstart, fend, gff.ftype], '\t'), '\n')
                    write(out, '\t'^3, "product\t$(attributes["Product"])", '\n')
                
                elseif gff.ftype == "rRNA"
                    write(out, join([fstart, fend, gff.ftype], '\t'), '\n')
                    write(out, '\t'^3, "product\t$(attributes["Product"])", '\n')
                
                elseif gff.ftype == "gene"
                    write(out, join([fstart, fend, gff.ftype], '\t'), '\n')
                    write(out, '\t'^3, "gene\t$(attributes["Name"])", '\n')
                
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
                    write(out, join([fstart, fend, gff.ftype], '\t'), '\n')
                    write(out, '\t'^3, "gene\t$(attributes["Name"])", '\n')
                    #= if note != nothing
                        write(out, '\t'^3, "note\t$note", '\n')
                    end =#
                end
            end
        #end
    end
end

#function write

function group_features(gffs::Vector{GFF})
    gff_dict = Dict{String, Vector{GFF}}()
    for gff in gffs
        name = match(r"Name=(\w+)", gff.attributes).captures[1]
        if haskey(gff_dict, name)
            push!(gff_dict[name], gff)
        else
            gff_dict[name] = [gff]
        end
    end
    sort!(collect(gff_dict), by = x->minimum(parse.(Int, getproperty.(last(x),:fstart))))
end