
gene2products = Dict("ND1"=>"NADH dehydrogenase subunit 1","nad1"=>"NADH dehydrogenase subunit 1","ND2"=>"NADH dehydrogenase subunit 2","nad2"=>"NADH dehydrogenase subunit 2","ND3"=>"NADH dehydrogenase subunit 3",
    "nad3"=>"NADH dehydrogenase subunit 3","ND4"=>"NADH dehydrogenase subunit 4","nad4"=>"NADH dehydrogenase subunit 4","ND4L"=>"NADH dehydrogenase subunit 4L","nad4L"=>"NADH dehydrogenase subunit 4L",
    "ND5"=>"NADH dehydrogenase subunit 5","nad5"=>"NADH dehydrogenase subunit 5","ND6"=>"NADH dehydrogenase subunit 6","nad6"=>"NADH dehydrogenase subunit 6","COX1"=>"cytochrome c oxidase subunit I",
    "cox1"=>"cytochrome c oxidase subunit I","COX2"=>"cytochrome c oxidase subunit II","cox2"=>"cytochrome c oxidase subunit II","COX3"=>"cytochrome c oxidase subunit III","cox3"=>"cytochrome c oxidase subunit III",
    "ATP6"=>"ATP synthase F0 subunit 6","atp6"=>"ATP synthase F0 subunit 6","ATP8"=>"ATP synthase F0 subunit 8","atp8"=>"ATP synthase F0 subunit 8","CYTB"=>"cytochrome b","cob"=>"cytochrome b",
    "trnA"=>"tRNA-Ala","trnR"=>"tRNA-Arg","trnN"=>"tRNA-Asn","trnD"=>"tRNA-Asp","trnC"=>"tRNA-Cys","trnQ"=>"tRNA-Gln","trnE"=>"tRNA-Glu","trnG"=>"tRNA-Gly","trnH"=>"tRNA-His","trnI"=>"tRNA-Ile","trnK"=>"tRNA-Lys",
    "trnM"=>"tRNA-Met","trnF"=>"tRNA-Phe","trnP"=>"tRNA-Pro","trnT"=>"tRNA-Thr","trnW"=>"tRNA-Trp","trnY"=>"tRNA-Tyr",
    "trnV"=>"tRNA-Val","trnS1"=>"tRNA-Ser","trnS2"=>"tRNA-Ser","trnL1"=>"tRNA-Leu","trnL2"=>"tRNA-Leu","12srna"=>"12S rRNA","rrn12"=>"12S rRNA","16srna"=>"16S rRNA","rrn16"=>"16S rRNA")

function feature_compare(x::GFF, y::GFF)
    feature_hierarchy = Dict("gene" => 1, "CDS" => 2, "mRNA" => 3, "tRNA" => 4, "rRNA" => 5)
    return feature_hierarchy[x.ftype] < feature_hierarchy[y.ftype]
end

function writeGB(id::AbstractString, gffs::Vector{GFF}, outfile_gb::String, glength::Integer)
    open(outfile_gb, "w") do out
        write(out, ">Feature $id\n")
        gffs = group_features(gffs) #Ensures that features of common locus are written together
        for (name, vector) in gffs
            sort!(vector, lt = feature_compare) #Ensures correct heirarchy in file
            for gff in vector
                gff.strand == '+' ? fstart = gff.fstart : fstart = gff.fend
                gff.strand == '+' ? fend = gff.fend : fend = gff.fstart
                attributes = eachmatch(r"Note=([^;]+)", gff.attributes) #Matches notes and adds them to attributes vector
                if gff.ftype == "CDS"
                    if "Note=CDS contains a frameshift where one nucleotide is skipped" in [m.match for m in attributes] #For special frameshift cases
                        write(out, join([fstart, fend, gff.ftype], '\t'), '\n')
                        cdsend = parse(Int64, fend)
                        second_cds = vector[findfirst(x->x.fstart==string(cdsend+2),vector)]
                        write(out, join([second_cds.fstart, second_cds.fend], '\t'), '\n')
                        write(out, '\t'^3, "exception\tribosomal slippage\n")
                        deleteat!(vector, findfirst(x->x==second_cds,vector))
                    else
                        write(out, join([fstart, fend, gff.ftype], '\t'), '\n')
                    end
                    write(out, '\t'^3, "gene\t$name", '\n')
                    write(out, '\t'^3, "product\t$(gene2products[name])", '\n')
                    write(out, '\t'^3, "transl_table\t2", '\n')
                    if any(startswith.([m.match for m in attributes], "Note=putative non-standard start codon"))
                        start_position = "$fstart..$(parse(Int,fstart)+2)"
                        if gff.strand == '-'
                            start_position = "$(parse(Int,fstart-2))..$fstart"
                        end
                        write(out, '\t'^3, "/transl_except\t(pos:$start_position,aa:Met)", '\n')
                    end
                    if !isempty(attributes)
                        note = join([last(split(m.match, '=')) for m in attributes], ", ") #Appends notes into single string
                        write(out, '\t'^3, "note\t$note", '\n')
                    end
                
                elseif gff.ftype == "tRNA"
                    anticodon = match(r"trn..-(\w+)", gff.attributes) #Matches trnS/L features to add anticodon note
                    write(out, join([fstart, fend, gff.ftype], '\t'), '\n')
                    write(out, '\t'^3, "product\t$(gene2products[name])", '\n')
                    if ~isnothing(anticodon)
                        write(out, '\t'^3, "note\tAnticodon: $(anticodon.captures[1])", '\n')
                    end
                
                elseif gff.ftype == "rRNA"
                    write(out, join([fstart, fend, gff.ftype], '\t'), '\n')
                    write(out, '\t'^3, "product\t$(gene2products[name])", '\n')
                
                elseif gff.ftype == "gene"
                    write(out, join([fstart, fend, gff.ftype], '\t'), '\n')
                    write(out, '\t'^3, "gene\t$name", '\n')
                
                elseif gff.ftype == "mRNA"
                    note = nothing
                    if !isempty(attributes)
                        if "Note=5' incomplete" in [m.match for m in attributes]
                            fstart = '<' * fstart
                        end
                        if "Note=3' incomplete" in [m.match for m in attributes]
                            fend = '>' * fend
                        end
                        note = join([last(split(m.match, '=')) for m in attributes], ", ")
                    end
                    write(out, join([fstart, fend, gff.ftype], '\t'), '\n')
                    write(out, '\t'^3, "gene\t$name", '\n')
                    write(out, '\t'^3, "product\t$(gene2products[name])", '\n')
                    if note != nothing
                        write(out, '\t'^3, "note\t$note", '\n')
                    end
                end
            end
        end
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