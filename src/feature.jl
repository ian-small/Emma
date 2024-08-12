
# unified struct to cover all matches of HMMs or CMs; target coordinates are from the 5' end of the strand.
struct FeatureMatch
    id::String
    query::String
    strand::Char
    model_from::Int
    model_to::Int
    target_from::Int
    target_length::Int
    evalue::Float64
end

function circularin(x::Integer, f::FeatureMatch, c::Integer)
    circularin(x, f.target_from, f.target_length, c)
end

function circularoverlap(f1::FeatureMatch, f2::FeatureMatch, c::Integer)
    if circularin(f1.target_from, f2, c)
        return circulardistance(f1.target_from, min(f1.target_from + f1.target_length -1, f2.target_from + f2.target_length - 1), c)
    elseif circularin(f2.target_from, f1, c)
        return circulardistance(f2.target_from, min(f1.target_from + f1.target_length -1, f2.target_from + f2.target_length - 1), c)
    end
    return 0
end

function circulardistance(m1::FeatureMatch, m2::FeatureMatch, c::Integer)
    return circulardistance(m1.target_from, m2.target_from, c)
end

function merge_matches(m1::FeatureMatch, m2::FeatureMatch, glength::Integer)
    return FeatureMatch(m1.id, m1.query, m1.strand, min(m1.model_from, m2.model_from), max(m1.model_to, m2.model_to), m1.target_from,
        max(m1.target_length, circulardistance(m1.target_from, m2.target_from + m2.target_length, glength)), min(m1.evalue, m2.evalue))
end

function get_ND3s(merged::FeatureMatch, glength::Integer, position::Integer)
    gene1 = FeatureMatch(merged.id, merged.query, merged.strand, merged.model_from, merged.model_from + position-2, merged.target_from,
    circulardistance(merged.target_from, merged.target_from + position-2, glength), merged.evalue)
    gene2 = FeatureMatch(merged.id, merged.query, merged.strand, merged.model_from + position + 2, merged.model_to, merged.target_from + position + 2,
    circulardistance(merged.target_from + position+2, merged.target_from + merged.target_length, glength), merged.evalue)
    return (gene1, gene2)
end

function five2three(x, y, glength)
    circulardistance(y, x, glength) > circulardistance(x, y, glength)
end

# merge adjacent matches to same model
function rationalise_matches!(matches::Vector{FeatureMatch}, glength::Integer)::Vector{FeatureMatch}
    length(matches) < 2 && return matches
    filter!(x -> x.evalue < 1e-5, matches)
    for i in matches, j in matches
        i == j && continue
        ij = [i,j]
        sort!(ij, lt=(x,y)->five2three(x,y,glength))
        i = first(ij)
        j = last(ij)      
        # if not same gene and substantially overlap, only keep the best one
        if i.query ≠ j.query
            if i.strand == j.strand && circularoverlap(i, j, glength) > 0.5 * min(i.target_length, j.target_length)
                @debug i, j, circularoverlap(i, j, glength), 0.5 * min(i.target_length, j.target_length)
                todelete = i.evalue > j.evalue ? i : j
                deleteat!(matches, findfirst(x->x==todelete,matches))
                return(rationalise_matches!(matches, glength))
            else
                continue
            end
        end
        # if they are the same gene and close, merge them; else delete the poorer match
        modeldistance = j.model_from - i.model_to
        modellength = max(i.model_to, j.model_to) - min(i.model_from, j.model_from)
        matchdistance = abs(closestdistance(i.target_from + i.target_length - 1, j.target_from, glength))
        tolerance = 0.1
        if (in_frame(i, j, glength) || i.query == "16srna" || i.query == "12srna") && (circularin(j.target_from, i, glength) || (matchdistance - modeldistance)^2 < tolerance * modellength^2) # arbitrary tolerance
            merged_match = merge_matches(i, j, glength)
            deleteat!(matches, findfirst(x->x==i,matches))
            deleteat!(matches, findfirst(x->x==j,matches))
            push!(matches, merged_match)
        else
            todelete = i.evalue > j.evalue ? i : j
            deleteat!(matches, findfirst(x->x==todelete,matches))
        end
        return(rationalise_matches!(matches, glength))
    end
    matches
end

function group_duplicates(matches::Vector{FeatureMatch})
    query_dict = Dict{String, Vector{FeatureMatch}}()
    for item in matches
        query_value = item.query
        if haskey(query_dict, query_value)
            push!(query_dict[query_value], item)
        else
            query_dict[query_value] = [item]
        end
    end
    return query_dict
end

#Check if two Features are in-frame. Assumes they are 5-to-3 sorted.
function in_frame(f1::FeatureMatch, f2::FeatureMatch, glength::Integer)
    if f2.target_from >= f1.target_from
        return mod1(f1.target_from, 3) == mod1(f2.target_from, 3)
    else
        return mod1(f1.target_from, 3) == mod1(glength + f2.target_from, 3)
    end
end

function frameshift_merge!(matches::Vector{FeatureMatch}, glength::Integer, genome::CircularSequence)
    #first identify duplicate cds features
    query_dict = group_duplicates(matches)
    for (key, vector) in query_dict
        if length(vector) > 1
            #Pairs duplicates for checks
            for i in vector, j in vector
                i == j && continue
                i.strand ≠ j.strand && continue
                ij = [i,j]
                #Ensures pairs are ordered 5' - 3'
                sort!(ij, lt=(x,y)->five2three(x,y,glength))
                i = first(ij)
                j = last(ij)
                modeldistance = j.model_from - i.model_to
                modellength = max(i.model_to, j.model_to) - min(i.model_from, j.model_from)
                matchdistance = abs(closestdistance(i.target_from + i.target_length - 1, j.target_from, glength))
                tolerance = 0.1
                if circularin(j.target_from, i, glength) || (matchdistance - modeldistance)^2 < tolerance * modellength^2 # arbitrary tolerance
                    merged_match = merge_matches(i, j, glength)
                    merge_start = merged_match.target_from
                    merge_end = merged_match.target_from + merged_match.target_length - 1
                    #Check merged hit for conserved frameshift sequence
                    merged_seq = i.strand == '+' ? genome.sequence[merge_start:merge_end] : reverse_complement(genome).sequence[merge_start:merge_end]
                    sequence_match = match(r".TT.CT.AGTAGC", String(merged_seq))
                    if sequence_match != nothing
                        position = sequence_match.offset+5
                        gene1, gene2 = get_ND3s(merged_match, glength, position)
                        deleteat!(matches, findfirst(x->x==i,matches))
                        deleteat!(matches, findfirst(x->x==j,matches))
                        filter!(x -> x != i, vector)
                        filter!(x -> x != j, vector)
                        @warn "frameshift detected at nt $position for $(i.query)"
                        push!(matches, merged_match)
                        push!(matches, gene1)
                        push!(matches, gene2)
                    end
                end
            end
        end
    end
end
