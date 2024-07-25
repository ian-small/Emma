struct CircularVector
    v::Vector{<:Integer}
end

@inline Base.length(cv::CircularVector) = length(cv.v)
@inline Base.getindex(cv::CircularVector, i::Int32) = @inbounds getindex(cv.v, mod1(i, length(cv)))
@inline Base.getindex(cv::CircularVector, is::Vector{<:Integer}) = getindex(cv.v, is .% length(cv.v))
function Base.getindex(cv::CircularVector, r::UnitRange{<:Integer})
    len = length(cv)
    start = mod1(r.start, len)
    stop = mod1(r.stop, len)
    if start < stop
        return cv.v[start:stop]
    else
        return vcat(cv.v[start:end], cv.v[1:stop])
    end
end
@inline Base.setindex!(cv::CircularVector, value::Integer, i::Integer) = @inbounds setindex!(cv.v, value, mod1(i, length(cv)))
@inline Base.push!(cv::CircularVector, value::Integer) = @inbounds push!(cv.v, value)
@inline Base.sum(cv::CircularVector) = sum(cv.v)
function Base.sum(cv::CircularVector, r::UnitRange{<:Integer})
    sum = 0
    for i in r
        sum += cv[i]
    end
    return sum
end

# clockwise distance
function circulardistance(x::Integer, y::Integer, c::Integer)
    x = mod1(x, c)
    y = mod1(y, c)
    return y >= x ? y - x : c + y - x
end

function closestdistance(x::Integer, y::Integer, c::Integer)
    clockwisedistance = circulardistance(x, y, c)
    return clockwisedistance < c/2 ? clockwisedistance : clockwisedistance - c
end

function circularin(x::Integer, start::Integer, length::Integer, c::Integer)
    r = range(mod1(start, c), length=length)
    mod1(x, c) ∈ r && return true
    mod1(x, c) + c ∈ r && return true
    false
end

function circularintersect(r1::UnitRange{<:Integer}, r2::UnitRange{<:Integer}, c::Int)::UnitRange{Int}
    maxintersect = UnitRange{Int32}(0, -1)
    maxintersectlength = 0
    i1 = intersect(r1, r2)
    if length(i1) > maxintersectlength
        maxintersect = i1
        maxintersectlength = length(i1)
    end
    if r2.stop > c
        i2 = intersect(r1 .+ c, r2)
        if length(i2) > maxintersectlength
            maxintersect = i2
            maxintersectlength = length(i2)
        end
    end
    if r1.stop > c
        i3 = intersect(r1, r2 .+ c)
        if length(i3) > maxintersectlength
            maxintersect = i3
            maxintersectlength = length(i3)
        end
    end
    if maxintersect.start > c
        maxintersect = maxintersect .- c
    end
    return maxintersect
end

struct CircularSequence
    length::Int
    sequence::LongDNA{4}

    function CircularSequence(seq::LongDNA{4})
        new(length(seq), append!(seq, LongSubSeq(seq, 1:length(seq)-1)))
    end
end

@inline Base.length(cs::CircularSequence) = cs.length
@inline Base.getindex(cs::CircularSequence, i::Integer) = @inbounds getindex(cs.sequence, mod1(i, cs.length))

function Base.getindex(cs::CircularSequence, r::UnitRange{<:Integer})
    if r.start > length(cs) || r.start < 1
        r = range(mod1(r.start, cs.length); length=length(r))
    end
    return LongSubSeq(cs.sequence, r)
end

function reverse_complement(cs::CircularSequence)
    rc = BioSequences.reverse_complement(cs.sequence[1:cs.length])
    return CircularSequence(rc)
end

function reverse_complement(i::Integer, glength::Integer)
    return glength - mod1(i,glength) + 1
end

@inline function getcodon(cs::CircularSequence, index::Integer)
    return (cs[index], cs[index+Int32(1)], cs[index+Int32(2)])
end

