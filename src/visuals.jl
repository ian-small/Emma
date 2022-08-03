using Luxor

const gradius = 200
const gcentre = Point(0, 0)
const arrowdepth = 50
const cdscolour = "gold"
const trncolour = "blue"
const rrncolour = "red"

function drawfeature(feature::GFF, glength::Integer)

    pos2radians(position) = 2π*position/glength

    start = parse(Int, feature.fstart)
    stop = parse(Int, feature.fend)
    fradius = 0
    flength = stop-start
    if feature.strand == '+'
        frame=mod(start,3)
        fradius = frame*10+gradius+5
        outside = fradius + 4
        inside = fradius - 4
        move(polar(outside, pos2radians(start)))
        arc(gcentre, outside, pos2radians(start), pos2radians(stop-min(arrowdepth,flength)))
        line(polar(fradius, pos2radians(stop)))
        line(polar(inside, pos2radians(stop-min(arrowdepth,flength))))
        carc(gcentre, inside, pos2radians(stop-min(arrowdepth,flength)), pos2radians(start))
    else
        frame=mod(reverse_complement(stop, glength),3)
        fradius = gradius - frame*10-5
        outside = fradius - 4 
        inside = fradius + 4
        move(polar(outside, pos2radians(stop)))
        carc(gcentre, outside, pos2radians(stop), pos2radians(start+min(arrowdepth,flength)))
        line(polar(fradius, pos2radians(start)))
        line(polar(inside, pos2radians(start+min(arrowdepth,flength))))
        arc(gcentre, inside, pos2radians(start+min(arrowdepth,flength)), pos2radians(stop))
    end
    closepath()
    if feature.ftype == "CDS"
        setcolor(cdscolour)
    elseif feature.ftype == "tRNA"
        setcolor(trncolour)
    elseif feature.ftype == "rRNA"
        setcolor(rrncolour)
    end
    fillpath()
    name = first(split(feature.attributes, ";"))[6:end]
    fontsize(6)
    setcolor("black")
    textcurvecentred(name, pos2radians((start+stop)/2), fradius, gcentre; baselineshift=-2)
end

function drawgenome(svgfile::String, id::String, glength::Integer, gffs::Vector{GFF})
    Drawing(500,500,svgfile)
    origin()
    background("white")
    setline(0.5)
    setcolor("black")
    circle(gcentre, gradius, action = :stroke)
    fontsize(12)
    fontface("Helvetica")
    text(id, Point(0, -20), halign=:center, valign=:center)
    text(string(glength) * " bp", halign=:center, valign=:center)
    text("Annotated by Emma", Point(0, 20), halign=:center, valign=:center)
    fontface("Helvetica Oblique")
    for gff in gffs
        drawfeature(gff, glength)
    end
    finish()
end