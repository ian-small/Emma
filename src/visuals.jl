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
        Luxor.line(polar(fradius, pos2radians(stop)))
        Luxor.line(polar(inside, pos2radians(stop-min(arrowdepth,flength))))
        carc(gcentre, inside, pos2radians(stop-min(arrowdepth,flength)), pos2radians(start))
    else
        frame=mod(reverse_complement(stop, glength),3)
        fradius = gradius - frame*10-5
        outside = fradius - 4 
        inside = fradius + 4
        move(polar(outside, pos2radians(stop)))
        carc(gcentre, outside, pos2radians(stop), pos2radians(start+min(arrowdepth,flength)))
        Luxor.line(polar(fradius, pos2radians(start)))
        Luxor.line(polar(inside, pos2radians(start+min(arrowdepth,flength))))
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
    prot_name = first(split(name, "."))
    fontsize(6)
    setcolor("black")

#    textcurvecentred(name, pos2radians((start+stop)/2), fradius, gcentre; baselineshift=-2)
    location=pos2radians((start+stop)/2)
#   place protine and rna names are right angles to circle
#   draw a line connecting them
    if (prot_name == "rrnS")
        (prot_name="12S rRNA")
    end
    if (prot_name == "rrnL")
        (prot_name="16S rRNA")
    end
    Luxor.text(prot_name, Point((fradius+20)*cos(location), (fradius+20)sin(location)), angle=location)
    Luxor.line(Point((fradius+5)*cos(location), (fradius+5)*sin(location)), Point((fradius+18)*cos(location), (fradius+18)sin(location)), :stroke)
end

function drawgenome(svgfile::String, id::AbstractString, glength::Integer, gffs::Vector{GFF})
    Drawing(600,600,svgfile)
    origin()
    background("white")
    setline(0.5)
    setcolor("black")
    Luxor.circle(gcentre, gradius, action = :stroke)
    fontsize(12)
    fontface("Helvetica")
    Luxor.text(id, Point(0, -20), halign=:center, valign=:center)
    Luxor.text(string(glength) * " bp", halign=:center, valign=:center)
    Luxor.text("Annotated by Emma", Point(0, 20), halign=:center, valign=:center)
    fontface("Helvetica Oblique")
    for gff in gffs
        drawfeature(gff, glength)
    end

    #draw circular bp map for full length of genome
    _counter() = (a = -1; () -> a += 1)
    counter = _counter()
        fontsize(6)
        arrow(O +  (0, 0), 140, 0, 2π,
            arrowheadlength=0,
            decoration=range(0, 1, length=glength),
            decorate = () -> begin
                    d = counter()
                    if d % 1000 == 0
                        Luxor.text(string(d), O + (0, -20), halign=:center)
                        setline(0.6)
                        Luxor.line(O - (0, 10), O + (0, 0), action = :stroke)
                    end
                    if d % 200 == 0
                        setline(0.4)
                        Luxor.line(O - (0, 3), O + (0, 0), action = :stroke)
                    end
                 end
            )

    finish()
end