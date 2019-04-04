using VoronoiDelaunay
using VoronoiDelaunay.GeometricalPredicates
using VegaLite
using ElectronDisplay
using LinearAlgebra
using SparseArrays

using Base.Iterators: product

import VoronoiDelaunay.GeometricalPredicates: getx, gety
import VoronoiDelaunay: min_coord, max_coord

import Base: show

abstract type ParticleType <: AbstractPoint2D end

struct CellWall
    a::Point2D
    b::Point2D
    ogen::Union{Nothing,ParticleType}
end

mutable struct Particle <: ParticleType
    _x::Float64
    _y::Float64
    index::Int64
    mass::Float64
    density::Float64
    pressure::Float64
    surface::Float64
    perimeter::Float64
    vx::Float64
    vy::Float64
    ax::Float64
    ay::Float64
    walls::Vector{CellWall}
end

function show(io::IO, part::Particle)
    show(io, "#$(part.index) ($(round(part._x, digits=2)),$(round(part._y, digits=2)))" *
        ", mass $(round(part.mass, sigdigits=2))" *
        ", pressure $(round(part.pressure, sigdigits=2))" *
        ", area $(round(part.surface, sigdigits=2))" *
        ", perim $(round(part.perimeter, sigdigits=2))")
end
# pts[4]

Particle(x::Float64, y::Float64) = Particle(x,y,0)
Particle(x::Float64, y::Float64, ind::Int64) =
    Particle(x, y, ind, 1., 1., 0., 0., 0., 0., 0., 0., 0., [])
getx(p::Particle) = p._x
gety(p::Particle) = p._y


##########################
# isentropic gas, H₂O γ ~=  1.33
const γ = 1.33



##################################################################################
function intersection(a1::Point2D, a2::Point2D, b1::Point2D, b2::Point2D)::Union{Nothing,Point2D}
    A = [ getx(a2)-getx(a1) getx(b1)-getx(b2) ;
          gety(a2)-gety(a1) gety(b1)-gety(b2)  ]
    (det(A) == 0.) && return nothing
    b = [ getx(b1) - getx(a1), gety(b1) - gety(a1) ]
    x = A \ b
    all(0 .<= x .<= 1.) && return Point2D( getx(a1) + x[1] * A[1,1], gety(a1) + x[1] * A[2,1] )
    nothing
end

function norm2(a1::Point2D, a2::Point2D)::Float64
    abs2(getx(a2) - getx(a1)) + abs2(gety(a2) - gety(a1))
end

# intersection(Point2D(0., 0.),Point2D(2., 0.),Point2D(0., -1.),Point2D(1., 2.))

function clip(a1::Point2D, a2::Point2D)::Point2D
    res = intersection(a1, a2, Point2D(1., 1.), Point2D(2., 1.))
    (res != nothing) && return res
    res = intersection(a1, a2, Point2D(1., 1.), Point2D(1., 2.))
    (res != nothing) && return res
    res = intersection(a1, a2, Point2D(2., 2.), Point2D(1., 2.))
    (res != nothing) && return res
    res = intersection(a1, a2, Point2D(2., 2.), Point2D(2., 1.))
    (res != nothing) && return res
    a2
end
# clip([1.5,1.5], [1.9, 1.9])
# clip([1.5,1.5], [2.9, 1.9])
# clip([1.5,1.5], [2.9, 2.9])
# clip([1.5,1.5], [1.9, 2.9])

## clips the voronoi cells to the domain
function clipcells(tess, pts)
    ### clipped segments
    segclipd = Dict{Tuple{Point2D,Point2D}, Any}()
    for ed in voronoiedges(tess)
        pta, ptb = geta(ed), getb(ed)
        pax, pay = getx(pta), gety(pta)
        aok = (min_coord <= pax <= max_coord) && (min_coord <= pay <= max_coord)
        pbx, pby = getx(ptb), gety(ptb)
        bok = (min_coord <= pbx <= max_coord) && (min_coord <= pby <= max_coord)

        if aok && !bok  # a in, b out
            segclipd[(pta, ptb)] = clip(pta, ptb)
        elseif !aok && bok # a out, b in
            segclipd[(ptb, pta)] = clip(ptb, pta)
        elseif !aok && !bok # a out and b out but they may intersect in 2 points
            ints = []
            res = intersection(pta, ptb, Point2D(1., 1.), Point2D(2., 1.))
            (res != nothing) && push!(ints, res)
            res = intersection(pta, ptb, Point2D(1., 1.), Point2D(1., 2.))
            (res != nothing) && push!(ints, res)
            res = intersection(pta, ptb, Point2D(2., 2.), Point2D(1., 2.))
            (res != nothing) && push!(ints, res)
            res = intersection(pta, ptb, Point2D(2., 2.), Point2D(2., 1.))
            (res != nothing) && push!(ints, res)

            if length(ints) == 0
                segclipd[(pta, ptb)] = segclipd[(ptb, pta)] = nothing
            elseif length(ints) == 2
                d1, d2 = norm2(pta,ints[1]), norm2(pta, ints[2])
                if d1 <= d2
                    segclipd[(pta, ptb)] = ints
                else
                    segclipd[(pta, ptb)] = reverse(ints)
                end
            else
                error("inconsistent number of intersections")
            end

        end
    end

    ### cell segment by particle
    foreach(pt -> pt.walls = [], pts)
    for ed in voronoiedges(tess)
        push!(getgena(ed).walls, CellWall(geta(ed), getb(ed), getgenb(ed)))
        push!(getgenb(ed).walls, CellWall(geta(ed), getb(ed), getgena(ed)))
    end

    ## order cell segment in a cycle
    # pt = pts[6]
    for pt in pts
        segs = pt.walls
        ui = [2:length(segs);]
        osegs = CellWall[ segs[1] ]
        eend = segs[1].b
        while length(ui) > 0
            noi = findfirst( [segs[ui[i]].a === eend for i in 1:length(ui)] )
            if noi != nothing
                nseg = segs[ui[noi]]
                push!(osegs, nseg)
                eend = nseg.b
            else
                noi = findfirst( [segs[ui[i]].b === eend for i in 1:length(ui)] )
                (noi == nothing) && error("not a complete cycle")
                nseg = segs[ui[noi]]
                push!(osegs, CellWall(nseg.b, nseg.a, nseg.ogen))
                eend = nseg.a
            end
            splice!(ui, noi)
        end
        pt.walls = osegs
    end

    ## clip cell segment
    for pt in pts
        segs = pt.walls
        csegs = CellWall[]
        # set lastpt if cycle starts with a point outside
        lastpt = nothing
        if haskey(segclipd, (segs[1].b, segs[1].a))
            if haskey(segclipd, (segs[end].a, segs[end].b))
                lastpt = segclipd[(segs[end].a, segs[end].b)]
            end
        end

        # clip cell walls
        for seg in segs
            if haskey(segclipd, (seg.a, seg.b))
                cpt = segclipd[(seg.a, seg.b)]
                if cpt == nothing
                elseif isa(cpt, Vector)
                    push!(csegs, CellWall(cpt[1], cpt[2], nothing))
                    lastpt = cpt[2]
                else
                    push!(csegs, CellWall(seg.a, cpt, seg.ogen))
                    lastpt = cpt
                end
            elseif haskey(segclipd, (seg.b, seg.a))
                cpt = segclipd[(seg.b, seg.a)]
                if cpt == nothing
                elseif isa(cpt, Vector)
                    push!(csegs, CellWall(cpt[2], cpt[1], nothing))
                    lastpt = cpt[1]
                else
                    (lastpt != nothing) && push!(csegs, CellWall(lastpt, cpt, nothing))
                    push!(csegs, CellWall(cpt, seg.b, seg.ogen))
                    lastpt = cpt
                end
            else
                push!(csegs, seg)
            end
        end
        pt.walls = csegs
    end
end

# sort(collect( [ (getgena(ed)._x, getgenb(ed)._y) for ed in voronoiedges((tess)) ] ))
# length(pts)

κ = 10.  # rounding force
ϵ = 0.01 # inner border

function oneloop(pts, dt)
    tess = DelaunayTessellation2D{Particle}(100)
    push!(tess, pts)

    ## prepare the cell around each particle
    clipcells(tess, pts)

    ### surface and perimeter calculation to evaluate P
    pt = pts[1];
    for pt in pts
        pos = Point2D(getx(pt), gety(pt))
        pt.surface = 0.
        pt.perimeter = 0.

        prevdx, prevdy = getx(pt.walls[1].a) - getx(pos), gety(pt.walls[1].a) - gety(pos)
        for ptc in pt.walls
            curdx, curdy = getx(ptc.b) - getx(pos), gety(ptc.b) - gety(pos)
            pt.surface   += abs(prevdx*curdy - curdx*prevdy) / 2
            pt.perimeter += sqrt(norm2(ptc.a, ptc.b))
            prevdx, prevdy = curdx, curdy
        end
        curdx, curdy = getx(pt.walls[1].a) - getx(pos), gety(pt.walls[1].a) - gety(pos)
        pt.surface   += abs(prevdx*curdy - curdx*prevdy) / 2
        pt.perimeter += sqrt(norm2(pt.walls[1].a, pt.walls[end].b))

        pt.pressure = (pt.mass / pt.surface) ^ γ # iso temp gas ?
    end
    maxp = maximum( pt.pressure for pt in pts)

    ### calculate acceleration
    # pt=pts[6]
    for pt in pts
        pti = Point2D(getx(pt), gety(pt))
        pt.ax = pt.ay = 0.

        # ptc = first(pt.walls)
        for ptc in pt.walls
            (ptc.ogen==nothing) && continue

            ptj = Point2D(getx(ptc.ogen), gety(ptc.ogen))
            pa, pb = ptc.a, ptc.b
            Aij = sqrt(norm2(pa, pb))

            Rij = sqrt(norm2(pti, ptj))
            eijx,   eijy   = (ptj._x - pti._x) / Rij, (ptj._y - pti._y) / Rij
            midAx,  midAy  = (pa._x + pb._x) / 2,     (pa._y + pb._y) / 2
            midijx, midijy = (ptj._x + pti._x) / 2,   (ptj._y + pti._y) / 2
            cijx,   cijy   = midAx - midijx,           midAy - midijy

            Pi, Pj = pt.pressure, ptc.ogen.pressure

            fx = -Aij * ( (Pi + Pj) * eijx / 2. + (Pj - Pi) * cijx / Rij )
            fy = -Aij * ( (Pi + Pj) * eijy / 2. + (Pj - Pi) * cijy / Rij )

            # PPO regularizing term
            Pij = 4 * pt.mass / ( Rij + pt.perimeter )
            Pji = 4 * pt.mass / ( Rij + ptc.ogen.perimeter )

            force = -Aij * κ * ( Pij - Pi + Pji - Pj) / 2.
            fx -= force * eijx
            fy -= force * eijy

            pt.ax += fx / pt.mass
            pt.ay += fy / pt.mass
        end
    end
    maxa = maximum( sqrt(pt.ax*pt.ax + pt.ay*pt.ay) for pt in pts)

    dt2 = min(dt, 0.01 / maxa)
    for pt in pts
        pt.vx += dt2 * pt.ax
        pt._x += dt2 * pt.vx
        if pt._x > max_coord - ϵ
           pt._x = max_coord - ϵ
           pt.vx = max(pt.vx, 0.)
        elseif pt._x < min_coord + ϵ
           pt._x = min_coord + ϵ
           pt.vx = min(pt.vx, 0.)
        end

        pt.vy += dt2 * pt.ay
        pt._y += dt2 * pt.vy
        pt._y = max(min_coord + ϵ, min(max_coord - ϵ, pt._y))
        if pt._y > max_coord - ϵ
           pt._y = max_coord - ϵ
           pt.vy = max(pt.vy, 0.)
        elseif pt._y < min_coord + ϵ
           pt._y = min_coord + ϵ
           pt.vy = min(pt.vy, 0.)
        end
    end
    maxv = maximum( sqrt(pt.vx*pt.vx + pt.vy*pt.vy) for pt in pts)

    (maxpress = round(maxp),              maxacc = round(maxa, sigdigits=3),
     maxvit   = round(maxv, sigdigits=3), dt     = round(dt2, sigdigits=3)), tess
end


###############################################################################
### plotting
###############################################################################

function plpoints(tess, pts)
    pd = [(x=getx(pt), y=gety(pt),
           id = pt.index,
           pressure  = round(pt.pressure),
           surface   = round(pt.surface, digits=4),
           perimeter = round(pt.perimeter, digits=3)) for pt in pts]

    pd2 = NamedTuple{(:x,:y,:idx), Tuple{Float64, Float64, Int64}}[]
    for (idx, ed) in enumerate(voronoiedges(tess))
        push!(pd2, (x=getx(geta(ed)), y=gety(geta(ed)), idx=idx))
        push!(pd2, (x=getx(getb(ed)), y=gety(getb(ed)), idx=idx))
    end

    @vlplot(background=:lightgrey, width=400, height=400) +
    @vlplot(data= pd,
            x={:x, typ=:quantitative, scale={domain=[1.,2.]}},
            y={:y, typ=:quantitative, scale={domain=[1.,2.]}},
            color={value=:black}, size={value=3},
            tooltip= [ {field=:pressure}, {field=:surface}, {field=:perimeter} ],
            mark=:point) +
    @vlplot(data= pd,
                    x={:x, typ=:quantitative, scale={domain=[1.,2.]}},
                    y={:y, typ=:quantitative, scale={domain=[1.,2.]}},
                    text="id:q",
                    tooltip= [ {field=:pressure}, {field=:surface}, {field=:perimeter} ],
                    mark={:text, dx=15, dy=4}) +
        @vlplot(data=pd2,
            x={:x, typ=:quantitative, scale={domain=[1.,2.]}},
            y={:y, typ=:quantitative, scale={domain=[1.,2.]}},
            detail="idx:n", mark={:line, clip=true})

end

function plpoints2(tess, pts)
    pd = [(x = getx(pt), y = gety(pt),
           id = pt.index,
           pressure  = round(pt.pressure),
           surface   = round(pt.surface, digits=4),
           perimeter = round(pt.perimeter, digits=3)) for pt in pts]

   maxa = maximum( sqrt(pt.ax*pt.ax + pt.ay*pt.ay) for pt in pts )
   pd3 = [(x = getx(pt), y = gety(pt), id = pt.index) for pt in pts]
   pd3 = vcat(pd3, [(x = getx(pt) + 0.1 * pt.ax / maxa,
                     y = gety(pt) + 0.1 * pt.ay / maxa, id = pt.index) for pt in pts])




    pd2 = NamedTuple{(:x,:y,:idx), Tuple{Float64, Float64, Int64}}[]
    for (idx, ed) in enumerate(voronoiedges(tess))
        push!(pd2, (x=getx(geta(ed)), y=gety(geta(ed)), idx=idx))
        push!(pd2, (x=getx(getb(ed)), y=gety(getb(ed)), idx=idx))
    end

    @vlplot(background=:lightgrey, width=400, height=400) +
    @vlplot(data= pd,
            x={:x, typ=:quantitative, scale={domain=[1.,2.]}},
            y={:y, typ=:quantitative, scale={domain=[1.,2.]}},
            color={value=:black}, size={value=3},
            mark=:point) +
    @vlplot(data= pd3,
            x={:x, typ=:quantitative, scale={domain=[1.,2.]}},
            y={:y, typ=:quantitative, scale={domain=[1.,2.]}},
            detail="id:n", color={value=:red},
            mark={:line}) +
    @vlplot(data=pd2,
        x={:x, typ=:quantitative, scale={domain=[1.,2.]}},
        y={:y, typ=:quantitative, scale={domain=[1.,2.]}},
        detail="idx:n", mark={:line, clip=true})
end

function plotdomains(tess, pts)

    pd = NamedTuple{(:x,:y,:id,:o), Tuple{Float64, Float64, Int64, Int64}}[]
    for (ipt, pt) in enumerate(pts)
        for (rk, cpt) in enumerate(pt.walls)
            push!(pd, (x = getx(cpt.a), y = gety(cpt.a), id = pt.index, o = rk))
        end
        push!(pd, (x = getx(pt.walls[end].b), y = gety(pt.walls[end].b), id=pt.index, o=length(pt.walls)+1))
    end

    pd2 = NamedTuple{(:x,:y,:idx), Tuple{Float64, Float64, Int64}}[]
    for (idx, ed) in enumerate(voronoiedges(tess))
        push!(pd2, (x=getx(geta(ed)), y=gety(geta(ed)), idx=idx))
        push!(pd2, (x=getx(getb(ed)), y=gety(getb(ed)), idx=idx))
    end

    pd3 = [(x = getx(pt), y = gety(pt), id=pt.index) for pt in pts]

    @vlplot(background=:lightgrey, width=1000, height=1000) +
        @vlplot(data= pd,
            x={:x, typ=:quantitative, scale={domain=[0.5,2.5]}},
            y={:y, typ=:quantitative, scale={domain=[0.5,2.5]}},
            color="id:n", order="o:n", mark={:area, interpolate="linear-closed"}, detail="id:n") +
        @vlplot(data= pd3,
            x={:x, typ=:quantitative, scale={domain=[1.,2.]}},
            y={:y, typ=:quantitative, scale={domain=[1.,2.]}},
            color={value=:black}, size={value=3},
            mark=:point) +
        @vlplot(data=pd2,
            x={:x, typ=:quantitative, scale={domain=[1.,2.]}},
            y={:y, typ=:quantitative, scale={domain=[1.,2.]}},
            detail="idx:n", mark={:line, clip=true}) +
        @vlplot(data= pd3,
            x={:x, typ=:quantitative, scale={domain=[1.,2.]}},
            y={:y, typ=:quantitative, scale={domain=[1.,2.]}},
            text="id:n",
            tooltip= [ {field=:pressure}, {field=:surface}, {field=:perimeter} ],
            mark={:text, dx=10, dy=4})
end
