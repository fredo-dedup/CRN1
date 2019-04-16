using VoronoiDelaunay
using VoronoiDelaunay.GeometricalPredicates
using VegaLite
using ElectronDisplay
using LinearAlgebra
using SparseArrays

using Base.Iterators: product

import VoronoiDelaunay.GeometricalPredicates: getx, gety
import VoronoiDelaunay: min_coord, max_coord

import Base: show, +, *, -, abs

abstract type ParticleType <: AbstractPoint2D end

mutable struct Particle <: ParticleType
    _x::Float64
    _y::Float64
    index::Int64
    mass::Float64
    density::Float64
    pressure::Float64
    entropy::Float64
    surface::Float64
    perimeter::Float64
    v::Point2D
    a::Point2D
    closest::Float64
end

function show(io::IO, part::Particle)
    show(io, "#$(part.index) ($(round(part._x, digits=2)),$(round(part._y, digits=2)))" *
        ", mass=$(round(part.mass, sigdigits=2))" *
        ", press.=$(round(part.pressure, sigdigits=2))" *
        ", surf.=$(round(part.surface, sigdigits=2))" *
        ", perim=$(round(part.perimeter, sigdigits=2))" *
        ", entropy=$(round(part.entropy, sigdigits=2))")
end
# pts[4]

Particle(x::Float64, y::Float64)             = Particle(x, y, 0)
Particle(x::Float64, y::Float64, ind::Int64) = Particle(x, y, ind, 1.)
Particle(x::Float64, y::Float64, ind::Int64, mass::Float64) =
    Particle(x, y, ind, mass, 0., 0., 0., 0., 0., Point2D(0., 0.), Point2D(0., 0.), 0.)
getx(p::Particle) = p._x
gety(p::Particle) = p._y

+(a::AbstractPoint2D, b::AbstractPoint2D) = Point2D(getx(a) + getx(b), gety(a) + gety(b))
-(a::AbstractPoint2D, b::AbstractPoint2D) = Point2D(getx(a) - getx(b), gety(a) - gety(b))
*(α::Float64, b::AbstractPoint2D) = Point2D(α * getx(b), α * gety(b))
∨(a::AbstractPoint2D, b::AbstractPoint2D) = getx(a) * gety(b) - getx(b) * gety(a)
⋄(a::AbstractPoint2D, b::AbstractPoint2D) = getx(a) * getx(b) + gety(a) * gety(b)
abs(a::AbstractPoint2D) = sqrt(a ⋄ a)


##########################
# isentropic gas, H₂O γ ~=  1.33
# const γ = 1.33
γ = 1.3


##################################################################################

function norm2(a1::T, a2::T)::Float64 where T <: AbstractPoint2D
    abs2(getx(a2) - getx(a1)) + abs2(gety(a2) - gety(a1))
end



###############################################################################
###   simulation init
###############################################################################

width = max_coord - min_coord
validrange = (min= min_coord + width / 3.0,
              max   = max_coord - width / 3.0)
validwidth = width / 3.
function part_gen1(Npart)
    [Particle(validrange.min + validwidth * rand(),
              validrange.min + validwidth * rand(),
              i) for i in 1:Npart]

end

function part_gen2(Npart)
    coordrg = min_coord+0.05:sqrt(1/Npart):max_coord-0.05
    vec( [ Particle(x, y + rand()*0.001,
                    i) for (i, (x,y)) in enumerate(product(coordrg, coordrg)) ] )
end

function part_gen4(Npart)
    unitlength = sqrt(sqrt(3) * validwidth * validwidth / Npart)
    uls3 = sqrt(3)*unitlength
    r1, r2 = validrange.min+2ϵ, validrange.max-2ϵ
    pts = Particle[]
    idx = 1
    for y in r1:uls3:r2
        for x in r1:unitlength:r2
            push!(pts, Particle( x, y, idx, 1. / Npart))
            idx += 1
        end
    end
    for y in r1+uls3/2:uls3:r2
        for x in r1+unitlength/2:unitlength:r2
            push!(pts, Particle( x, y, idx, 1. / Npart))
            idx += 1
        end
    end
    pts
end

# set entropy such that all pressure are equal
function equalize_pressure(pts, Pref)
    tess = borderize(pts)

    ### surface calculation
    foreach(pt -> pt.surface = 0., pts) # only on original points
    for ve in voronoiedges(tess)
        pa, pb = getgena(ve), getgenb(ve)
        (pa.index <= 0) && (pb.index <= 0) && continue

        p1, p2 = geta(ve), getb(ve)

        pa1, pa2 = p1 - pa, p2 - pa
        surf = abs( pa1 ∨ pa2 ) / 2

        if pa.index > 0
            pa.surface   += surf
        end
        if pb.index > 0
            pb.surface   += surf
        end
    end

    for pt in pts
        ρ = pt.mass / pt.surface
        pt.entropy = Pref * ρ ^ (-γ)
    end

end



# set Particle mass using density and current surface
function calc_mass(pts, Pref)
    tess = borderize(pts)

    ### surface and perimeter calculation to evaluate P
    foreach(pt -> pt.surface = 0., pts) # only on original points
    for ve in voronoiedges(tess)
        pa, pb = getgena(ve), getgenb(ve)
        (pa.index <= 0) && (pb.index <= 0) && continue

        p1, p2 = geta(ve), getb(ve)

        pa1, pa2 = p1 - pa, p2 - pa
        surf = abs( pa1 ∨ pa2 ) / 2

        if pa.index > 0
            pa.surface   += surf
        end
        if pb.index > 0
            pb.surface   += surf
        end
    end
    for pt in pts
        pt.mass = pt.density * pt.surface / (validwidth * validwidth)
        pt.entropy = Pref * pt.density ^ (-γ)
    end
end


# create ghost points to have a neat square border
function borderize(pts)
    tess = DelaunayTessellation2D{Particle}(2 * length(pts))
    push!(tess, pts)

    # identify points close to the border
    ptsflip = falses(maximum( pt.index for pt in pts ), 4)
    for ve in voronoiedges(tess)
        up    = (ve._a._y > validrange.max) || (ve._b._y > validrange.max)
        down  = (ve._a._y < validrange.min) || (ve._b._y < validrange.min)
        right = (ve._a._x > validrange.max) || (ve._b._x > validrange.max)
        left  = (ve._a._x < validrange.min) || (ve._b._x < validrange.min)

        if getgena(ve).index != 0
            up    && (ptsflip[getgena(ve).index, 1] = true)
            down  && (ptsflip[getgena(ve).index, 2] = true)
            right && (ptsflip[getgena(ve).index, 3] = true)
            left  && (ptsflip[getgena(ve).index, 4] = true)
        end
        if getgenb(ve).index != 0
            up    && (ptsflip[getgenb(ve).index, 1] = true)
            down  && (ptsflip[getgenb(ve).index, 2] = true)
            right && (ptsflip[getgenb(ve).index, 3] = true)
            left  && (ptsflip[getgenb(ve).index, 4] = true)
        end
    end

    ghostpts = Particle[]
    gidx = -1
    for pt in pts
        if ptsflip[pt.index,1]  # up
            push!(ghostpts, Particle(pt._x, 2validrange.max - pt._y, gidx))
            gidx -= 1
        end
        if ptsflip[pt.index,2]  # down
            push!(ghostpts, Particle(pt._x, 2validrange.min - pt._y, gidx))
            gidx -= 1
        end
        if ptsflip[pt.index,3]  # right
            push!(ghostpts, Particle(2validrange.max - pt._x, pt._y, gidx))
            gidx -= 1
        end
        if ptsflip[pt.index,4]  # left
            push!(ghostpts, Particle(2validrange.min - pt._x, pt._y, gidx))
            gidx -= 1
        end
    end

    # push!(tess, ghostpts)

    allpts = vcat(pts, ghostpts)
    tess = DelaunayTessellation2D{Particle}(length(allpts))
    push!(tess, allpts)

    tess
end

###############################################################################
### simulation run
###############################################################################

# sort(collect( [ (getgena(ed)._x, getgenb(ed)._y) for ed in voronoiedges((tess)) ] ))
# length(pts)
κ = 0.01  # rounding force
ϵ = 0.001 # inner border
g = Point2D(0., -0.1)

function oneloop(pts, dt)
    tess = borderize(pts)
    ves = collect(voronoiedges(tess))

    ### surface and perimeter calculation to evaluate P
    # ttt = filter(ed -> getgenb(ed).index == 1, collect(voronoiedges(tess)))
    # ve = ttt[1]
    for pt in pts # only on original points, ghost points can be ignored
        pt.surface = pt.perimeter = 0.
        pt.closest = Inf
    end

    for ve in ves
        pa, pb = getgena(ve), getgenb(ve)
        (pa.index <= 0) && (pb.index <= 0) && continue

        p1, p2 = geta(ve), getb(ve)

        pa1, pa2 = p1 - pa, p2 - pa

        surf = abs( pa1 ∨ pa2 ) / 2
        peri = abs(p1 - p2)
        dist = abs(pa - pb)

        if pa.index > 0
            pa.surface   += surf
            pa.perimeter += peri
        end
        if pb.index > 0
            pb.surface   += surf
            pb.perimeter += peri
        end
        if (pa.index > 0) && (pb.index > 0)
            pa.closest = min(pa.closest, dist)
            pb.closest = min(pb.closest, dist)
        end
    end

    for pt in pts
        pt.density = pt.mass / pt.surface
        pt.pressure = pt.entropy * pt.density ^ γ
    end
    # foreach(pt -> pt.pressure = (pt.mass / pt.surface / pt.density) ^ γ, pts)
    maxp = maximum( pt.pressure for pt in pts)
    minR = minimum( pt.closest for pt in pts)

    ### calculate acceleration
    foreach(pt -> pt.a = Point2D(0., 0.), pts) # only on original points
    for ve in ves
        pa, pb = getgena(ve), getgenb(ve)
        (pa.index <= 0) && (pb.index <= 0) && continue

        p1, p2 = geta(ve), getb(ve)

        Aij = abs(p2 - p1)
        midA, midij = 0.5 * (p1 + p2), 0.5 * (pa + pb)

        if pa.index <= 0
            Pa = Pb = pb.pressure
            Pera, Perb = Inf, pb.perimeter
        elseif pb.index <= 0
            Pa = Pb = pa.pressure
            Pera, Perb = pa.perimeter, Inf
        else
            Pa, Pb = pa.pressure, pb.pressure
            Pera, Perb = pa.perimeter, pb.perimeter
        end

        Rij = abs(pb - pa)
        eij = 1. / Rij * (pb - pa)
        cij = midA - midij

        force = -Aij * ( (Pa + Pb) / 2. * eij + (Pb - Pa) / Rij * cij )

        # PPO regularizing term
        Pab = (4 * pa.mass / ( Rij * Pera )) ^ γ
        Pba = (4 * pb.mass / ( Rij * Perb )) ^ γ

        regforce = -Aij * κ * ( Pab - Pa + Pba - Pb) / 2.
        force += regforce * eij

        if pa.index > 0
            pa.a += force
        end
        if pb.index > 0
            pb.a -= force
        end
    end

    for pt in pts
        pt.a = (1. / pt.mass) * pt.a + g
    end
    maxa = maximum( abs(pt.a) for pt in pts)

    # calculate time step so as to avoid large moves
    dt2 = dt
    for pt in pts
        nv = pt.v + dt2 * pt.a
        while (abs(nv) * dt2) > 0.1 * pt.closest
            dt2 /= 2.
            nv = pt.v + dt2 * pt.a
        end
    end

    for pt in pts
        pt.v += dt2 * pt.a
        pt._x += dt2 * pt.v._x
        pt._y += dt2 * pt.v._y

        if pt._x > validrange.max - ϵ
           pt._x = validrange.max - ϵ
           pt.v = Point2D(min(pt.v._x, 0.), pt.v._y)
        elseif pt._x < validrange.min + ϵ
           pt._x = validrange.min + ϵ
           pt.v = Point2D(max(pt.v._x, 0.), pt.v._y)
        end

        if pt._y > validrange.max - ϵ
           pt._y = validrange.max - ϵ
           pt.v = Point2D(pt.v._x, min(pt.v._y, 0.))
        elseif pt._y < validrange.min + ϵ
           pt._y = validrange.min + ϵ
           pt.v = Point2D(pt.v._x, max(pt.v._y, 0.))
        end
    end
    maxv = maximum( abs(pt.v) for pt in pts )

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

function plotspeed(pts, pltrange =[validrange.min-0.01,validrange.max+0.01] )
    tess = borderize(pts)

    pd = [(x = getx(pt), y = gety(pt), id = pt.index) for pt in pts]

    maxv = maximum( abs(pt.v) for pt in pts )
    pd3 = [(x = getx(pt), y = gety(pt), id = pt.index) for pt in pts]
    pd3 = vcat(pd3, [(x = getx(pt) + 0.03 * pt.v._x / maxv,
                      y = gety(pt) + 0.03 * pt.v._y / maxv, id = pt.index) for pt in pts])

    pd2 = NamedTuple{(:x,:y,:idx), Tuple{Float64, Float64, Int64}}[]
    for (idx, ed) in enumerate(voronoiedges(tess))
        pa, pb = getgena(ed), getgenb(ed)
        (pa.index <= 0) && (pb.index <= 0) && continue

        push!(pd2, (x=getx(geta(ed)), y=gety(geta(ed)), idx=idx))
        push!(pd2, (x=getx(getb(ed)), y=gety(getb(ed)), idx=idx))
    end

    # pltrange = [validrange.min-0.01,validrange.max+0.01]

    @vlplot(background=:lightgrey, width=400, height=400) +
    @vlplot(data= pd,
            x={:x, typ=:quantitative, scale={domain=pltrange}},
            y={:y, typ=:quantitative, scale={domain=pltrange}},
            color={value=:black}, size={value=3},
            mark=:point) +
    @vlplot(data= pd3,
            x={:x, typ=:quantitative, scale={domain=pltrange}},
            y={:y, typ=:quantitative, scale={domain=pltrange}},
            detail="id:n", color={value=:red},
            mark={:line}) +
    @vlplot(data=pd2,
        x={:x, typ=:quantitative, scale={domain=pltrange}},
        y={:y, typ=:quantitative, scale={domain=pltrange}},
        detail="idx:n", mark={:line, clip=true})
end


function plotval(pts; pltrange =[validrange.min-0.01,validrange.max+0.01], field=:pressure)
    tess = borderize(pts)

    pd = NamedTuple{(:x,:y, :idx, :o, :val), Tuple{Float64, Float64, Int64, Int64, Float64}}[]
    for (idx, ed) in enumerate(voronoiedges(tess))
        pa, pb = getgena(ed), getgenb(ed)
        p1x, p1y = getx(geta(ed)), gety(geta(ed))
        p2x, p2y = getx(getb(ed)), gety(getb(ed))
        if pa.index > 0
            pax, pay = getx(pa), gety(pa)
            t1 = (x=p1x, y=p1y, idx=idx, o=1, val=getfield(pa, field))
            t2 = (x=p2x, y=p2y, idx=idx, o=2, val=getfield(pa, field))
            t3 = (x=pax, y=pay, idx=idx, o=3, val=getfield(pa, field))
            push!(pd, [t1, t2, t3]...)
        end
        if pb.index > 0
            pbx, pby = getx(pb), gety(pb)
            t1 = (x=p1x, y=p1y, idx=-idx, o=1, val=getfield(pb, field))
            t2 = (x=p2x, y=p2y, idx=-idx, o=2, val=getfield(pb, field))
            t3 = (x=pbx, y=pby, idx=-idx, o=3, val=getfield(pb, field))
            push!(pd, [t1, t2, t3]...)
        end
    end

    @vlplot(background=:lightgrey, width=400, height=400) +
    @vlplot(data= pd,
        x={:x, typ=:quantitative, scale={domain=pltrange}},
        y={:y, typ=:quantitative, scale={domain=pltrange}},
        color="val:q", order="o:o", detail="idx:n",
        mark={:area, interpolate="linear-closed"})
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
