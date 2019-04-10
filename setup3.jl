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

mutable struct Particle <: ParticleType
    _x::Float64
    _y::Float64
    index::Int64
    mass::Float64
    density::Float64
    entropy::Float64
    pressure::Float64
    surface::Float64
    perimeter::Float64
    vx::Float64
    vy::Float64
    ax::Float64
    ay::Float64
    closest::Float64
end

P = s * ρ ^ γ


function show(io::IO, part::Particle)
    show(io, "#$(part.index) ($(round(part._x, digits=2)),$(round(part._y, digits=2)))" *
        ", mass $(round(part.mass, sigdigits=2))" *
        ", pressure $(round(part.pressure, sigdigits=2))" *
        ", area $(round(part.surface, sigdigits=2))" *
        ", perim $(round(part.perimeter, sigdigits=2))")
end
# pts[4]

Particle(x::Float64, y::Float64)             = Particle(x, y, 0)
Particle(x::Float64, y::Float64, ind::Int64) = Particle(x, y, ind, 1., 1.)
Particle(x::Float64, y::Float64, ind::Int64, mass::Float64, density::Float64) =
    Particle(x, y, ind, mass, density, 0., 0., 0., 0., 0., 0., 0., Inf)
getx(p::Particle) = p._x
gety(p::Particle) = p._y


##########################
# isentropic gas, H₂O γ ~=  1.33
const γ = 1.33

##### part generators   #################################

width = max_coord - min_coord
validrange = (min= min_coord + width / 3.0,
              max   = max_coord - width / 3.0)
validwidth = width / 3.
function part_gen1(Npart)
    [Particle(validrange.min + validwidth * rand(),
              validrange.min + validwidth * rand(),
              i, 1. / Npart, 1., NaN, NaN, NaN,
              0., 0., 0., 0., []) for i in 1:Npart]

end

function part_gen2(Npart)
    coordrg = min_coord+0.05:sqrt(1/Npart):max_coord-0.05
    vec( [ Particle(x, y + rand()*0.001,
                    i, 1. /Npart, 1., NaN, NaN, NaN,
                    0., 0., 0., 0., []) for
                    (i, (x,y)) in enumerate(product(coordrg, coordrg)) ] )
end

function part_gen3(Npart)
    unitarea     = validwidth * validwidth * 0.5 / Npart
    unittriangle = sqrt(unitarea * 4 / sqrt(3))
    htri = unittriangle * sqrt(3) / 2
    pts = Particle[]
    for t in 1:Npart
        l = t * unittriangle
        re, di = rem(l, validwidth), div(l , validwidth)
        x = validrange.min + 0.01 + re
        y = validrange.min + 0.01 + di * htri
        (x > validrange.max - 0.01) && continue
        (y > validrange.max - 0.01) && continue
        push!(pts,
            Particle( x, y,
                      t, 1. / Npart, 1., NaN, NaN, NaN,
                      0., 0., 0., 0., []) )
    end
    pts
end

function part_gen4(Npart)
    unitlength = sqrt(sqrt(3) * validwidth * validwidth / Npart)
    uls3 = sqrt(3)*unitlength
    r1, r2 = validrange.min+2ϵ, validrange.max-2ϵ
    pts = Particle[]
    idx = 1
    for y in r1:uls3:r2
        for x in r1:unitlength:r2
            push!(pts, Particle( x, y, idx, 1. / Npart, 1.))
            idx += 1
        end
    end
    for y in r1+uls3/2:uls3:r2
        for x in r1+unitlength/2:unitlength:r2
            push!(pts, Particle( x, y, idx, 1. / Npart, 1.))
            idx += 1
        end
    end
    pts
end


function calc_mass(pts)
    tess = borderize(pts)
    # tess = DelaunayTessellation2D{Particle}(length(pts2))
    # push!(tess, pts2)

    ### surface and perimeter calculation to evaluate P
    foreach(pt -> pt.surface = 0., pts) # only on original points
    for ve in voronoiedges(tess)
        pa, pb = getgena(ve), getgenb(ve)
        (pa.index <= 0) && (pb.index <= 0) && continue

        p1, p2 = geta(ve), getb(ve)

        pa1x, pa1y = getx(p1) - getx(pa), gety(p1) - gety(pa)
        pa2x, pa2y = getx(p2) - getx(pa), gety(p2) - gety(pa)

        surf = abs(pa1x*pa2y - pa2x*pa1y) / 2

        if pa.index > 0
            pa.surface   += surf
        end
        if pb.index > 0
            pb.surface   += surf
        end
    end
    foreach(pt -> pt.mass = pt.density * pt.surface / (validwidth * validwidth), pts)
end

##################################################################################
function norm2(a1::T, a2::T)::Float64 where T <: AbstractPoint2D
    abs2(getx(a2) - getx(a1)) + abs2(gety(a2) - gety(a1))
end


# create ghost points to have a neat border around regular points
function borderize(pts)
    tess = DelaunayTessellation2D{Particle}(2 * length(pts))
    push!(tess, pts)

    # identify points close to the border
    # pts
    ptsflip = falses(maximum( pt.index for pt in pts ), 4)
    # ttt = collect(voronoiedges(tess))
    # ve = ttt[16]
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

function oneloop(pts, dt)
    # pts2 = borderize(pts)
    #
    tess = borderize(pts)

    ves = collect(voronoiedges(tess))

    ### surface and perimeter calculation to evaluate P
    # ttt = filter(ed -> getgenb(ed).index == 1, collect(voronoiedges(tess)))
    # ve = ttt[1]
    for pt in pts
        pt.surface = pt.perimeter = 0.
        pt.closest = Inf
    end

    for ve in ves
        pa, pb = getgena(ve), getgenb(ve)
        (pa.index <= 0) && (pb.index <= 0) && continue

        p1, p2 = geta(ve), getb(ve)

        pa1x, pa1y = getx(p1) - getx(pa), gety(p1) - gety(pa)
        pa2x, pa2y = getx(p2) - getx(pa), gety(p2) - gety(pa)

        surf = abs(pa1x*pa2y - pa2x*pa1y) / 2
        peri = sqrt(norm2(p1, p2))
        dist = sqrt(norm2(pa, pb))

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
        pt.pressure = (pt.mass / pt.surface) ^ γ / pt.density
    end
    maxp = maximum( pt.pressure for pt in pts)
    minR = minimum( pt.closest for pt in pts)

    ### calculate acceleration
    for pt in pts
        pt.ax = pt.ay = 0.
    end

    for ve in ves
        pa, pb = getgena(ve), getgenb(ve)
        (pa.index <= 0) && (pb.index <= 0) && continue

        p1, p2 = geta(ve), getb(ve)

        Aij = sqrt(norm2(p1, p2))
        midAx,  midAy  = (p1._x + p2._x) / 2, (p1._y + p2._y) / 2
        midijx, midijy = (pa._x + pb._x) / 2, (pa._y + pb._y) / 2

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

        Rij = sqrt(norm2(pa, pb))
        eijx,   eijy   = (pb._x - pa._x) / Rij, (pb._y - pa._y) / Rij
        cijx,   cijy   = midAx - midijx,           midAy - midijy

        fx = -Aij * ( (Pa + Pb) * eijx / 2. + (Pb - Pa) * cijx / Rij )
        fy = -Aij * ( (Pa + Pb) * eijy / 2. + (Pb - Pa) * cijy / Rij )

        # PPO regularizing term
        Pab = (4 * pa.mass / ( Rij * Pera )) ^ γ
        Pba = (4 * pb.mass / ( Rij * Perb )) ^ γ

        force = -Aij * κ * ( Pab - Pa + Pba - Pb) / 2.
        fx += force * eijx
        fy += force * eijy

        if pa.index > 0
            pa.ax += fx
            pa.ay += fy
        end

        if pb.index > 0
            pb.ax -= fx
            pb.ay -= fy
        end
    end

    for pt in pts
        pt.ax /= pt.mass
        pt.ay = pt.ay / pt.mass - 1.
    end
    maxa = maximum( sqrt(pt.ax*pt.ax + pt.ay*pt.ay) for pt in pts)

    # calculate time step so as to avoid large moves
    dt2 = dt
    for pt in pts
        nvx, nvy = pt.vx + dt2 * pt.ax, pt.vy + dt2 * pt.ay
        while (sqrt(nvx*nvx+nvy*nvy) * dt2) > 0.1 * pt.closest
            dt2 /= 2.
            nvx, nvy = pt.vx + dt2 * pt.ax, pt.vy + dt2 * pt.ay
        end
    end

    for pt in pts
        pt.vx += dt2 * pt.ax
        pt._x += dt2 * pt.vx
        if pt._x > validrange.max - ϵ
           pt._x = validrange.max - ϵ
           pt.vx = max(pt.vx, 0.)
        elseif pt._x < validrange.min + ϵ
           pt._x = validrange.min + ϵ
           pt.vx = min(pt.vx, 0.)
        end

        pt.vy += dt2 * pt.ay
        pt._y += dt2 * pt.vy
        pt._y = max(validrange.min + ϵ, min(validrange.max - ϵ, pt._y))
        if pt._y > validrange.max - ϵ
           pt._y = validrange.max - ϵ
           pt.vy = max(pt.vy, 0.)
        elseif pt._y < validrange.min + ϵ
           pt._y = validrange.min + ϵ
           pt.vy = min(pt.vy, 0.)
        end
    end
    maxv = maximum( sqrt(pt.vx*pt.vx + pt.vy*pt.vy) for pt in pts)

    (maxpress = round(maxp),              maxacc = round(maxa, sigdigits=3),
     maxvit   = round(maxv, sigdigits=3), minr   = minR,
     dt     = round(dt2, sigdigits=3)), tess
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

function plpoints2(tess, pts, pltrange =[validrange.min-0.01,validrange.max+0.01] )
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

function plpoints2v(tess, pts, pltrange =[validrange.min-0.01,validrange.max+0.01] )
    pd = [(x = getx(pt), y = gety(pt),
           id = pt.index,
           pressure  = round(pt.pressure),
           surface   = round(pt.surface, digits=4),
           perimeter = round(pt.perimeter, digits=3)) for pt in pts]

    maxv = maximum( sqrt(pt.vx*pt.vx + pt.vy*pt.vy) for pt in pts )
    pd3 = [(x = getx(pt), y = gety(pt), id = pt.index) for pt in pts]
    pd3 = vcat(pd3, [(x = getx(pt) + 0.02 * pt.vx / maxv,
                     y = gety(pt) + 0.02 * pt.vy / maxv, id = pt.index) for pt in pts])

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


function plpoints3(pts, pltrange =[validrange.min-0.01,validrange.max+0.01])
    tess = borderize(pts)

    pd = NamedTuple{(:x,:y, :idx, :o, :dens), Tuple{Float64, Float64, Int64, Int64, Float64}}[]
    for (idx, ed) in enumerate(voronoiedges(tess))
        pa, pb = getgena(ed), getgenb(ed)
        p1x, p1y = getx(geta(ed)), gety(geta(ed))
        p2x, p2y = getx(getb(ed)), gety(getb(ed))
        if pa.index > 0
            pax, pay = getx(pa), gety(pa)
            t1 = (x=p1x, y=p1y, idx=idx, o=1, dens=pa.density)
            t2 = (x=p2x, y=p2y, idx=idx, o=2, dens=pa.density)
            t3 = (x=pax, y=pay, idx=idx, o=3, dens=pa.density)
            push!(pd, [t1, t2, t3]...)
        end
        if pb.index > 0
            pbx, pby = getx(pb), gety(pb)
            t1 = (x=p1x, y=p1y, idx=-idx, o=1, dens=pb.density)
            t2 = (x=p2x, y=p2y, idx=-idx, o=2, dens=pb.density)
            t3 = (x=pbx, y=pby, idx=-idx, o=3, dens=pb.density)
            push!(pd, [t1, t2, t3]...)
        end
    end

    @vlplot(background=:lightgrey, width=400, height=400) +
    @vlplot(data= pd,
        x={:x, typ=:quantitative, scale={domain=pltrange}},
        y={:y, typ=:quantitative, scale={domain=pltrange}},
        color="dens:q", order="o:o", detail="idx:n",
        mark={:area, interpolate="linear-closed"})

end

function plpoints3p(pts, pltrange =[validrange.min-0.01,validrange.max+0.01])
    tess = borderize(pts)

    pd = NamedTuple{(:x,:y, :idx, :o, :val), Tuple{Float64, Float64, Int64, Int64, Float64}}[]
    for (idx, ed) in enumerate(voronoiedges(tess))
        pa, pb = getgena(ed), getgenb(ed)
        p1x, p1y = getx(geta(ed)), gety(geta(ed))
        p2x, p2y = getx(getb(ed)), gety(getb(ed))
        if pa.index > 0
            pax, pay = getx(pa), gety(pa)
            t1 = (x=p1x, y=p1y, idx=idx, o=1, val=pa.pressure)
            t2 = (x=p2x, y=p2y, idx=idx, o=2, val=pa.pressure)
            t3 = (x=pax, y=pay, idx=idx, o=3, val=pa.pressure)
            push!(pd, [t1, t2, t3]...)
        end
        if pb.index > 0
            pbx, pby = getx(pb), gety(pb)
            t1 = (x=p1x, y=p1y, idx=-idx, o=1, val=pb.pressure)
            t2 = (x=p2x, y=p2y, idx=-idx, o=2, val=pb.pressure)
            t3 = (x=pbx, y=pby, idx=-idx, o=3, val=pb.pressure)
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
