using VoronoiDelaunay
using VoronoiDelaunay.GeometricalPredicates
using VegaLite
using ElectronDisplay
using LinearAlgebra

import VoronoiDelaunay.GeometricalPredicates: getx, gety
import VoronoiDelaunay: min_coord, max_coord

mutable struct MyPointType4 <: VoronoiDelaunay.AbstractPoint2D
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
end

MyPointType4(x::Float64, y::Float64) = MyPointType4(x,y,0)
MyPointType4(x::Float64, y::Float64, ind::Int64) =
    MyPointType4(x, y, ind, 1., 1., 0., 0., 0., 0., 0., 0., 0.)
getx(p::MyPointType4) = p._x
gety(p::MyPointType4) = p._y



######################################
dt = 1e-4
width = max_coord - min_coord
pts = [MyPointType4(min_coord + rand() * width,
                    min_coord + rand() * width, i) for i in 1:10]

# pts = vec( [ MyPointType4(x, y) for x in min_coord+0.05:0.08:max_coord-0.05,
#                 y in min_coord+0.05:0.08:max_coord-0.05] )

foreach(pt -> pt.mass = 100, pts)

tess = DelaunayTessellation2D{MyPointType4}(100)
push!(tess, pts)
eds = collect(voronoiedges(tess))

_, tess = oneloop2(pts, dt)
plpoints(tess, pts)


ppplot = plpoints(tess, pts);

tmppath = VegaLite.writehtml_full(VegaLite.JSON.json(ppplot.params))
show(VegaLite.JSON.json(ppplot.params))
VegaLite.launch_browser(tmppath) # Open the browser


 |> display

for i in 1:100
    global tess
    _, tess = oneloop2(pts, dt)
end
plpoints(tess, pts)













##################################################################################


function intersection(a1, a2, b1, b2)
    # a1, a2, b1, b2 = [0., 0.], [2., 0.], [0., -1.], [1., 4.]
    A = hcat( a2 .- a1 , b1 .- b2 )
    (det(A) == 0.) && return nothing
    b = b1 .- a1
    x = A \ b
    all(0 .<= x .<= 1.) && return a1 + x[1] .* (a2 .- a1)
    nothing
end

function clip(a1, a2)
    res = intersection(a1, a2, [1., 1.], [2., 1.])
    (res != nothing) && return res
    res =  intersection(a1, a2, [1., 1.], [1., 2.])
    (res != nothing) && return res
    res =  intersection(a1, a2, [2., 2.], [1., 2.])
    (res != nothing) && return res
    res =  intersection(a1, a2, [2., 2.], [2., 1.])
    (res != nothing) && return res
    a2
end


clip([1.5,1.5], [1.9, 1.9])
clip([1.5,1.5], [2.9, 1.9])
clip([1.5,1.5], [2.9, 2.9])
clip([1.5,1.5], [1.9, 2.9])


κ = 0.0
function oneloop2(pts, dt)
    tess = DelaunayTessellation2D{MyPointType4}(100)
    push!(tess, pts)

    ### surface and perimeter calculation to evaluate pressure
    foreach(pt -> pt.surface = 0., pts)
    foreach(pt -> pt.perimeter = 0., pts)
    # ed = first(voronoiedges(tess))
    # eds = collect(voronoiedges(tess))
    # ed = eds[1]
    for ed in voronoiedges(tess)
        gens = MyPointType4[]
        pax, pay = getx(getgena(ed)), gety(getgena(ed))
        if getgena(ed).index > 0  # not a fake point
            push!(gens, getgena(ed))
            (getgenb(ed).index > 0) && push!(gens, getgenb(ed))  # other gen is not fake too
        elseif getgenb(ed).index > 0  # not a fake point
            push!(gens, getgenb(ed))
        else
            error("both fakes")
        end

        pax, pay = getx(gens[1]), gety(gens[1])

        p1x, p1y = getx(geta(ed)), gety(geta(ed))
        aout = (min_coord <= p1x <= max_coord) &&
               (min_coord <= p1y <= max_coord)
        if !aout
            p1x, p1y = clip([pax, pay], [p1x, p1y])
        end

        p2x, p2y = getx(getb(ed)), gety(getb(ed))
        bout = (min_coord <= p2x <= max_coord) &&
               (min_coord <= p2y <= max_coord)
        if !bout
            p2x, p2y = clip([pax, pay], [p1x, p1y])
        end

        d1x = p1x - pax
        d1y = p1y - pay
        d2x = p2x - pax
        d2y = p2y - pay
        d12x = p2x - p1x
        d12y = p2y - p1y

        S = abs(d1x*d2y - d2x*d1y) / 2
        foreach(gpt -> gpt.surface += S, gens)

        L = sqrt(d12x*d12x + d12y*d12y)
        foreach(gpt -> gpt.perimeter += L, gens)
    end
    foreach(pt -> pt.pressure = pt.density / pt.surface, pts)
    maxp = maximum( pt.pressure for pt in pts)

    # extrema( pt.perimeter for pt in pts)

    ### calculate acceleration
    foreach(pt -> pt.ax = pt.ay = 0., pts)
    eds = collect(voronoiedges(tess))
    # ed = eds[1]
    for ed in voronoiedges(tess)
        ((getgena(ed).index == 0) || (getgenb(ed).index == 0)) && continue
        Aij = sqrt(abs2(getx(geta(ed))-getx(getb(ed))) +
                   abs2(gety(geta(ed))-gety(getb(ed))))
        Rij = sqrt(abs2(getx(getgena(ed))-getx(getgenb(ed))) +
                   abs2(gety(getgena(ed))-gety(getgenb(ed))))
        eijx = (getx(getgena(ed)) - getx(getgenb(ed))) / Rij
        eijy = (gety(getgena(ed)) - gety(getgenb(ed))) / Rij
        midAx = (getx(geta(ed)) + getx(getb(ed))) / 2
        midAy = (gety(geta(ed)) + gety(getb(ed))) / 2
        midijx = (getx(getgena(ed)) + getx(getgenb(ed))) / 2
        midijy = (gety(getgena(ed)) + gety(getgenb(ed))) / 2

        cijx, cijy = midAx - midijx, midAy - midijy
        Pi, Pj = getgena(ed).pressure, getgenb(ed).pressure

        fx = -Aij * ( (Pi + Pj) * eijx / 2. + (Pj - Pi) * cijx / Rij )
        fy = -Aij * ( (Pi + Pj) * eijy / 2. + (Pj - Pi) * cijy / Rij )

        # PPO regularizing term
        Pij = 4 * getgena(ed).mass / ( Rij + getgena(ed).perimeter )
        Pji = 4 * getgenb(ed).mass / ( Rij + getgenb(ed).perimeter )

        force = -Aij * κ * ( Pij - Pi + Pji - Pj) / 2.
        fx -= force * eijx
        fy -= force * eijy

        getgena(ed).ax += fx / getgena(ed).mass
        getgena(ed).ay += fy / getgena(ed).mass
        getgenb(ed).ax -= fx / getgenb(ed).mass
        getgenb(ed).ay -= fy / getgenb(ed).mass
    end
    maxa = maximum( sqrt(pt.ax*pt.ax + pt.ay*pt.ay) for pt in pts)

    dt2 = dt * min(1.0, 0.01 / maxa)
    for pt in pts
        pt.vx += dt2 * pt.ax
        pt._x += dt2 * pt.vx
        pt._x = max(min_coord, min(max_coord, pt._x))
        pt.vy += dt2 * pt.ay
        pt._y += dt2 * pt.vy
        pt._y = max(min_coord, min(max_coord, pt._y))
    end
    maxv = maximum( sqrt(pt.vx*pt.vx + pt.vy*pt.vy) for pt in pts)

    (maxpress=maxp, maxacc=maxa, maxvit=maxv, dt=dt2), tess
end

round(1.45465,digits=2)

function plpoints(tess, pts)
    pd = [(x=getx(pt), y=gety(pt),
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
    @vlplot(data=pd2,
            x={:x, typ=:quantitative, scale={domain=[1.,2.]}},
            y={:y, typ=:quantitative, scale={domain=[1.,2.]}},
            detail="idx:n", mark={:line, clip=true})

end
