
include("setup.jl")



######################################
const MPT = MyPointType4
dt = 0.1
width = max_coord - min_coord
pts = [MPT(min_coord + rand() * width, min_coord + rand() * width,
           i, 1., 1., NaN, NaN, NaN,
           0., 0., 0., 0.) for i in 1:10]

coordrg = min_coord+0.05:0.08:max_coord-0.05
pts = vec( [ MPT(min_coord + rand() * width, min_coord + rand() * width,
                 i, 1., 1., NaN, NaN, NaN, 0., 0., 0., 0.) for
                 (i, (x,y)) in enumerate(product(coordrg, coordrg)) ] )

tess = DelaunayTessellation2D{MyPointType4}(100)
push!(tess, pts)
eds = collect(voronoiedges(tess))

_, tess = oneloop2(pts, dt)
plpoints2(tess, pts)
plpoints(tess, pts)


ppplot = plpoints(tess, pts);
tmppath = VegaLite.writehtml_full(VegaLite.JSON.json(ppplot.params))
show(VegaLite.JSON.json(ppplot.params))
VegaLite.launch_browser(tmppath) # Open the browser

t = 0.
for i in 1:1000
    global tess, t
    (mp, ma, mvit, dt2), tess = oneloop2(pts, dt)
    t += dt2
    isfinite(ma) || break
end
t
plpoints2(tess, pts)

##### test   #################################
width = max_coord - min_coord
pts = [Particle(min_coord + rand() * width, min_coord + rand() * width,
           i, 1., 1., NaN, NaN, NaN,
           0., 0., 0., 0., []) for i in 1:6]

tess = DelaunayTessellation2D{Particle}(10)
push!(tess, pts)


eds = collect(voronoiedges(tess))
ied, ed = 12, eds[12]

### clipped segments
segclipd = Dict{Tuple{Point2D,Point2D}, Union{Point2D, Nothing}}()
for (ied, ed) in enumerate(voronoiedges(tess))
    pta, ptb = geta(ed), getb(ed)
    pax, pay = getx(pta), gety(pta)
    aok = (min_coord <= pax <= max_coord) &&
          (min_coord <= pay <= max_coord)
    pbx, pby = getx(ptb), gety(ptb)
    bok = (min_coord <= pbx <= max_coord) &&
          (min_coord <= pby <= max_coord)

    if aok && !bok
        segclipd[(pta, ptb)] = Point2D(clip([pax, pay], [pbx, pby])...)
    elseif !aok && bok
        segclipd[(ptb, pta)] = Point2D(clip([pbx, pby], [pax, pay])...)
    elseif !aok && !bok
        segclipd[(pta, ptb)] = segclipd[(ptb, pta)] = nothing
    end
end
segclipd

### cell segment by particle
for ed in voronoiedges(tess)
    push!(getgena(ed).cell, (a=geta(ed), b=getb(ed), ogen=1 )) #getgenb(ed)))
    push!(getgenb(ed).cell, (a=geta(ed), b=getb(ed), ogen=2 )) #getgena(ed)))
end

## order cell segment in a cycle
pt = pts[2]
for pt in pts
    segs = pt.cell
    ui = [2:length(segs);]
    osegs = segs[1:1]
    eend = segs[1].b
    while length(ui) > 0
        noi = findfirst( [segs[ui[i]].a == eend for i in 1:length(ui)] )
        if noi != nothing
            nseg = segs[ui[noi]]
            push!(osegs, nseg)
            eend = nseg.b
        else
            noi = findfirst( [segs[ui[i]].b == eend for i in 1:length(ui)] )
            (noi == nothing) && error("not a complete cycle")
            nseg = segs[ui[noi]]
            push!(osegs, (a=nseg.b, b=nseg.a, ogen=nseg.ogen))
            eend = nseg.a
        end
        splice!(ui, noi)
    end
    pt.cell = osegs
end

## clip cell segment
pt = pts[2]
for pt in pts
    segs = pt.cell
    csegs = Any[]
    # set lastpt if cycle starts with a point outside
    if haskey(segclipd, (segs[1].b, segs[1].a))
        lastpt = segclipd[(segs[end].a, segs[end].b)]
    end

    # clip cell walls
    for seg in segs
        if haskey(segclipd, (seg.a, seg.b))
            cpt = segclipd[(seg.a, seg.b)]
            if cpt == nothing
            else
                push!(csegs, (a=seg.a, b=cpt, ogen=seg.ogen))
                lastpt = cpt
            end
        elseif haskey(segclipd, (seg.b, seg.a))
            cpt = segclipd[(seg.b, seg.a)]
            if cpt == nothing
            else
                push!(csegs, (a=lastpt, b=cpt, ogen=nothing))
                push!(csegs, (a=cpt, b=seg.b, ogen=seg.ogen))
                lastpt = cpt
            end
        else
            push!(csegs, seg)
        end
    end
    pt.cell = csegs
end


pts

pd = [ [ [(x = getx(getfield(cpt,sym)), y = gety(getfield(cpt,sym)), id = pt.index) for sym in [:a, :b] ]
        for cpt in pt.cell]
            for pt in pts]
pd = vcat([vcat(pde...) for pde in pd]...)
pd = [ (x=p.x, y=p.y, id=p.id, o=i) for (i,p) in enumerate(pd) ]

pd2 = NamedTuple{(:x,:y,:idx), Tuple{Float64, Float64, Int64}}[]
for (idx, ed) in enumerate(voronoiedges(tess))
    push!(pd2, (x=getx(geta(ed)), y=gety(geta(ed)), idx=idx))
    push!(pd2, (x=getx(getb(ed)), y=gety(getb(ed)), idx=idx))
end

pd3 = [(x = getx(pt), y = gety(pt)) for pt in pts]

@vlplot(background=:lightgrey, width=400, height=400) +
    @vlplot(data= pd,
        x={:x, typ=:quantitative, scale={domain=[0.,3.]}},
        y={:y, typ=:quantitative, scale={domain=[0.,3.]}},
        color="id:n", order="o:n", mark=:area, detail="id:n") +
    @vlplot(data= pd3,
            x={:x, typ=:quantitative, scale={domain=[1.,2.]}},
            y={:y, typ=:quantitative, scale={domain=[1.,2.]}},
            color={value=:black}, size={value=3},
            mark=:point) +
    @vlplot(data=pd2,
        x={:x, typ=:quantitative, scale={domain=[1.,2.]}},
        y={:y, typ=:quantitative, scale={domain=[1.,2.]}},
        detail="idx:n", mark={:line, clip=true})
        


#######



########################################################################



### calculate acceleration
foreach(pt -> pt.ax = pt.ay = 0., pts)
# eds = collect(voronoiedges(tess))
# findall( [ getgena(ed).index == 2 for ed in eds] )
# findall( [ getgenb(ed).index == 2 for ed in eds] )



ied,ed = 12, eds[12]

for (ied, ed) in enumerate(voronoiedges(tess))
    ((getgena(ed).index == 0) || (getgenb(ed).index == 0)) && continue

    p1x, p1y, p2x, p2y = segcli[ied].p1x, segcli[ied].p1y, segcli[ied].p2x, segcli[ied].p2y
    Aij = segcli[ied].l

    pax, pay = getx(getgena(ed)), gety(getgena(ed))
    pbx, pby = getx(getgenb(ed)), gety(getgenb(ed))

    Rij = sqrt(abs2(pax-pbx) + abs2(pay-pby))
    eijx = (pbx - pax) / Rij
    eijy = (pby - pay) / Rij
    midAx = (p1x + p2x) / 2
    midAy = (p1y + p2y) / 2
    midijx = (pax + pbx) / 2
    midijy = (pay + pby) / 2

    cijx, cijy = midAx - midijx, midAy - midijy
    Pi, Pj = getgena(ed).pressure, getgenb(ed).pressure

    fx = -Aij * ( (Pi + Pj) * eijx / 2. + (Pj - Pi) * cijx / Rij )
    fy = -Aij * ( (Pi + Pj) * eijy / 2. + (Pj - Pi) * cijy / Rij )

    # PPO regularizing term
    Pij = 4 * getgena(ed).mass / ( Rij + getgena(ed).perimeter )
    Pji = 4 * getgenb(ed).mass / ( Rij + getgenb(ed).perimeter )
    #
    force = -Aij * κ * ( Pij - Pi + Pji - Pj) / 2.
    fx -= force * eijx
    fy -= force * eijy

    getgena(ed).ax += fx / getgena(ed).mass
    getgena(ed).ay += fy / getgena(ed).mass
    getgenb(ed).ax -= fx / getgenb(ed).mass
    getgenb(ed).ay -= fy / getgenb(ed).mass
end
maxa = maximum( sqrt(pt.ax*pt.ax + pt.ay*pt.ay) for pt in pts)
