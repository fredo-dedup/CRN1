include("setup4.jl")


################################################################################

using Compose, Colors

cmap = colormap("RdBu")

Comp2 = Tuple{Float64, Float64}

function quickplot(pts, tess)
    gonsbycol = [ Vector{Comp2}[] for i in 1:100 ]

    for ed in voronoiedges(tess)
        pa, pb = getgena(ed), getgenb(ed)
        p1x, p1y = getx(geta(ed)), gety(geta(ed))
        p2x, p2y = getx(getb(ed)), gety(getb(ed))
        if pa.index > 0
            pax, pay = getx(pa), gety(pa)
            t1, t2, t3 = (p1x, p1y), (p2x, p2y), (pax, pay)
            col = max(1, 100 - round(Int64, pa.density))
            push!(gonsbycol[col], [t1, t2, t3])
        end
        if pb.index > 0
            pbx, pby = getx(pb), gety(pb)
            t1, t2, t3 = (p1x, p1y), (p2x, p2y), (pbx, pby)
            col = max(1, 100 - round(Int64, pb.density))
            push!(gonsbycol[col], [t1, t2, t3])
        end
    end

    # col = cmap[max(1, 100 - round(Int64, pb.density))]

    set_default_graphic_size(10cm,10cm)

    gonset = []
    for c in 1:100
        (length(gonsbycol[c]) == 0) && continue
        push!( gonset, (context(), polygon(gonsbycol[c]), fill(cmap[c])) )
    end

    compose(context(units=UnitBox(1.33, 1.67, 0.34, -0.34)),
            gonset..., fillopacity(1.0))
end


###############################################################################
# Run loop
###############################################################################
# TODO : adapter pression au départ poiur éviter l'onde de presion ensuite
# TODO : meilleure echelle de couleur

Npart = 10000
dt = 0.01
κ, g = 0.01, Point2D(0., -1.0)

pts = part_gen4(Npart)
for pt in pts
    pt.density = abs(pt - Point2D(1.5, 1.0)) < 0.4 ? 5. : 10.
end
calc_mass(pts, 10.)

tess = borderize(pts)
plpoints2(tess, pts)

info, tess = oneloop(pts, 0.01)
plotval(pts, field=:density)
sum(pt.mass for pt in pts)
sum(pt.surface for pt in pts)
sum(pt.pressure for pt in pts) / length(pts)

t = 0.
δ = 0.02
for j in 4001:10000
    global tess, t
    (ma,_,_,dt2), tess = oneloop(pts, dt)
    t += dt2
    if div(t-dt2, δ) != div(t, δ)  # image point passed ?
        qp = quickplot(pts, tess);
        # display(qp)
        draw(PNG("/tmp/rt2/rt$j.png", 10cm, 10cm), qp)
        println("/tmp/rt$j.png, t = $t")
    end
    isfinite(ma) || break
end
t

plotval(pts, field=:density)
plotval(pts, field=:mass)
plotval(pts, field=:pressure)
plotval(pts, field=:entropy)
plotval(pts, field=:surface)
plotval(pts, field=:density)

plotspeed(pts)
plpoints2(tess, pts)
plpoints2v(tess, pts)

###############################################################################
# Energy conservation
###############################################################################

energy(pt::Particle) =
    pt.mass * (
        pt.v ⋄ pt.v +
        pt.entropy * (pt.mass / pt.surface) ^ (γ - 1.) / (γ - 1.) -
        g ⋄ pt
    )

ekin(pt::Particle)       =  pt.mass * (pt.v ⋄ pt.v)
ethermal(pt::Particle)   =  pt.mass * pt.entropy *
    (pt.mass / pt.surface) ^ (γ - 1.) / (γ - 1.)
epotential(pt::Particle) =  pt.mass * (g ⋄ pt)

Npart = 1000
dt = 0.01

pts = part_gen4(Npart)
for pt in pts
    pt.mass = 1.0
end
equalize_pressure(pts, 10.)

for pt in pts
    pt.density = 1.0
end
calc_mass(pts)



κ, g = 0.01, Point2D(0., -1.0)
# g = Point2D(0., 0.)

meas = typeof((t=0.2, name=:press, val=12.45))[]

t = 0.
for i in 1:1000
    global tess, t
    (ma,_,_,dt2), tess = oneloop(pts, dt)
    isfinite(ma) || break
    t += dt2
    push!(meas, (t=t, name=:avgpre, val = sum(pt.pressure for pt in pts) / length(pts)))
    push!(meas, (t=t, name=:mass, val = sum(pt.mass for pt in pts)))
    push!(meas, (t=t, name=:surf, val = sum(pt.surface for pt in pts)))
    push!(meas, (t=t, name=:energy, val = sum(energy,pts)))
    push!(meas, (t=t, name=:ekin, val = sum(ekin,pts)))
    push!(meas, (t=t, name=:ethermal, val = sum(ethermal,pts)))
    push!(meas, (t=t, name=:epotential, val = sum(epotential,pts)))
end

meas
@vlplot(
    data=filter(el -> el.name in [:ekin, :energy, :ethermal, :epotential], meas),
    background=:lightgrey,width=400, height=200,
    mark=:line, color="name:n",
    x={:t, typ=:quantitative, scale={zero=false}},
    y="val:q", resolve={axis={y="independent"}})

@vlplot(data=meas,
    background=:lightgrey,width=400, height=200,
    mark=:line, color="name:n", row="name:n",
    x={:t, typ=:quantitative, scale={zero=false}},
    y="val:q", resolve={axis={y="independent"}})
