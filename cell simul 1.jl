################################################################################
# premier prototype organisme simulé
################################################################################

using ProgressMeter

include("setup cell simul 1.jl")

###############################################################################
# Physical constants
###############################################################################
# isentropic gas, H₂O γ ~=  1.33
γ = 1.3
# cell rounding force
κ = 0.01
# gravity
g = Point2D(0., -1.0)
# indicative time step
dt = 0.01

###############################################################################
# Run loop
###############################################################################
# TODO : adapter pression au départ pour éviter l'onde de pression ensuite
# TODO : meilleure echelle de couleur

Npart = 1000

pts = part_gen4(Npart)
for pt in pts
    if abs(pt - Point2D(1.5, 1.5)) < 0.02
        pt.density = 5.
        pt.orgindex = 1
    else
        pt.density = 7.
        pt.orgindex = 0
    end
end
calc_mass(pts, 10.)

tess = borderize(pts)
qp = quickplot3(pts, tess)

sum(pt.mass for pt in pts)
sum(pt.surface for pt in pts)
sum(pt.pressure for pt in pts) / length(pts)
extrema(pt.density for pt in pts)

t = 0.
δ = 0.02
@showprogress for j in 1:1000
    global tess, t
    (ma,_,_,dt2), tess = oneloop(pts, dt)
    t += dt2
    if div(t-dt2, δ) != div(t, δ)  # image point passed ?
        qp = quickplot3(pts, tess);
        # display(quickplot2(pts, tess))
        draw(PNG("/tmp/csim1/cs$j.png", 10cm, 10cm), qp)
        # println("/tmp/rt$j.png, t = $t")
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
