include("setup2.jl")

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

##### test   #################################
Npart = 100
Npart = 6

pts = part_gen1(Npart)
calc_mass(pts)
pts2 = borderize(pts)
tess2 = DelaunayTessellation2D{Particle}(length(pts2))
push!(tess2, pts2)
plpoints2(tess2, pts2)

tess = DelaunayTessellation2D{Particle}(length(pts))
push!(tess, pts)
plpoints2(tess, pts, [0.5,2.5])




(mp, ma, mvit, dt2), tess = oneloop(pts, 0.01)
plpoints2(tess, pts)


pts2 = vcat(pts, ghostpts)
tess = DelaunayTessellation2D{Particle}(length(pts2))
push!(tess, pts2)
plpoints2(tess, pts2)

###############################################################################
# Run loop
###############################################################################
Npart = 250

pts = part_gen3(Npart)
for pt in pts
    if pt._y < 1.45
        pt.density = 0.5
    else
        pt.density = 1.0
    end
end
calc_mass(pts)


(mp, ma, mvit, dt2), tess = oneloop(pts, 0.01)
plpoints2(tess, pts)
plpoints3(pts)

dt = 0.01
t = 0.
for i in 1:10000
    global tess, t
    (mp, ma, mvit, dt2), tess = oneloop(pts, dt)
    t += dt2
    isfinite(ma) || break
end
t
plpoints3(pts)
plpoints2(tess, pts)

tess2 = DelaunayTessellation2D{Particle}(length(pts))
push!(tess, pts)

tess0 = DelaunayTessellation2D{Particle}(length(pts))
push!(tess, pts)

plpoints2(tess, pts)


###############################################################################
# Profiling
###############################################################################
using Profile

Npart = 100
pts = part_gen3(Npart)
calc_mass(pts)

Profile.clear()
@profile for i in 1:100; oneloop(pts,0.01); end

open("/tmp/prof.txt", "w") do s
    Profile.print(IOContext(s, :displaysize => (24, 500)), maxdepth=18)
end

Profile.print()

###############################################################################
# ref time
###############################################################################
Npart = 1000
pts = part_gen3(Npart)
calc_mass(pts)

@time for i in 1:100; oneloop(pts,0.01); end
# 3.184856 seconds (58.36 M allocations: 1.609 GiB, 12.30% gc time)
# setup2 : 2.837733 seconds (49.59 M allocations: 893.880 MiB, 4.23% gc time)
# 3.095939 seconds (48.88 M allocations: 880.696 MiB, 3.32% gc time)
# 2.291700 seconds (45.27 M allocations: 820.231 MiB, 4.21% gc time)
# 2.263285 seconds (45.00 M allocations: 805.174 MiB, 4.64% gc time)


@time for i in 1:1000; borderize(pts); end
# 8.112397 seconds (80.21 M allocations: 1.862 GiB, 4.02% gc time)
# 14.298300 seconds (198.71 M allocations: 4.010 GiB, 15.74% gc time) push!(ghost)
# 13.999349 seconds (193.52 M allocations: 3.919 GiB, 15.32% gc time) foreach(push!(pt))
# 19.000086 seconds (282.49 M allocations: 5.562 GiB, 18.88% gc time) new tess
