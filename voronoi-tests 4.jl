include("setup3.jl")

###############################################################################
# Run loop
###############################################################################
Npart = 1000

pts = part_gen4(Npart)
for pt in pts
    if pt._y < 1.45
        pt.density = 0.5
    else
        pt.density = 1.
    end
end
calc_mass(pts)

tess = borderize(pts)
plpoints2(tess, pts)

info, tess = oneloop(pts, 0.01)
plpoints3p(pts)
plpoints3(pts)
sum(pt.mass for pt in pts)
sum(pt.surface for pt in pts)
sum(pt.pressure for pt in pts) / length(pts)

dt = 0.01
t = 0.
for i in 1:1000
    global tess, t
    (ma,_,_,_,dt2), tess = oneloop(pts, dt)
    t += dt2
    isfinite(ma) || break
end
t
plpoints3(pts)
plpoints2(tess, pts)
plpoints2v(tess, pts)

plpoints3p(pts)


tess2 = DelaunayTessellation2D{Particle}(length(pts))
push!(tess, pts)

tess0 = DelaunayTessellation2D{Particle}(length(pts))
push!(tess, pts)

plpoints2(tess, pts)


###############################################################################
# Profiling
###############################################################################
using Profile

Npart = 1000
pts = part_gen3(Npart)
calc_mass(pts)

Profile.clear()
@profile for i in 1:100; oneloop(pts,0.01); end

@profile oneloop(pts,0.01)

open("c:/temp/prof.txt", "w") do s
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

# portable
# 4.468091 seconds (45.91 M allocations: 820.515 MiB, 3.33% gc time)



@time for i in 1:1000; borderize(pts); end
# 8.112397 seconds (80.21 M allocations: 1.862 GiB, 4.02% gc time)
# 14.298300 seconds (198.71 M allocations: 4.010 GiB, 15.74% gc time) push!(ghost)
# 13.999349 seconds (193.52 M allocations: 3.919 GiB, 15.32% gc time) foreach(push!(pt))
# 19.000086 seconds (282.49 M allocations: 5.562 GiB, 18.88% gc time) new tess
