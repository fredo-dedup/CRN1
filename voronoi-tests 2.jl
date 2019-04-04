include("setup.jl")

##### test   #################################
width = max_coord - min_coord
Npart = 5
pts = [Particle(min_coord + rand() * width, min_coord + rand() * width,
           i, 1. / Npart, 1., NaN, NaN, NaN,
           0., 0., 0., 0., []) for i in 1:Npart]

coordrg = min_coord+0.05:sqrt(1/Npart):max_coord-0.05
pts = vec( [ Particle(x, y + rand()*0.001,
                i, 1. /Npart, 1., NaN, NaN, NaN,
                0., 0., 0., 0., []) for
                (i, (x,y)) in enumerate(product(coordrg, coordrg)) ] )

# balance mass
tess = DelaunayTessellation2D{Particle}(100)
push!(tess, pts)

## prepare the cell around each particle
clipcells(tess, pts)

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
    # last segment closing the loop
    curdx, curdy = getx(pt.walls[1].a) - getx(pos), gety(pt.walls[1].a) - gety(pos)
    pt.surface   += abs(prevdx*curdy - curdx*prevdy) / 2

    pt.mass = pt.surface
end

res, tess = oneloop(pts, 0.01)
plpoints2(tess, pts)


dt = 0.01
t = 0.
for i in 1:1000
    global tess, t
    (mp, ma, mvit, dt2), tess = oneloop(pts, dt)
    t += dt2
    isfinite(ma) || break
end
t
plpoints2(tess, pts)

eds = collect(voronoiedges(tess))
ied, ed = 12, eds[12]
