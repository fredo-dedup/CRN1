include("setup.jl")

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



function borderize(pts)
    tess = DelaunayTessellation2D{Particle}(length(pts))
    push!(tess, pts)

    filter!(pts) do pt

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
end


function calc_mass(pts)
    tess = DelaunayTessellation2D{Particle}(length(pts))
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
end

##### test   #################################
Npart = 100

pts = part_gen3(Npart)
calc_mass(pts)

pts = part_gen1(5)

extrema( getx(pt) for pt in pts )
extrema( gety(pt) for pt in pts )

tess = DelaunayTessellation2D{Particle}(length(pts))
push!(tess, pts)
plpoints2(tess, pts)


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
    elseif getgenb(ve).index != 0
        up    && (ptsflip[getgenb(ve).index, 1] = true)
        down  && (ptsflip[getgenb(ve).index, 2] = true)
        right && (ptsflip[getgenb(ve).index, 3] = true)
        left  && (ptsflip[getgenb(ve).index, 4] = true)
    end
end

ghostpts = Particle[]
gidx = -1
for pt in pts
    global gidx
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
length(ghostpts)
extrema( getx(pt) for pt in pts )
extrema( gety(pt) for pt in pts )

pts2 = vcat(pts, ghostpts)
tess = DelaunayTessellation2D{Particle}(length(pts2))
push!(tess, pts2)
plpoints2(tess, pts2)



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

###############################################################################
# Profiling
###############################################################################
using Profile

@profile for i in 1:100; oneloop(pts,0.01); end

Profile.print()

###############################################################################
# ref time
###############################################################################
Npart = 1000
pts = part_gen3(Npart)
calc_mass(pts)

@time for i in 1:100; oneloop(pts,0.01); end
# 3.184856 seconds (58.36 M allocations: 1.609 GiB, 12.30% gc time)
