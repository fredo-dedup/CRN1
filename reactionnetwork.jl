using Distributions
using VegaLite
using SparseArrays
using StatsBase

using Base.Iterators
import Base.Iterators: product

# lotke volterra
Ne = 2
Nr = 3
A = [1 0 ; 1 1 ; 0 1]
B = [2 0 ; 0 2 ; 0 0]
v = [1., 0.5, 1.]


it = product([0:2;], [0:2;], [0:2;], [0:2;], [0:2;], [0:2;],
             [0:2;], [0:2;], [0:2;], [0:2;], [0:2;], [0:2;])
length(it) # 531441
# réduit 60.000 cas en retirant les symétries

# initialisation
(coefs, state) = iterate(it)

ngrp = 5000
Ne = 2*ngrp
Nr = 3*ngrp
A = spzeros(Int64, Nr, Ne) ;  B = spzeros(Int64, Nr, Ne)
v = zeros(Nr)
# cnt = 0
for i in 1:ngrp
    global coefs, state
    switched = (coefs[2], coefs[1], coefs[4], coefs[3],
                coefs[6], coefs[5], coefs[8], coefs[7],
                coefs[10], coefs[9], coefs[12], coefs[11])

    if coefs >= switched

        r1 = (coefs[1], coefs[2], coefs[3], coefs[4])
        r2 = (coefs[5], coefs[6], coefs[7], coefs[8])
        r3 = (coefs[9], coefs[10], coefs[11], coefs[12])

        if (r1 >= r2) && (r2 >= r3)
            # global cnt
            # cnt += 1
            ir = 3 * i - 2 ; ie = 2*i-1
            A[ir,   ie:ie+1] = [coefs[1]  coefs[2] ]
            B[ir,   ie:ie+1] = [coefs[3]  coefs[4] ]
            A[ir+1, ie:ie+1] = [coefs[5]  coefs[6] ]
            B[ir+1, ie:ie+1] = [coefs[7]  coefs[8] ]
            A[ir+2, ie:ie+1] = [coefs[9]  coefs[10]]
            B[ir+2, ie:ie+1] = [coefs[11] coefs[12]]

            v[ir:ir+2] = [1., 0.5, 0.33]
        end
    end
    (coefs, state) = iterate(it, state)
end
# cnt
state

Nstep = 5000
cs = calc(A, B, v, dt=0.01, Nstep=Nstep)
# plotcs(cs)

## recherche des réactants qui oscillent
tol = 1e-3
ttt = diff(cs[:,1000:end], dims=2)
ttt = map(v -> (v > tol) ? 1.0 : ((v < -tol) ? -1.0 : 0.0), ttt)
ttt = mapslices(ttt, dims=2) do v
    cnt, stt = 0, v[1]
    for x in v[2:end]
        x == 0 && continue
        if x != stt
            cnt +=1
            stt = x
        end
    end
    cnt
end
specselect = findall(n -> n > 5, vec(ttt))


ttt[9413,:]
countmap(vec(ttt))

cs2 = vec([ (i=i, j=j, val=cs[i,j]) for i in specselect, j in 1:20:Nstep])
# cs2 = vec([ (i=i, j=j, val=cs[i,j]) for i in [rand(1:Ne,50);], j in 100:10:Nstep])
cs2 |> @vlplot(background=:lightgrey, height=200, config={view={strokeWidth=0},scale={rangeStep=5} }) +
 @vlplot(mark=:line, x="j:ordinal", y="val:q", color="i:n")

cs2 = vec([ (i=i, j=j, val=cs[i,j]) for i in [8925,8926], j in 1:20:Nstep])


###  5140 / 5139
idx = round(Int64, 6060/2)
Array([A[ idx*3-2:idx*3, idx*2-1:idx*2] B[ idx*3-2:idx*3, idx*2-1:idx*2]] )





# re test
Ne = 2
Nr = 3
A = sparse([2 0 ; 1 2 ; 0 1])
B = sparse([2 2 ; 2 1 ; 0 0])
v = [1., 0.5, 0.33]

cs = calc(A, B, v, dt=0.01, Nstep=Nstep)

## recherche des réactants qui oscillent
ttt = abs.(diff(sign.(diff(cs[:,1000:end], dims=2)), dims=2))
specselect = findall(n -> n > 1, vec(sum(map(x -> x==2, ttt), dims=2)))

cs2 = vec([ (i=i, j=j, val=cs[i,j]) for i in 1:Ne, j in 1:20:Nstep])
# cs2 = vec([ (i=i, j=j, val=cs[i,j]) for i in [rand(1:Ne,50);], j in 100:10:Nstep])
cs2 |> @vlplot(background=:lightgrey, height=200, config={view={strokeWidth=0},scale={rangeStep=5} }) +
 @vlplot(mark=:line, x="j:ordinal", y="val:q", color="i:n")

cs2 = vec([ (i=i, j=j, val=cs[i,j]) for i in [7940,7939], j in 1000:20:Nstep])



### large system more constrained
Ne = 100
Nr = 200


####
sum(A)
sum(B)
mapslices(sum, A, dims=1)
mapslices(sum, A, dims=1)
mapslices(sum, A, dims=2) - mapslices(sum, B, dims=2)

### every if a -> b then ∃  b -> a
Cb = Array(A' * B)
findall(!iszero, A[:,3])

Cb = min.(1, Cb * Cb)
count(!iszero, Cb)
count(v -> v < 0, Cb - Cb')
count(v -> v > 0, Cb - Cb')


Nstep = 3000
cs = calc(A, B, v, dt=0.01, Nstep=Nstep)
plotcs(cs)

## recherche des réactants qui oscillent
ttt = abs.(diff(sign.(diff(cs[:,100:end], dims=2)), dims=2))
specselect = findall(n -> n > 1, vec(sum(map(x -> x==2, ttt), dims=2)))
ttt2 = vec(sum(map(x -> x==2, ttt), dims=2))
ttt2[72]

cs2 = vec([ (i=i, j=j, val=cs[i,j]) for i in specselect, j in 100:20:Nstep])
cs2 = vec([ (i=i, j=j, val=cs[i,j]) for i in [rand(1:Ne,50);], j in 100:10:Nstep])
cs2 |> @vlplot(background=:lightgrey, height=500, config={view={strokeWidth=0},scale={rangeStep=5} }) +
 @vlplot(mark=:line, x="j:ordinal", y="val:q", color="i:n")



############ calc and plot


function calc(A2, B2, v; Nstep=1000, dt=0.1)
    cs = Array{Float64}(undef, Ne, Nstep)
    cs[:,1] = rand(Ne)
    for i in 2:Nstep
       s = exp.(A2 * log.(cs[:,i-1])) .* v
       dμ = (B2' - A2') * s
       cs[:,i] = min.(10., max.(cs[:,i-1] + dμ * dt, 1e-10))
    end
    # 100x200 : 50-66 ms  (x5)
    # 500x1000 : 260-400 ms (x20)

    cs
end

function calc2(A2, B2, v; Nstep=1000, dt=0.1)
    cs = Array{Float64}(undef, Ne, Nstep)
    cs[:,1] = rand(Ne)
    for i in 2:Nstep
       s = exp.(A2 * log.(cs[:,i-1])) .* v
       dμ = (B2' - A2') * s
       cs[:,i] = min.(1., max.(cs[:,i-1] + dμ * dt, 1e-10))
    end
    # 100x200 : 50-66 ms  (x5)
    # 500x1000 : 260-400 ms (x20)

    cs
end


function plotcs(cs)
    cs2 = vec([ (i=i, j=j, val=cs[i,j]) for i in 1:Ne, j in 1:30:Nstep])
    # cs2 = vec([ (i=i, j=j, val=cs[i,j]) for i in 1:Ne, j in 500:600])
    cs2 |> @vlplot(background=:lightgrey, config={view={strokeWidth=0},scale={rangeStep=5} }) +
       @vlplot(mark=:rect, x="j:ordinal", y="i:ordinal",
          color="val:q")

    # cs2 = vec([ (i=i, j=j, val=cs[i,j]) for i in 1:50, j in 1:10:Nstep])
    # cs2 |> @vlplot(background=:lightgrey, config={view={strokeWidth=0},scale={rangeStep=5} }) +
    #  @vlplot(mark=:line, x="j:ordinal", y="val:q", color="i:o")
end

# prepare # 1
function prepare1()
    A = spzeros(Int64, Nr, Ne); B = spzeros(Int64, Nr, Ne)
    for i in 1:Nr
        nb = rand(Poisson(2))
        r = rand(1:Ne, nb)
        foreach(ri ->  A[i, ri] = rand(Categorical([0.6, 0.2, 0.1, 0.1])), r)
        p = rand(1:Ne, nb)
        foreach(ri ->  B[i, ri] = rand(Categorical([0.6, 0.2, 0.1, 0.1])), p)
    end
    v = rand(Nr)
    A, B, v
end



# prepare # 2
function prepare2()
    A = spzeros(Int64, Nr, Ne)
    B = spzeros(Int64, Nr, Ne)
    for i in 1:Nr
        nb = 1 + rand(Poisson(2))  # number of species (both sides)
        center = rand(1:Ne) + Ne
        species = ((center .+ round.(Int64, rand(Normal(0,5), nb))) .% Ne) .+ 1
        println(species)
        splitidx = nb == 1 ? 0 : rand(1:nb-1)
        rand(1:1)
        for (j,s) in enumerate(species)
            if j > splitidx
                B[i, s] = 1
            else
                A[i, s] = 1
            end
        end
        print(join(species[1:splitidx], "+"))
        print(" -> ")
        println(join(species[splitidx+1:end], "+"))

    end
    v = rand(Nr)
    A, B, v
end


# prepare # 3
function prepare3(;poisson=3, normal=2)
    A = spzeros(Int64, Nr, Ne); B = spzeros(Int64, Nr, Ne)
    for i in 1:Nr
        nb = (1 + rand(Poisson(poisson))) * 2  # number of species (both sides)
        center = rand(1:Ne) + Ne
        species = ((center .+ round.(Int64, rand(Normal(0,normal), nb))) .% Ne) .+ 1
        for (j,s) in enumerate(species)
            if j > nb / 2
                B[i, s] += 1
            else
                A[i, s] += 1
            end
        end
        print(join(species[1:splitidx], "+"))
        print(" -> ")
        println(join(species[splitidx+1:end], "+"))
    end
    v = rand(Nr)
    v = ones(Nr)
    A, B, v
end
