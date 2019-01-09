using Random
using StatsBase
using KernelDensity
using Plots

struct Point{T}
	x::T
	y::T
	z::T
end


rng = MersenneTwister(0)

punti = [Point{Float64}(rand(rng), rand(rng), rand(rng))  for i = 1:500000 ]


ix = 1/64.
iy = 1/64.
iz = 1/64.

p_id = []

for i = 1:length(punti)
	if punti[i].x < ix
		if punti[i].y < iy
			if punti[i].z < iz
				push!(p_id, punti[i])
			end
		end
	end
end

println(p_id)

#println(punti[p_id])

#println(length(p_id))




