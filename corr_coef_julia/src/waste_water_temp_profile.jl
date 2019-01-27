using Random
using Statistics
using StatsBase
using Plots
import Base.+
import Base.-


include("Point_Struct.jl")
include("fun.jl")
include("constants.jl")


# This script reproduce the Figure 3 in the Errard's paper

r = 0:20000

punti = [Point{Float64}(i, π / 2 , 0., undef, undef, undef)()   for i in r] #remember the functor!() else it is a cartesian poin
z = [i.z for i in punti]
w = [ω(o) for o in punti]
T = [Tₚ(o) for o in punti]
X = [χ₂(o,o) for o in punti]

plot( w[2:length(r)] ./ w[2]  , z[2:length(r)], yaxis=:log, xlims=(0,2), ylims = (1,2E4))
plot!( T[2:length(r)] ./ T[2], z[2:length(r)], yaxis=:log, xlims=(0,2), ylims = (1,2E4))
plot!( X[2:length(r)] ./ X[2], z[2:length(r)], yaxis=:log, xlims=(0,2), ylims = (1,2E4))
