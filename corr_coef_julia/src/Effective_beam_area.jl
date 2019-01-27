using Random
using Statistics
using StatsBase
using Plots
import Base.+
import Base.-


include("Point_Struct.jl")
include("fun.jl")
include("constants.jl")

# This script reproduce the turbulance distribution function (a gaussian)
# Figure 4 in Errard's article

p_0 = [Point{Float64}(20000, 0.9, i, undef, undef, undef)() for i = (1:25000)*( (2*π) / 25000  )  ] # must be constant
p_1 = [Point{Float64}(20000, i, 0.9, undef, undef, undef)() for i = (1:25000)*( (2*π) / 25000 )  ] # must be constant
p_2 = [Point{Float64}(i, 0.9, 0.9, undef, undef, undef)() for i = (10:1000000)/1000000 * 100000  ] # must be constant


w   = [ω(o) for o in p_2]
BB = [B(i) for i in p_2]

z   = [i.z for i in p_2]

plot(  BB[2:length(BB)] / (2*π*(ω₀^2)) , z[2:length(z)], yaxis=:log)#, xlims=(0,1050))
#plot( w[2:length(z)] .* ( (ω₀ * π) / ( λ * 100000 )), z[2:length(z)])#, yaxis=:log)#, xlims=(0,1050))
