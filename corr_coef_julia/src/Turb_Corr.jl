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

p_0 = Point{Float64}(4000, 0.8, 0.8, undef, undef, undef)()
p_1 = [Point{Float64}(4000, 0.8, i, undef, undef, undef)() for i = (1:1000)*( (π/2) / 1000 )  ]

dr_vec = [ diff_cart(i,p_0) for i in p_1  ]

DR = [( i.x^2 + i.y^2 + i.z^2 )^0.5  for i in dr_vec ]

w = [ χ₁(p_0, i) for i in p_1 ]

plot(DR, w, seriestype=:scatter, xlims = (0,1000))
