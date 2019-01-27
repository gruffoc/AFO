include("/home/algebrato/.julia/dev/VegasInt/src/VegasInt.jl")
using .VegasInt


function gaus(p::Point{Float64})
    return @. exp( - (p.x)^2 - (p.y)^2 - (p.z)^2 )
end


chiamate   = 1000000
dimensioni = 3

xu = Point{Float64}(1.,1.,1.)
xl = Point{Float64}(-1.,-1.,-1.)

#result = Integrate(gaus, xl, xu, chiamate, dimensioni)

#println("Final Integral = ", result)
