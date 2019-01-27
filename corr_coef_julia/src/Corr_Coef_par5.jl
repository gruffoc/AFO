using Random
using Statistics
using StatsBase
using Plots
import Dates
import Base.+
import Base.-
using Base.Threads

include("Point_Struct.jl")
include("fun.jl")
include("constants.jl")
include("integrate.jl")


NoT = 1


#################### Integration #####################

calls = Int64(1E3)
it    = 5
dimensions = 6.

Φ₀ = 0.

Φₛ = 0.0         # start point
vₛ = 2*π / 60.   # speed rad/sec
dt = 1           # Time interval

Start_Time = 100800.0*4
Stop_Time  = 100800.0*5

Step_start = Start_Time / dt
Step_stop  = Stop_Time / dt

Delta_Step = (Step_stop - Step_start)/NoT

time = Array{Float64, 1}(undef, 0)
Coo  = Array{Float64, 1}(undef, 0)




xu = Point{Float64}(-L₀,     (2. / .5)*π + θᵦ, Φ₀ + θᵦ, undef, undef, undef)
xl = Point{Float64}(L₀, (2. / .5)*π - θᵦ, Φ₀ - θᵦ , undef, undef, undef )

start = Dates.value(Dates.now())




@threads for i = 1:NoT
    for t = (Step_start+(i-1)*Delta_Step):(Step_start+i*Delta_Step)
        xu_p = Point{Float64}(-L₀, (2. / .5)*π + θᵦ, (Φ₀ + (t*dt) * vₛ) + θᵦ, undef, undef, undef)
        xl_p = Point{Float64}(L₀, (2. / .5)*π - θᵦ, (Φ₀ + (t*dt) * vₛ) - θᵦ, undef, undef ,undef)
        result = Integrate(corr, xl, xu, xl_p, xu_p, calls, it, dimensions)
        #println(t*dt," ", result)
        append!(time, t*dt)
        append!(Coo, result)
    end
end

ord = sortperm(time)
time = time[ord]
Coo  = Coo[ord]

for j = 1:min(length(Coo), length(time))
    println(time[j], " ", Coo[j])
end


stop = Dates.value(Dates.now())
println("Performance: ", ( stop - start ) )
