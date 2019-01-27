using Random
using Statistics
using StatsBase
using Plots
import Dates
import Base.+
import Base.-


include("Point_Struct.jl")
include("fun.jl")
include("constants.jl")
include("integrate.jl")









#################### Integration #####################

calls = Int64(1E3)
it    = 5
dimensions = 6.

Φ₀ = 0.

Φₛ = 0.0         # start point
vₛ = 2*π / 60.   # speed rad/sec
dt = 1           # Time interval

Start_Time = 0.0
Stop_Time  = 10.0

time = Array{Float64, 1}(undef, 0)
Coo  = Array{Float64, 1}(undef, 0)

Step_start = Start_Time / dt
Step_stop  = Stop_Time / dt


xu = Point{Float64}(-L₀,     (2. / .5)*π + θᵦ, Φ₀ + θᵦ, undef, undef, undef)
xl = Point{Float64}(L₀, (2. / .5)*π - θᵦ, Φ₀ - θᵦ , undef, undef, undef )

start = Dates.value(Dates.now())

for t = Step_start:Step_stop
    xu_p = Point{Float64}(-L₀, (2. / .5)*π + θᵦ, (Φ₀ + (t*dt) * vₛ) + θᵦ, undef, undef, undef)
    xl_p = Point{Float64}(L₀, (2. / .5)*π - θᵦ, (Φ₀ + (t*dt) * vₛ) - θᵦ, undef, undef ,undef)
    result = Integrate(corr, xl, xu, xl_p, xu_p, calls, it, dimensions)
    println(t*dt," ", result)
    append!(time, t*dt)
    append!(Coo, result)
end

stop = Dates.value(Dates.now())

println("Performance: ", ( stop - start ) / (1000*(Step_stop - Step_start)) , " sec/int" )

plot(time, Coo, seriestype=:scatter)
