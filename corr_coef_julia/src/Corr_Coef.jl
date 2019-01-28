using Random
using Statistics
using StatsBase
#using Plots
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

Start_Time = parse(Float64, ARGS[1])
Stop_Time  = parse(Float64, ARGS[2])
file_name  = ARGS[3]

Step_start = Start_Time / dt
Step_stop  = Stop_Time / dt


time_s = Array{Float64, 1}(undef, 0)
Coo  = Array{Float64, 1}(undef, 0)




xu = Point{Float64}(-L₀,     (2. / .5)*π + θᵦ, Φ₀ + θᵦ, undef, undef, undef)
xl = Point{Float64}(L₀, (2. / .5)*π - θᵦ, Φ₀ - θᵦ , undef, undef, undef )

start = Dates.value(Dates.now())




for t = Step_start:Step_stop
    xu_p = Point{Float64}(-L₀, (2. / .5)*π + θᵦ, (Φ₀ + (t*dt) * vₛ) + θᵦ, undef, undef, undef)
    xl_p = Point{Float64}(L₀, (2. / .5)*π - θᵦ, (Φ₀ + (t*dt) * vₛ) - θᵦ, undef, undef ,undef)
    result = Integrate(corr, xl, xu, xl_p, xu_p, calls, it, dimensions, 0.0)
    println(t*dt," ", result)
    append!(time_s, t*dt)
    append!(Coo, result)
end

ord = sortperm(time_s)
time_s = time_s[ord]
Coo  = Coo[ord]

# Print to file
f = open(file_name, "w")
for j = 1:min(length(Coo), length(time_s))
    println(f, time_s[j], " ", Coo[j])
end
close(f)


#stop = Dates.value(Dates.now())
#println("Performance: ", ( stop - start ) / (1000 * (Step_stop - Step_start)), " sec/step" )
