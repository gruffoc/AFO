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



Φₛ = π/2         # start point
vₛ = 2*π / 60.   # speed rad/sec
dt = 0.3         # Time interval

Start_Time = 0   #parse(Float64, ARGS[1])
Stop_Time  = 30 #parse(Float64, ARGS[2])
file_name  = "ciao.txt" #ARGS[3]

Step_start = Start_Time / dt
Step_stop  = Stop_Time / dt


time_s = Array{Float64, 1}(undef, 0)
Coo  = Array{Float64, 1}(undef, 0)





xu = Point{Float64}(4000.0,  π/4 + θᵦ, Φₛ + θᵦ, undef, undef, undef)
xl = Point{Float64}(0.,      π/4 - θᵦ,  Φₛ -  θᵦ, undef, undef, undef )

start = Dates.value(Dates.now())




for t = Step_start:Step_stop
    tempo = Float64(t*dt)
    # Rivedere gli estremi di integrazione.
    # xu_p = Point{Float64}(5035.,   2.8902652413026098,  (Φₛ + vₛ*t*dt) + θᵦ, undef, undef, undef)
    # xl_p = Point{Float64}(0.,     0.12566370614359174,  (Φₛ + vₛ*t*dt) - θᵦ , undef, undef ,undef)
    xu_p = Point{Float64}(4000.,   π/4 + θᵦ,  (Φₛ + 0.2*sin(2*π*0.05*t*dt)) + θᵦ, undef, undef, undef)
    xl_p = Point{Float64}(0.,      π/4 - θᵦ,  (Φₛ + 0.2*sin(2*π*0.05*t*dt)) - θᵦ , undef, undef ,undef)
    result = Integrate(corr, xl, xu, xl_p, xu_p, calls, it, dimensions, tempo )
    println(t*dt," ", result)
    append!(time_s, t*dt)
    append!(Coo, result)
end

ord = sortperm(time_s)
time_s = time_s[ord]
Coo  = Coo[ord]

#
# # Print to file
# f = open(file_name, "w")
# for j = 1:min(length(Coo), length(time_s))
#     println(f, time_s[j], " ", Coo[j])
# end
# close(f)

plot(time_s, Coo/findmax(Coo)[1], seriestype=:scatter)
#
# using DSP
# period = periodogram(Coo, fs=1)
# Freq   = freq(period)
# A      = power(period)





#stop = Dates.value(Dates.now())
#println("Performance: ", ( stop - start ) / (1000 * (Step_stop - Step_start)), " sec/step" )
