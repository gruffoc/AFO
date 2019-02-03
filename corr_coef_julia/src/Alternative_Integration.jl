using Random
using Statistics
using StatsBase
using Plots
import Base.+
import Base.-


include("Point_Struct.jl")
include("fun.jl")
include("constants.jl")

Φ₀ = 0.

Φₛ = 0.0         # start point
vₛ = 2*π / 60.   # STRIP speed rad/sec
dt = 0.1         # Time interval

Start_Time = 0   #parse(Float64, ARGS[1])
Stop_Time  = 50 #parse(Float64, ARGS[2])
file_name  = "ciao.txt" #ARGS[3]

Step_start = Start_Time / dt
Step_stop  = Stop_Time / dt


time_s = Array{Float64, 1}(undef, 0)
Coo  = Array{Float64, 1}(undef, 0)


r = 1:40000

for i in r
    θₐ = L_₀ / Float64(i)
    if θₐ > θᵦ
