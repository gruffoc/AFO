using Random
using Statistics
using StatsBase
using Plots
import Base.+
import Base.-


include("Point_Struct.jl")
include("fun.jl")
include("constants.jl")
include("integrate.jl")









#################### Integration #####################

calls = Int64(5E4)
dimensions = 6.

Φ₀ = 0.

Φₛ = 0.0         # start point
vₛ = 2*π / 60.   # speed rad/sec
dt = 0.1         # Time interval
global time = 0.0       # start time

Start_Time = 241
Stop_Time  = 300 # 10 minuti di osservazione

Step_start = Start_Time / dt
Step_stop  = Stop_Time / dt


xu = Point{Float64}(-L₀,     (2. / .5)*π + θᵦ, Φ₀ + θᵦ, undef, undef, undef)
xl = Point{Float64}(L₀, (2. / .5)*π - θᵦ, Φ₀ - θᵦ , undef, undef, undef )

for t = Step_start:Step_stop
    xu_p = Point{Float64}(-L₀, (2. / .5)*π + θᵦ, (Φ₀ + (t*dt) * vₛ) + θᵦ, undef, undef, undef)
    xl_p = Point{Float64}(L₀, (2. / .5)*π - θᵦ, (Φ₀ + (t*dt) * vₛ) - θᵦ, undef, undef ,undef)
    result = Integrate(corr, xl, xu, xl_p, xu_p, calls, dimensions)
    println(t*dt," ", result)
end
