using Random
using Statistics
using StatsBase
using Plots
import Base.+
import Base.-


const λ  = 0.002 #in mm => 150GHz
const θᵦ = (π / 180.) * 0.17   #0.17 deg for STRIP W-band
const ω₀ = λ / ( π * θᵦ ) # θᵦ is the beam opening angle in radians
const z_atm = 3E4 #depends on the observation site
const Tᵧ = 290 # ground temperature in kelvin
const χ₁_0 = 1 # in m^{-2}
const χ₂_0 = 1 # in m^{-2}
const z_0  = 2E3 # the quote where water vapor vanish
const L₀   = 100 # in m, the turbulance correlation length



struct Point{T}
	x::T
 	y::T
	z::T
end


function(p::Point)()
	x = p.x * cos(p.y) * cos(p.z)
	y = p.x * cos(p.y) * sin(p.z)
	z = p.x * sin(p.y)
	return Point(x, y, z)
end


Base.:+(a::Point{Float64}, b::Point{Float64}) = Point{Float64}(a.x+b.x , a.y+b.y, a.z+b.z)
Base.:-(a::Point{Float64}, b::Point{Float64}) = Point{Float64}(a.x-b.x , a.y-b.y, a.z-b.z)
Base.:*(a::Point{Float64}, b::Point{Float64}) = Point{Float64}(a.x*b.x , a.y*b.y, a.z*b.z)
Base.:/(a::Point{Float64}, b::Float64) = Point{Float64}(a.x/b , a.y/b, a.z/b)


r = 0:20000

punti = [Point{Float64}(i, π / 2 , 0.)()   for i in r]
z = [i.z for i in punti]
ww = [ ω₀ * ( 1 + ( (λ * z) / (π * ω₀^2)  )^2)^0.5  for z in r]

plot( ww[2:length(r)] ./ ww[2]  , z[2:length(r)], yaxis=:log, xlims=(0,2), ylims = (1,2E4))
