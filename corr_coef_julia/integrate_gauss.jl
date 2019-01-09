using Random
using StatsBase
using Plots
import Base.+
import Base.-



struct Point{T}
	x::T
	y::T
	z::T
end

Base.:+(a::Point{Float64}, b::Point{Float64}) = Point{Float64}(a.x+b.x , a.y+b.y, a.z+b.z)
Base.:-(a::Point{Float64}, b::Point{Float64}) = Point{Float64}(a.x-b.x , a.y-b.y, a.z-b.z)
Base.:*(a::Point{Float64}, b::Point{Float64}) = Point{Float64}(a.x*b.x , a.y*b.y, a.z*b.z)



function gauss_3D(p::Point{Float64})
	return exp.(-(p.x) .^ 2 .- (p.y) .^ 2 .- (p.z) .^ 2)
end


function Integrate(fun::Function, xl::Point{Float64}, xu::Point{Float64}, cal::Int64)
	
	rng = MersenneTwister(0)
	pt = [Point{Float64}(rand(rng), rand(rng), rand(rng)) for n = 1:Int64(floor(calls/3))]
	
	
	#; Cambio coordinate
	p = ( pt .* Ref(xu-xl) ) .+ Ref(xl)

	#; Valuto l'integrale
	w = [fun(point) * (xu.x-xl.x) * (xu.y-xl.y) * (xu.z-xl.z)  for point in p ]
	w_abs = abs.(w)
	
	Integral      = sum(w)/(cal/3.)
	Integral_abs  = sum(w_abs) / (cal/3.) 

	#; Valore dell'integrale di prova
	Int_Vero = 5.56833
	println("Prima Stima **Uniforme** = ",Integral-Int_Vero)
	
	#; Calcolo la nuova densita` di probabilita`
	g = (w_abs ./ Integral_abs) / (cal/3)
	
	println("Controllo integrale pdf = ",sum(g))

	NSIDE = 3

	#; posso usare i pt non i p
	#; devo trovare un modo sensato di mettere le cose
	density = 0
	idx_sel = []
	for ix = (1:NSIDE)
		for iy = (1:NSIDE)
			for iz = (1:NSIDE)

				for idx_g = 1:length(g)
					#; in teoria l'ordinamento puo` avvenire anche facendo
					#; NSIDE * pt / i[x,y,z] = N[x,y,z] , dove N[x,y,z] = indice del BOX in
					#; x, y e z. Questo ridurrebbe al minimo le operazioni di check da fare
					#; per posizionare un punto nella giusta BOX.
					if pt[idx_g].x < (ix * (1/NSIDE)) && pt[idx_g].x >= ((ix-1)*(1/NSIDE))
						if pt[idx_g].y < (iy * (1/NSIDE)) && pt[idx_g].y >= ((iy-1)*(1/NSIDE))
							if pt[idx_g].z < (iz * (1/NSIDE)) && pt[idx_g].z >= ((iz-1)*(1/NSIDE))
								push!(idx_sel,Int64(idx_g))
							end
						end
					end
				end
				g_box      = g[idx_sel]
				density    = sum(g_box)
				hit_in_box = density * (cal/3.)
				println(density," ", Int64(floor(hit_in_box)))
				
				#p_2  = ( pt_2 .* Ref(xu-xl) ) .+ Ref(xl)
				#w_2 = [fun(point) * (xu.x-xl.x) * (xu.y-xl.y) * (xu.z-xl.z)  for point in p_2 ]
				idx_sel = []
			end
		end
	end
	
	#append!(w,w_2)
	
	Integral = (sum(w) )/(cal/3.)
	println(Integral-Int_Vero)
		

	
	return Integral
end


#################### Integration #####################

calls = 1000000

xu = Point{Float64}(10., 10., 10.)
xl = Point{Float64}(-10.,-10.,-10.)



result = Integrate(gauss_3D, xl, xu, calls)

println(result)




