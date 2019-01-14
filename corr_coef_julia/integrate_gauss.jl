using Random
using Statistics
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
	
	#; Primo step, cal/3. valutazioni della funzione
	rng = MersenneTwister(0)
	pt = [Point{Float64}(rand(rng), rand(rng), rand(rng)) for n = 1:Int64(floor(calls/3))]

	#; Cambio coordinate, per mettermi nel range di integrazione scelto
	p = ( pt .* Ref(xu-xl) ) .+ Ref(xl)

	#; Valuto l'integrale
	f = [fun(point) * (xu.x-xl.x) * (xu.y-xl.y) * (xu.z-xl.z)  for point in p ]
	f_abs = abs.(f)
	
	Integral      = sum(f .* Int64(floor(calls/3)) )/(cal/3.)
	Integral_abs  = sum(f_abs .* Int64(floor(calls/3))) / (cal/3.) 

	#; Valore dell'integrale di prova
	Int_Vero = 3.33231 * cal / 3.
	println("Prima Stima **Uniforme** = ",(Integral-Int_Vero)/cal * 3.)
	
	#; Calcolo la nuova densita` di probabilita`
	g = (f_abs ./ Integral_abs)# / length(f_abs)
	
	#println("Controllo integrale pdf = ",sum(g))

	NSIDE =40 

	#; posso fare qui la separazione dei punti nelle varie BOX
	
	w_m	 = zeros(Int64(NSIDE^3))
	#coords_x = zeros(Int64(NSIDE^3))
	Int_box  = zeros(Int64(NSIDE^3))
	Int_abs_box = zeros(Int64(NSIDE^3))
	control = cal / 3.	
	dix = diy = diz = 1.0/NSIDE
	for iter = 1:20
		#; Organizzo la pdf nel BOX
		for idx_g = 1:length(g)
			box_id_x = Int64(ceil(pt[idx_g].x / dix)) -1
			box_id_y = Int64(ceil(pt[idx_g].y / diy)) -1
			box_id_z = Int64(ceil(pt[idx_g].z / diz)) -1
			ID = (box_id_x) + (box_id_y)*NSIDE + (box_id_z)*NSIDE*NSIDE
			#coords_x[ID+1] = box_id_x  * dix
			w_m[ID+1] += g[idx_g]  
		end
	
		#println("Densita` in BOX")
		#for i in 1:length(w_m)
			#println(w_m[i])
		#end
		#println("Controllo Integrale PDF = ", sum(w_m))
		
		#; Ma come cazzo si inizializzano????
		#pt_2_all = Array{Point{Float64}}
		#f_2_all  = Array{Float64,1}
		
		for ix = (1:NSIDE)
			for iy = (1:NSIDE)
				for iz = (1:NSIDE)
					p_start = Point{Float64}((ix-1) * dix , (iy-1) * diy, (iz-1) * diz)
					p_final = Point{Float64}(ix * dix, iy * diy, iz * diz)

					ID_BOX = ( (ix-1) + ( (iy-1) * NSIDE ) + ( (iz-1) * NSIDE * NSIDE ) ) + 1
					
					#; Densita` Uniforme
					density    = 1/(NSIDE*NSIDE*NSIDE)

					#; Densita` Pesata
					#density    = w_m[ID_BOX]


					hit_in_box = Int64(floor(density * (cal/3.)))
					control    += hit_in_box
					#println("Hit = ", hit_in_box)
					if hit_in_box != 0
						#println("Punti X da ", p_start.x," fino a ",p_final.x)
						#println("Punti Y da ", p_start.y," fino a ",p_final.y)
						#println("Punti Z da ", p_start.z," fino a ",p_final.z)

						pt_2 = [Point{Float64}((p_final.x-p_start.x)*rand(rng)+p_start.x, (p_final.y-p_start.y)*rand(rng)+p_start.y,(p_final.z-p_start.z)*rand(rng)+p_start.z) for n = 1:hit_in_box]
						#println(pt_2)
						#; Cambio coordinate
						p_2 = ( pt_2 .* Ref(xu-xl) ) .+ Ref(xl)
						
						f_2= [fun(point) * (xu.x-xl.x) * (xu.y-xl.y) * (xu.z-xl.z)  for point in p_2 ]
							
						append!(pt,pt_2)
						append!(f,f_2)

					else
						continue
					end
					
				end
			end
		end
		#println("CONTROLLO: ",control)
		Integral = sum(f) #/ ( cal / 3.)
		Integral_abs = sum(abs.(f)) #/ ( cal / 3.  )

		g = (abs.(f) ./ Integral_abs) #/ ( cal / 3. )
		w_m     = zeros(Int64(NSIDE^3))
	
	
		println("Stima *Err*\t \t = ", ( Integral / (control) ) - (Int_Vero * 3. / cal)    ) 	
		println(length(f))
		println(length(pt) / control)  

		
	end
		

	
	return Integral / (control)
end


#################### Integration #####################

calls = 1000000

xu = Point{Float64}(1., 1., 1.)
xl = Point{Float64}(-1.,-1.,-1.)



result = Integrate(gauss_3D, xl, xu, calls)

println("Valore dell'integrale = ", result)




