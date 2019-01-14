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
	return @. exp(-(p.x)^ 2 - (p.y)^ 2 - (p.z)^ 2)
end


function Integrate(fun::Function, xl::Point{Float64}, xu::Point{Float64}, cal::Int64)
	
	#; Primo step, cal/3. valutazioni della funzione
	rng = MersenneTwister(0)
	pt = [Point{Float64}(rand(rng), rand(rng), rand(rng)) for n = 1:Int64(floor(calls/3.))]

	#; Cambio coordinate, per mettermi nel range di integrazione scelto
	p = ( pt .* Ref(xu-xl) ) .+ Ref(xl)
	density_unif = 1 / Int64(floor(calls/3.))

	#; Valuto l'integrale
	f = [fun(point) * (xu.x-xl.x) * (xu.y-xl.y) * (xu.z-xl.z) * density_unif  for point in p ]
	f_abs = abs.(f)
	
	Integral      = sum(f)
	Integral_abs  = sum(f_abs) 

	#; Valore dell'integrale di prova
	Int_Vero = 3.33231
	println("0. Integrale = ",(Integral-Int_Vero)," ", Integral)
	
	#; Calcolo la nuova densita` di probabilita`
	g = (f_abs ./ Integral_abs)
	
	#println("Controllo integrale pdf = ",sum(g))
	print("Insert bin = " )
	NSIDE = 14  #Int64(floor(1+log(2,cal/3.) ))# + log(2, 1+skewness(g)/((6*(cal/3. -2)/((cal/3. + 1)*(cal/3. + 2)))^0.5)  )) )
	println(NSIDE)

	#; posso fare qui la separazione dei punti nelle varie BOX
	
	w_m	 = zeros(Int64(NSIDE^3))
	Int_box  = zeros(Int64(NSIDE^3))
	Int_abs_box = zeros(Int64(NSIDE^3))
	dix = diy = diz = 1.0/NSIDE
	for iter = 1:20
		#; Organizzo la pdf nel BOX
		for idx_g = 1:length(g)
			box_id_x = Int64(ceil(pt[idx_g].x / dix)) -1
			box_id_y = Int64(ceil(pt[idx_g].y / diy)) -1
			box_id_z = Int64(ceil(pt[idx_g].z / diz)) -1
			ID = (box_id_x) + (box_id_y)*NSIDE + (box_id_z)*NSIDE*NSIDE
			w_m[ID+1] += g[idx_g]  
		end
	
		#println("Densita` in BOX")
		#for i in 1:length(w_m)
			#println(w_m[i])
		#end
		#println("Controllo Integrale PDF = ", sum(w_m))
		
		control=0
		pt_2_all = Array{Point{Float64}}(undef, 0)
		f_2_all  = Array{Float64,1}(undef, 0)
		
		for ix = (1:NSIDE)
			for iy = (1:NSIDE)
				for iz = (1:NSIDE)
					p_start = Point{Float64}((ix-1) * dix , (iy-1) * diy, (iz-1) * diz)
					p_final = Point{Float64}(ix * dix, iy * diy, iz * diz)

					ID_BOX = ( (ix-1) + ( (iy-1) * NSIDE ) + ( (iz-1) * NSIDE * NSIDE ) ) + 1
					
					#; Densita` Uniforme
					#density    = 1/(NSIDE*NSIDE*NSIDE)

					#; Densita` Pesata
					density    = w_m[ID_BOX]


					hit_in_box = Int64(floor(density * (cal/3.)))
					control    += hit_in_box
					#println("Hit = ", hit_in_box)
					if hit_in_box != 0
						#; Qui va costruita la nuova densita` di punti
						pt_2 = [Point{Float64}((p_final.x-p_start.x)*rand(rng)+p_start.x, (p_final.y-p_start.y)*rand(rng)+p_start.y,(p_final.z-p_start.z)*rand(rng)+p_start.z) for n = 1:hit_in_box]
						#; Cambio coordinate
						p_2 = ( pt_2 .* Ref(xu-xl) ) .+ Ref(xl)
						
						f_2= [fun(point) * (xu.x-xl.x) * (xu.y-xl.y) * (xu.z-xl.z) / density  for point in p_2 ]
							
						append!(pt_2_all, pt_2)
						append!(f_2_all, f_2)

					else
						continue
					end
					
				end
			end
		end
		
		Integral += sum(f_2_all) / (control * NSIDE * NSIDE * NSIDE)
		Integral_abs = sum(abs.(f_2_all))
		g = (abs.(f_2_all) ./ Integral_abs) 
		w_m     = zeros(Int64(NSIDE^3))
	
	
		println(iter,". Integrale = ", (Integral / (iter+1)) - Int_Vero, " ", Integral / (iter + 1)  ) 	
		
	end
		

	
	return Integral
end


#################### Integration #####################

calls = 1000000

xu = Point{Float64}(1., 1., 1.)
xl = Point{Float64}(-1.,-1.,-1.)



result = Integrate(gauss_3D, xl, xu, calls)





