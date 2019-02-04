using Random
using Statistics
using StatsBase
using Plots
import Base.+
import Base.-



struct Point{T}
	x::T
	y::T
end

Base.:+(a::Point{Float64}, b::Point{Float64}) = Point{Float64}(a.x+b.x , a.y+b.y)
Base.:-(a::Point{Float64}, b::Point{Float64}) = Point{Float64}(a.x-b.x , a.y-b.y)
Base.:*(a::Point{Float64}, b::Point{Float64}) = Point{Float64}(a.x*b.x , a.y*b.y)



function gauss_3D(p::Point{Float64})
	return exp( -(p.x)^2 - (p.y)^2 )
end

function gradino(p::Point{Float64})
	if (p.x > -1 && p.x < 1 && p.y > -1 && p.y < 1)
		return 1.0
	else
		return 0.0
	end
end

function cilindro(p::Point{Float64})
	if (p.x^2 + p.y^2 < 1 )
		return 1.0
	else
		return 0.0
	end
end


function Integrate(fun::Function, xl::Point{Float64}, xu::Point{Float64}, cal::Int64)

	#; Primo step, cal/3. valutazioni della funzione
	rng = MersenneTwister(0)
	pt = [Point{Float64}(rand(rng), rand(rng)) for n = 1:Int64(floor(calls/2.))]

	#; Cambio coordinate, per mettermi nel range di integrazione scelto
	p = ( pt .* Ref(xu-xl) ) .+ Ref(xl)
	density_unif = 1 / Int64(floor(calls/2.))

	#; Valuto l'integrale
	f = [ fun(point) * (xu.x-xl.x) * (xu.y-xl.y) * density_unif  for point in p ]
	f_abs = abs.(f)

	Integral      = sum(f)
	Integral_abs  = sum(f_abs)

	#; Valore dell'integrale di prova
	Int_Vero = 2.23099
	println("0. Integrale = ",(Integral - Int_Vero)," ", Integral)

	#; Calcolo la nuova densita` di probabilita`
	g = (f_abs ./ Integral_abs)

	#for i = 1:length(g)
	#	println(pt[i].x," ",pt[i].y," ",g[i])
	#end

	#println("1. Controllo integrale pdf = ",sum(g))

	#println(1+log(2,cal/2.) + log(2, 1+skewness(g)/((6*(cal/2. -2)/((cal/2. + 1)*(cal/2. + 2)))^0.5)  ) )

	NSIDE = 1000

	#; posso fare qui la separazione dei punti nelle varie BOX

	w_m	 = zeros(Int64(NSIDE^2))
	Int_box  = zeros(Int64(NSIDE^2))
	Int_abs_box = zeros(Int64(NSIDE^2))
	dix = diy = diz = 1.0/NSIDE
	xx = Array{Float64}(undef, 0)
	yy = Array{Float64}(undef, 0)

	Box_X = Array{Float64}(undef, 0)
	Box_Y = Array{Float64}(undef, 0)


	for iter = 1:20

		#; Organizzo la pdf nel BOX
		for idx_g = 1:length(g)
			box_id_x = Int64(ceil(pt[idx_g].x / dix)) -1
			box_id_y = Int64(ceil(pt[idx_g].y / diy)) -1
			ID = (box_id_x) + (box_id_y)*NSIDE
			w_m[ID+1] += g[idx_g]
		end

		#println("Densita` in BOX")
		#for i in 1:length(w_m)
		#	println(w_m[i])
		#end
		#println(iter+1,". ","Controllo Integrale PDF = ", sum(w_m))
		control = 0
		pt_2_all = Array{Point{Float64}}(undef, 0)
		f_2_all  = Array{Float64}(undef, 0)
		p_2 = Array{Point{Float64}}(undef, 0)

		for ix = (1:NSIDE)
			for iy = (1:NSIDE)
				p_start = Point{Float64}((ix-1) * dix , (iy-1) * diy)
				p_final = Point{Float64}(ix * dix, iy * diy)
				ID_BOX = ( (ix-1) + ( (iy-1) * NSIDE ) ) + 1

				#; Densita` Uniforme
				#density    = 1/(NSIDE*NSIDE)

				#; Densita` Pesata
				density    = w_m[ID_BOX]


				hit_in_box = Int64(floor(density * (cal/2.)))
				control    += hit_in_box
				#println("Hit = ", hit_in_box)
				if hit_in_box != 0
					#; Qui va costruita la nuova densita` di punti
					pt_2 = [Point{Float64}((p_final.x-p_start.x)*rand(rng)+p_start.x, (p_final.y-p_start.y)*rand(rng)+p_start.y) for n = 1:hit_in_box]

					#; Cambio Coordinate
					p_2 = ( pt_2 .* Ref(xu-xl) ) .+ Ref(xl)


					#; Valuto l'integrale
					f_2 = [fun(point) * (xu.x-xl.x) * (xu.y-xl.y) / density   for point in p_2 ]

					#println(p_start.x, " ", p_final.x, " ", p_start.y," ", p_final.y, " ", sum(f_2)/hit_in_box)

					append!(pt_2_all, pt_2)
					append!(f_2_all, f_2)

					append!(Box_X, ix)
					append!(Box_Y, iy)

					#println("Box = ", ix, " ", iy)



				else
					continue
				end
			end
		end

		# xu = Point{Float64}(x_max, y_max)
		# xl = Point{Float64}(x_min, y_min)

		# println( xmin, " ", xmax, " ", ymin, " ", ymax )

		pt = pt_2_all



		Integral += (sum(f_2_all) / (control*NSIDE*NSIDE) )
		Integral_abs = sum(abs.(f_2_all))
		g = (abs.(f_2_all) ./ Integral_abs)
		w_m     = zeros(Int64(NSIDE^2))

		println(iter,". Integrale = ", (Integral / (iter + 1) ) - Int_Vero, " ", Integral/(iter + 1) )

		prova = (pt_2_all .* Ref(xu-xl) ) .+ Ref(xl)

		for i = 1:length(g)
			xx=append!(xx, prova[i].x)
			yy=append!(yy, prova[i].y)
		end



		# println("Per le X")
		# println(      ((findmin(Box_X)[1]-1)/100 * (xu.x-xl.x)) + xl.x   )
		# println(      ((findmax(Box_X)[1])/100 * (xu.x-xl.x)) + xl.x  )
		#
		# println("Per le Y")
		# println(((findmin(Box_Y)[1]-1)/100 * (xu.x-xl.x)) + xl.x)
		# println(((findmax(Box_Y)[1])/100 * (xu.x-xl.x)) + xl.x )
		#
		# # xu_n = Point{Float64}(((findmax(Box_X)[1])/100 * (xu.x-xl.x)) + xl.x, ((findmax(Box_Y)[1])/100 * (xu.x-xl.x)) + xl.x     )
		# # xl_n = Point{Float64}(((findmin(Box_X)[1])/100 * (xu.x-xl.x)) + xl.x, ((findmin(Box_Y)[1])/100 * (xu.x-xl.x)) + xl.x     )
		# #
		# # xu = xu_n
		# # xl = xl_n
		#
		# println(xu)
		# Box_X = Array{Float64}(undef, 0)
		# Box_Y = Array{Float64}(undef, 0)


	end #end iteration

	return Integral,xx,yy
end


#################### Integration #####################

calls = 100000

xu = Point{Float64}(1., 1.)
xl = Point{Float64}(-1.,-1.)



result,x,y = Integrate(cilindro, xl, xu, calls)
plot(x,y, seriestype=:scatter)

#println("Valore dell'integrale = ", result)


# il punto e` che gli estremi di integrazione devono essere rimodellati
# dove la funzione non e` diversa da zero!
