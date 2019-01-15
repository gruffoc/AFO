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


function Integrate(fun::Function, xl::Point{Float64}, xu::Point{Float64}, cal::Int64, dim::Float64)

	# First Step: cal evaluations of the integrand using uniform distributed random number
	rng = MersenneTwister(0)
	pt = [Point{Float64}(rand(rng), rand(rng), rand(rng)) for n = 1:Int64(floor(dim * cal))]

	# Coordinate transformation. ∫ₐᵇ f(x) dx  = ∫₀¹ f(y *(b-a) +a )*(b-a) dy
	p = ( pt .* Ref(xu-xl) ) .+ Ref(xl)
	density_unif = 1 / Int64(floor(dim * calls))

	# Integrand evaluations
	f = [fun(point) * (xu.x-xl.x) * (xu.y-xl.y) * (xu.z-xl.z) * density_unif  for point in p ]
	f_abs = abs.(f)

	# Perform the summation of intergrand and its absolute value
	Integral      = sum(f)
	Integral_abs  = sum(f_abs)

	# Check value - \Int_{-1}^{1} gauss_3D(x) dx
	Int_True = 3.33231
	#println("0. Integral = ",(Integral-Int_True)," ", Integral)

	# Estimation of the new weight function
	g = (f_abs ./ Integral_abs)

	# println("Check PDF Integral = ",sum(g))
	# I have to study hard to find a more suitable method to estimate the binning
	NBIN = 15

	# PDF approximation with stepper function
	w_m	 = zeros(Int64(NBIN^3))
	Int_box  = zeros(Int64(NBIN^3))
	Int_abs_box = zeros(Int64(NBIN^3))
	dix = diy = diz = 1.0/NBIN
	Integral_avg = 0
	for iter = 1:20
		# Evaluate the weigth on each BOX
		for idx_g = 1:length(g)
			box_id_x = Int64(ceil(pt[idx_g].x / dix)) -1
			box_id_y = Int64(ceil(pt[idx_g].y / diy)) -1
			box_id_z = Int64(ceil(pt[idx_g].z / diz)) -1
			ID = (box_id_x) + (box_id_y)*NBIN + (box_id_z)*NBIN*NBIN
			w_m[ID+1] += g[idx_g]
		end

		# Check the weights distribution on each box and check if the
		# PDF requirements is satisfied \Int_{\Omega} pdf(x) dx = 1

		#println("Densita` in BOX")
		#for i in 1:length(w_m)
			#println(w_m[i])
		#end
		#println("Check PDF Integral = ", sum(w_m))

		# Evaluate the Integral using a non uniform random number distribution
		total_hits = 0
		pt_2_all   = Array{Point{Float64}}(undef, 0)
		f_2_all    = Array{Float64,1}(undef, 0)

		for ix = (1:NBIN)
			for iy = (1:NBIN)
				for iz = (1:NBIN)
					p_start = Point{Float64}((ix-1) * dix , (iy-1) * diy, (iz-1) * diz)
					p_final = Point{Float64}(ix * dix, iy * diy, iz * diz)

					ID_BOX = ( (ix-1) + ( (iy-1) * NBIN ) + ( (iz-1) * NBIN * NBIN ) ) + 1

					# Uniform density
					#density    = 1/(NBIN*NBIN*NBIN)

					# Non uniform weight PDF
					density    = w_m[ID_BOX]


					hit_in_box = Int64(floor(density * (dim * cal)))
					total_hits += hit_in_box
					#println("Hit = ", hit_in_box)
					if hit_in_box != 0

						pt_2 = [Point{Float64}((p_final.x-p_start.x)*rand(rng)+p_start.x, (p_final.y-p_start.y)*rand(rng)+p_start.y,(p_final.z-p_start.z)*rand(rng)+p_start.z) for n = 1:hit_in_box]
						p_2 = ( pt_2 .* Ref(xu-xl) ) .+ Ref(xl)
						# Aggiustare quel /density!!!!!! Non e` conforme con quello sopra!!!!
						f_2= [fun(point) * (xu.x-xl.x) * (xu.y-xl.y) * (xu.z-xl.z) / density  for point in p_2 ]

						append!(pt_2_all, pt_2)
						append!(f_2_all, f_2)

					else
						continue
					end

				end
			end
		end

		Integral += sum(f_2_all) / (total_hits * NBIN * NBIN * NBIN)
		Integral_abs = sum(abs.(f_2_all))
		g = (abs.(f_2_all) ./ Integral_abs)
		w_m     = zeros(Int64(NBIN^3))


		println(iter,". Integral = ", (Integral / (iter+1)) - Int_True, " ", Integral / (iter + 1)  )
		Integral_avg = Integral / (iter+1)


		# check the point distribution in the cube
		#plotly()
		#x = [i.x for i in pt_2_all]
		#y = [i.y for i in pt_2_all]
		#z = [i.z for i in pt_2_all]
		#scatter3d(x,y,z, markersize = 0.2, show = true, dpi=200, size=(1024,768))

	end



	return Integral_avg
end


#################### Integration #####################

calls = 1000000
dimensions = 3.

xu = Point{Float64}(1., 1., 1.)
xl = Point{Float64}(-1.,-1.,-1.)



result = Integrate(gauss_3D, xl, xu, calls, dimensions)

println("Final Integral = ",result)
