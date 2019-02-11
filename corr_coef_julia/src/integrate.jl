function Integrate(fun::Function, xl::Point{Float64}, xu::Point{Float64}, xl_p::Point{Float64}, xu_p::Point{Float64}, cal::Int64, iterations::Int64, dim::Float64, time_stamp::Float64)

	# First Step: cal evaluations of the integrand using uniform distributed random number
	rng = MersenneTwister(0)
	#la dimensione non va bene
	pt = [Point{Float64}(rand(rng), rand(rng), rand(rng), undef, undef, undef)() for n = 1:Int64(floor(3. * cal))]
	pt_p = [Point{Float64}(rand(rng), rand(rng), rand(rng), undef, undef, undef)() for n = 1:Int64(floor(3. * cal))]

	#Deve farmi la differenza degli r, theta e phi e calcolarmi le x,y,z di conseguenza
	# Coordinate transformation. ∫ₐᵇ f(x) dx  = ∫₀¹ f(y *(b-a) +a )*(b-a) dy
	p   = ( pt .* Ref(xu-xl) ) .+ Ref(xl)
	p_p = ( pt_p .* Ref(xu_p -xl_p) ) .+ Ref(xl)
	# Jacobian of the translation
	J_t = jac_trasl(xl, xu, xl_p, xu_p)

	density_unif = 1. / (Int64(floor(3. * calls)))

	# Integrand evaluations
	f = [fun(p[i], p_p[i], time_stamp) * J_t * density_unif  for i = 1:Int64(floor(3. * cal)) ]
	f_abs = abs.(f)

	# Perform the summation of intergrand and its absolute value
	Integral      = sum(f)
	Integral_abs  = sum(f_abs)

	# Check value - \Int_{-1}^{1} gauss_3D(x) dx
	#Int_True = 0.00378269
	#println("0. Integral = ", Integral)

	# Estimation of the new weight function
	g = (f_abs ./ Integral_abs)

	#println("Check PDF Integral = ",sum(g))
	# I have to study hard to find a more suitable method to estimate the binning
	NBIN = 3

	# PDF approximation with stepper function
	w_m	 = zeros(Int64(NBIN^dim))
	Int_box  = zeros(Int64(NBIN^dim))
	Int_abs_box = zeros(Int64(NBIN^dim))
	dix = diy = diz = div = diw = dik = 1.0/NBIN
	Integral_avg = 0



	# Box_R = Array{Float64}(undef, 0)
	# Box_T = Array{Float64}(undef, 0)
	# Box_P = Array{Float64}(undef, 0)
	#
	# Box_R_P = Array{Float64}(undef, 0)
	# Box_T_P = Array{Float64}(undef, 0)
	# Box_P_P = Array{Float64}(undef, 0)

	# HIT = Array{Float64}(undef, 0)

	chi_square = 0


	for iter = 1:iterations
		# Evaluate the weigth on each BOX
		for idx_g = 1:length(g)
			box_id_x = Int64(ceil(pt[idx_g].r / dix)) -1
			box_id_y = Int64(ceil(pt[idx_g].θ / diy)) -1
			box_id_z = Int64(ceil(pt[idx_g].ϕ / diz)) -1

			box_id_v = Int64(ceil(pt_p[idx_g].r / dix)) -1
			box_id_w = Int64(ceil(pt_p[idx_g].θ / diy)) -1
			box_id_k = Int64(ceil(pt_p[idx_g].ϕ / diz)) -1
			# Qui si puo` forse scrivere meglio
			ID = (box_id_x) + (box_id_y)*NBIN + (box_id_z)*NBIN*NBIN + (box_id_v)*NBIN*NBIN*NBIN + (box_id_w)*NBIN*NBIN*NBIN*NBIN + (box_id_k)*NBIN*NBIN*NBIN*NBIN*NBIN


			w_m[ID+1] += g[idx_g]
		end

		# Fino a qui mi sembra tutto in ordine. Gli indici mi tornano tutti.

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
		pt_2_p_all = Array{Point{Float64}}(undef, 0)
		f_2_all    = Array{Float64,1}(undef, 0)

		for ix = (1:NBIN)
			for iy = (1:NBIN)
				for iz = (1:NBIN)
					for iv = (1:NBIN)
						for iw = (1:NBIN)
							for ik = (1:NBIN)
								p_start = Point{Float64}((ix-1) * dix , (iy-1) * diy, (iz-1) * diz, undef, undef, undef)
								p_final = Point{Float64}(ix * dix, iy * diy, iz * diz, undef, undef, undef)
								p_p_start = Point{Float64}((iv-1) * div , (iw-1) * diw, (ik-1) * dik, undef, undef, undef)
								p_p_final = Point{Float64}(iv * div, iw * diw, ik * dik, undef, undef, undef)

								# Si puo` semplificare la scrittura?
								ID_BOX = ( (ix-1) + ( (iy-1) * NBIN ) + ( (iz-1) * NBIN * NBIN ) +  ( (iv-1) * NBIN * NBIN * NBIN ) + ( (iw-1) * NBIN * NBIN * NBIN * NBIN ) + ( (ik-1) * NBIN * NBIN * NBIN * NBIN * NBIN)    ) + 1

								# Uniform density
								#density    = 1/(NBIN*NBIN*NBIN)

								# Non uniform weight PDF
								density    = w_m[ID_BOX]


								hit_in_box = Int64(floor(density * (3. * cal)))

								total_hits += hit_in_box
								#println("Hit = ", hit_in_box)
								if hit_in_box != 0


									x_ran = (p_final.r-p_start.r)
									y_ran = (p_final.θ-p_start.θ)
									z_ran = (p_final.ϕ-p_start.ϕ)

									v_ran = (p_p_final.r-p_p_start.r)
									w_ran = (p_p_final.θ-p_p_start.θ)
									k_ran = (p_p_final.ϕ-p_p_start.ϕ)

									pt_2 = [Point{Float64}(x_ran*rand(rng)+p_start.r, y_ran*rand(rng)+p_start.θ, z_ran*rand(rng)+p_start.ϕ, undef,undef,undef)() for n = 1:hit_in_box]
									pt_p_2 = [Point{Float64}(v_ran*rand(rng)+p_p_start.r, w_ran*rand(rng)+p_p_start.θ, k_ran*rand(rng)+p_p_start.ϕ,undef,undef,undef)() for n = 1:hit_in_box]


									p_2 = ( pt_2 .* Ref(xu-xl) ) .+ Ref(xl)
									p_p_2 = ( pt_p_2 .* Ref(xu_p-xl_p) ) .+ Ref(xl_p)
									# Aggiustare quel /density!!!!!! Non e` conforme con quello sopra!!!!
									f_2 = [fun(p_2[i], p_p_2[i], time_stamp) * J_t / density  for i = 1:hit_in_box ]


									append!(pt_2_all, pt_2)
									append!(pt_2_p_all, pt_p_2)

									append!(f_2_all, f_2)

									# append!(Box_R, ix)
									# append!(Box_T, iy)
									# append!(Box_P, iz)
									#
									# append!(Box_R_P, iv)
									# append!(Box_T_P, iw)
									# append!(Box_P_P, ik)
									# append!(HIT,hit_in_box)


								else
									continue
								end
							end
						end
					end
				end
			end
		end
		Integral += sum(f_2_all) / (total_hits * NBIN^6 )
		Integral_abs = sum(abs.(f_2_all))
		g = (abs.(f_2_all) ./ Integral_abs)
		w_m = zeros(Int64(NBIN^6))

		pt = pt_2_all
		pt_p = pt_2_p_all

		#println(iter,". Integral = ", Integral / (iter + 1)  )
		Integral_avg = Integral / (iter+1)
		chi_square += (sum(f_2_all) / (total_hits * NBIN^6 ) - Integral_avg)^2


		# println("#######################################################################")
		# println( "r_min = ",     ((findmin(Box_R)[1]-1)/NBIN * (xu.r-xl.r)) + xl.r   )
		# println( "r_max = ",    ((findmax(Box_R)[1])/NBIN * (xu.r-xl.r)) + xl.r  )
		#
		# println( "θ_min = ",   ((findmin(Box_T)[1]-1)/NBIN * (xu.θ-xl.θ)) + xl.θ   )
		# println( "θ_max = ",    ((findmax(Box_T)[1])/NBIN * (xu.θ-xl.θ)) + xl.θ  )
		#
		# println( "ϕ_min = ",    ((findmin(Box_P)[1]-1)/NBIN * (xu.ϕ-xl.ϕ)) + xl.ϕ   )
		# println( "ϕ_min = ",   ((findmax(Box_P)[1])/NBIN * (xu.ϕ-xl.ϕ)) + xl.ϕ  )
		#
		# println("################################## PRIMED #############################")
		#
		# println( "r_min = ",   ((findmin(Box_R_P)[1]-1)/NBIN * (xu_p.r-xl_p.r)) + xl_p.r   )
		# println( "r_max = ",    ((findmax(Box_R_P)[1])/NBIN * (xu_p.r-xl_p.r)) + xl_p.r  )
		#
		# println( "θ_min = ",   ((findmin(Box_T_P)[1]-1)/NBIN * (xu_p.θ-xl_p.θ)) + xl_p.θ   )
		# println( "θ_max = ",   ((findmax(Box_T_P)[1])/NBIN * (xu_p.θ-xl_p.θ)) + xl_p.θ  )
		#
		# println( "ϕ_min = ",   ((findmin(Box_P_P)[1]-1)/NBIN * (xu_p.ϕ-xl_p.ϕ)) + xl_p.ϕ   )
		# println( "ϕ_max = ",   ((findmax(Box_P_P)[1])/NBIN * (xu_p.ϕ-xl_p.ϕ)) + xl_p.ϕ  )
		# println("#######################################################################")
		#
		#
		# Box_R = Array{Float64}(undef, 0)
		# Box_T = Array{Float64}(undef, 0)
		# Box_P = Array{Float64}(undef, 0)
		#
		# Box_R_P = Array{Float64}(undef, 0)
		# Box_T_P = Array{Float64}(undef, 0)
		# Box_P_P = Array{Float64}(undef, 0)



	end



	return Integral_avg, chi_square #, HIT
end
