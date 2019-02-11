# Test Function
function gauss_3D(p::Point{Float64}, p_p::Point{Float64})
	return @. exp(-(p.x)^ 2 - (p.y)^ 2 - (p.z)^ 2 - (p_p.x)^2 - (p_p.y)^2 - (p_p.z)^2)
end

# Jacobian of the translation
function jac_trasl(xl::Point{Float64}, xu::Point{Float64}, xl_p::Point{Float64}, xu_p::Point{Float64})
	return (xu.r-xl.r) * (xu.θ-xl.θ) * (xu.ϕ-xl.ϕ) * (xu_p.r-xl_p.r) * (xu_p.θ-xl_p.θ) * (xu_p.ϕ-xl_p.ϕ)
end

# A naif scalar product between two Point
function inner_prod(p_one::Point{Float64}, p_two::Point{Float64})
	return (p_one.x*p_two.x + p_one.y*p_two.y + p_one.z*p_two.z)
end

# The BEAM WAIST
function ω(p::Point{Float64})
	return (ω₀) * (1+(λ*inner_prod(get_versor(p), p) / (π * ω₀^2)  )^2)^0.5
end

# Atmospheric vertical temperature profile
function Tₚ(p::Point{Float64})
	return Tᵧ * ( 1 - (p.z / z_atm) )
end

# Distribution of turbulent structures on the atmosphere
function χ₁(p::Point{Float64}, p_p::Point{ Float64}, t_s::Float64)
	W        = Point{Float64}(45, 0, (π/2)+(π/(7.0)), undef, undef, undef)()
	W_incr   = Point{Float64}(0.0, 0.0, 0.0, W.x*t_s, W.y*t_s, 0)
	dr = diff_cart(p, p_p)
	dr_mod = ((dr.x-W_incr.x)^2 + (dr.y-W_incr.y)^2 + dr.z^2)^0.5
	return χ₁_0 * exp(- ( (dr_mod)^2 / (2*(L₀)^2) ) )
end

# Water Vapor vertical density profile
function χ₂(p::Point{Float64}, p_p::Point{Float64})
	return χ₂_0 * exp(- (p.z+p_p.z) / (2*z_0) )
end


# Effective area assiaciated to the beam projection in the sky
function B(p::Point{Float64})
	r_mod = (p.x^2 + p.y^2 +p.z^2)^0.5
	first_part = ( 2 * λ^2 * abs(inner_prod(get_versor(p), p))^2 ) / (π * ( ω(p) )^2 )
	#second_part = exp( - (r_mod^2 - (inner_prod(get_versor(p), p))^2) / (ω(p)^2 ) )
	return first_part #* second_part
end

# Jacobiano coordinate sferiche
function Jac_Sp(p::Point{Float64})
	return  cos(p.θ)
end

# Integrand of correlation coef

function corr(p::Point{Float64}, p_p::Point{Float64}, time_stamp::Float64)
 	return B(p) * B(p_p) * χ₁(p, p_p, time_stamp) * χ₂(p, p_p) * Tₚ(p) * Tₚ(p_p) * Jac_Sp(p) * Jac_Sp(p_p)
end
