struct Point{T}
	r::T
	θ::T
	ϕ::T
	x
 	y
	z
end


function(p::Point)()
	x = p.r * cos(p.θ) * cos(p.ϕ)
	y = p.r * cos(p.θ) * sin(p.ϕ)
	z = p.r * sin(p.θ)
	return Point(p.r, p.θ, p.ϕ, x, y, z)
end

function get_versor(p::Point{Float64})
	return Point{Float64}(1, p.θ, p.ϕ, p.x / p.r, p.y / p.r, p.z / p.r)
end

function diff_cart(p::Point{Float64}, p_2::Point{Float64})
	return Point{Float64}(0, 0, 0, p.x-p_2.x, p.y-p_2.y, p.z-p_2.z)
end




Base.:+(a::Point{Float64}, b::Point{Float64}) = Point{Float64}(a.r+b.r, a.θ+b.θ, a.ϕ+b.ϕ, undef, undef, undef)()
Base.:-(a::Point{Float64}, b::Point{Float64}) = Point{Float64}(a.r-b.r, a.θ-b.θ, a.ϕ-b.ϕ, undef, undef, undef)()
Base.:*(a::Point{Float64}, b::Point{Float64}) = Point{Float64}(a.r*b.r, a.θ*b.θ, a.ϕ*b.ϕ, undef, undef, undef)()

# Non e` una vera divisione e` quello che mi serve per ottenere il versore partendo dal vettore
# Forse e` meglio fare una funzione del tipo get_versor -> a.r / a.r , a.θ, a.ϕ, a.x / a.r, a.y /a.r , a.z / a.r
# in questo modo si ottiene direttamente il versore nella direzione di a
