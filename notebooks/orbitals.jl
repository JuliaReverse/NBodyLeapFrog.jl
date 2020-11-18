### A Pluto.jl notebook ###
# v0.12.10

using Markdown
using InteractiveUtils

# ╔═╡ 83a2c2cc-296b-11eb-38bf-fd97aa3c6fda
using PlutoUI, NBodyLeapFrog, Plots, NiLang

# ╔═╡ fdafd8a6-29d0-11eb-0d05-73cb4a836a20
using NiLang.AD

# ╔═╡ c08d80cc-29d0-11eb-002d-d1e08da6a7c1
md"# Simulating Celestial Bodies"

# ╔═╡ 814c93d0-296c-11eb-333d-216ce303ce34
plotly();

# ╔═╡ bd337dc6-29d0-11eb-1a58-b354615e11e1
md"## Computing the orbitals with time-reversible Forest-Ruth integrator"

# ╔═╡ 8e419bc2-296b-11eb-0e0f-afdab78b1392
function showres(vs; plt=plot([], []))
    a = map(x->x.x, vs)
    b = map(x->x.y, vs)
	labels = ["sun", "mercury", "venus", "earth", "mars",
		"jupyter", "saturn", "uranus", "neptune", "pluto"]
    for i in 1:size(a, 1)
        plot!(plt,a[i,:],b[i,:]; aspect_ratio=:equal, label=labels[i])
    end
    plt
end

# ╔═╡ 4afeee16-29d4-11eb-0a02-ff725633872c
planets = Bodies.chunit_day2year.(Bodies.set)

# ╔═╡ a4110f0e-29d0-11eb-1cf8-a950830d0cfa
md"## Differentiate the simulation"

# ╔═╡ bbe14adc-296b-11eb-1cb2-c149ac6e6321
@i function loss!(y::T, r0::AbstractVector{V3{T}}, v0::AbstractVector{V3{T}}, m; n, dt) where T
	for i=1:n
		NBodyLeapFrog.i_fr_step!(r0, v0, m; dt=dt)
	end
	
	# compute the squared distance between Mercury (2) and earth (4)
	@routine begin
		sqy ← zero(T)
		sqy += NBodyLeapFrog.sqdistance(r0[9], r0[10])
	end
	y += sqrt(sqy)

	~@routine # uncompute sqy
end

# ╔═╡ de67315c-296b-11eb-0ae2-93757c2e2331
r0, v0, m = let
	planets = Bodies.chunit_day2year.(Bodies.set)
	getfield.(planets, :r), getfield.(planets, :v), getfield.(planets, :m)
end;

# ╔═╡ 36acacfe-29ce-11eb-18b3-0fc163ea0dd5
"""simulate and keep history"""
@i function simulate!(rmat::AbstractMatrix{V3{T}}, vmat::AbstractMatrix{V3{T}}, planets; dt) where T
	@safe @assert size(rmat) == size(vmat) && size(rmat, 1) == length(planets)
	
	# store the mass into a vector
	@routine begin
		m ← zeros(T, length(planets))
		for j=1:length(planets)
			m[j] += planets[j].m
		end
	end
	# setup initial values
	for j=1:length(planets)
		rmat[j,1] += planets[j].r
		vmat[j,1] += planets[j].v
	end

	# simulate!
	for i=1:size(rmat, 2)-1
		for j=1:length(planets)
			rmat[j,i+1] += rmat[j,i]
			vmat[j,i+1] += vmat[j,i]
		end
		# single step Forest-Ruth algorithm
		@inbounds NBodyLeapFrog.i_fr_step!(view(rmat,:,i+1), view(vmat,:,i+1), m; dt=dt)
	end
	
	~@routine #uncompute m
end

# ╔═╡ 157f8d8e-29cf-11eb-3926-01d791d6d703
rmat, vmat = let
	rmat = zeros(V3{Float64}, length(planets), 4001)
	vmat = zeros(V3{Float64}, length(planets), 4001)
	simulate!(rmat, vmat, planets; dt=0.01)
end;

# ╔═╡ 6529ba3a-29d4-11eb-1123-59951753b134
showres(rmat)

# ╔═╡ 6b217206-296c-11eb-0e56-a742db172508
y, r_, v_, m_ = loss!(0.0, copy(r0), copy(v0), copy(m); dt=0.01, n=4000);

# ╔═╡ ba521b86-296b-11eb-2bd4-ebf29f78b87a
_, _, gr, gv, gm = Grad(loss!)(Val(1), 0.0, copy(r0), copy(v0), copy(m); dt=0.01, n=4000); # Val(1) means, the the first argument stores the loss

# ╔═╡ 3863e6ea-29d1-11eb-3d9c-8baf25f581ca
gr

# ╔═╡ 0fa5ad5a-29d2-11eb-258f-590d950d565b
function showgrad(r0, gr)
    x0 = map(x->x.x, r0)
    y0 = map(x->x.y, r0)
    gx0 = map(x->x.x.g, gr)
    gy0 = map(x->x.y.g, gr)
	labels = ["sun", "mercury", "venus", "earth", "mars",
		"jupyter", "saturn", "uranus", "neptune", "pluto"]
    quiver(x0, y0; quiver=collect(zip(gx0, gy0)))
end

# ╔═╡ 88441af8-29d2-11eb-3e77-5bace9a3a1cf
let
	plt = showgrad(r0, gr)
	showres(rmat, plt=plt)
end

# ╔═╡ Cell order:
# ╟─c08d80cc-29d0-11eb-002d-d1e08da6a7c1
# ╠═83a2c2cc-296b-11eb-38bf-fd97aa3c6fda
# ╠═814c93d0-296c-11eb-333d-216ce303ce34
# ╟─bd337dc6-29d0-11eb-1a58-b354615e11e1
# ╠═8e419bc2-296b-11eb-0e0f-afdab78b1392
# ╠═36acacfe-29ce-11eb-18b3-0fc163ea0dd5
# ╠═4afeee16-29d4-11eb-0a02-ff725633872c
# ╠═157f8d8e-29cf-11eb-3926-01d791d6d703
# ╠═6529ba3a-29d4-11eb-1123-59951753b134
# ╟─a4110f0e-29d0-11eb-1cf8-a950830d0cfa
# ╠═fdafd8a6-29d0-11eb-0d05-73cb4a836a20
# ╠═bbe14adc-296b-11eb-1cb2-c149ac6e6321
# ╠═de67315c-296b-11eb-0ae2-93757c2e2331
# ╠═6b217206-296c-11eb-0e56-a742db172508
# ╠═ba521b86-296b-11eb-2bd4-ebf29f78b87a
# ╠═3863e6ea-29d1-11eb-3d9c-8baf25f581ca
# ╠═0fa5ad5a-29d2-11eb-258f-590d950d565b
# ╠═88441af8-29d2-11eb-3e77-5bace9a3a1cf
