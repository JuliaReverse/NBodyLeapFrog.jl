### A Pluto.jl notebook ###
# v0.12.10

using Markdown
using InteractiveUtils

# ╔═╡ e5bec8b4-2823-11eb-2cc8-45152db443b1
using Revise, NBodyLeapFrog, NiLang, BenchmarkTools, DoubleFloats

# ╔═╡ 66f740f0-2824-11eb-2722-bfb0d8f9a3ed
using Plots

# ╔═╡ fc600196-282d-11eb-2cbc-bbaefeb02467
# patch
for OP in [:PlusEq, :Mis]
    function (f::PlusEq{typeof(/)})(x::Double64, y::Double64, z::Double64)
        (x+y/z), y, z
    end
    function (f::MinusEq{typeof(/)})(x::Double64, y::Double64, z::Double64)
        (x-y/z), y, z
    end
end

# ╔═╡ cc384744-2823-11eb-2d52-41d7fd7cca29
planets = Bodies.chunit_day2year.(Bodies.set);

# ╔═╡ 0cf50074-2824-11eb-115e-f39e8ad092d2
function convertelems(::Type{T}, b::Body) where T
    Body(
        V3(T(b.r.x), T(b.r.y), T(b.r.z)),
        V3(T(b.v.x), T(b.v.y), T(b.v.z)),
        T(b.m)
    )
end

# ╔═╡ 2c6337d0-282e-11eb-1c64-259ac62c2cc0
false && let
	T = Double64
	nplanets = length(planets)
	ps = convertelems.(T1, planets)
	f(getfield.(ps, :r), getfield.(ps, :v), getfield.(ps, :m); n = n, dt = 0.01)
end

# ╔═╡ 3fc1ac08-2830-11eb-3f66-5bc4491ae8fe
function simulate(::Type{T1}, f, planets, k::Int) where T1
	nplanets = length(planets)
	n = 1<<k
	ps = convertelems.(T1, planets)
	r, v = f(getfield.(ps, :r), getfield.(ps, :v), getfield.(ps, :m); n = n, dt = 0.01)
	r[:,mod1(n+1,2)]
end

# ╔═╡ fd6c7b96-2823-11eb-1727-3b0cc6d50cf9
function errors(r, rref)
	sum(NBodyLeapFrog.distance.(r, rref))/length(r)
end

# ╔═╡ d37c5db4-2824-11eb-1e9c-9715ce7268eb
logns = 2:12

# ╔═╡ 88b3bb56-2824-11eb-2de1-e37dd1d015f4
let
	r_Double64 = simulate.(Double64, NBodyLeapFrog.i_leapfrog!, Ref(planets), logns)
	plt = plot([], [], yscale=:log10, xscale=:log10, label="")
	for T in [Float32, Float64]
		rt = simulate.(T, NBodyLeapFrog.i_leapfrog_reuse!, Ref(planets), logns)
		rt_clean = simulate.(T, NBodyLeapFrog.i_leapfrog!, Ref(planets), logns)
		errs1 = errors.(r_Double64, rt)
		errs2 = errors.(r_Double64, rt_clean)
		plot!(plt, 1 .<< logns, errs1; label="$T, reusing")
		plot!(plt, 1 .<< logns, errs2; label="$T")
		#errs1 = errors.(r_Double64, r_Double64_clean)
		#plot!(plt, 1 .<< logns, errs1; label="double")
	end
	plt
end

# ╔═╡ Cell order:
# ╠═e5bec8b4-2823-11eb-2cc8-45152db443b1
# ╠═fc600196-282d-11eb-2cbc-bbaefeb02467
# ╠═cc384744-2823-11eb-2d52-41d7fd7cca29
# ╠═0cf50074-2824-11eb-115e-f39e8ad092d2
# ╠═2c6337d0-282e-11eb-1c64-259ac62c2cc0
# ╠═3fc1ac08-2830-11eb-3f66-5bc4491ae8fe
# ╠═fd6c7b96-2823-11eb-1727-3b0cc6d50cf9
# ╠═d37c5db4-2824-11eb-1e9c-9715ce7268eb
# ╠═66f740f0-2824-11eb-2722-bfb0d8f9a3ed
# ╠═88b3bb56-2824-11eb-2de1-e37dd1d015f4
