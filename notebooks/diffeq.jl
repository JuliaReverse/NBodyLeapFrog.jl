### A Pluto.jl notebook ###
# v0.12.11

using Markdown
using InteractiveUtils

# ╔═╡ 92ae949c-2a00-11eb-0f11-6570ad25144b
using Revise, DifferentialEquations

# ╔═╡ 897e68d2-2a3f-11eb-32e4-0b56dec432e8
using NBodyLeapFrog

# ╔═╡ 4aaf7264-2a45-11eb-1a8e-7976e5bbc634
using StaticArrays

# ╔═╡ 4ba29d18-2a63-11eb-16cf-d13492677a59
using PlutoUI

# ╔═╡ 4f659e96-2a45-11eb-03ab-01dbf5e82eb4
tosvec(x::V3) = SVector(x.x, x.y, x.z)

# ╔═╡ 7e0214da-2a45-11eb-26a5-9bb0c6b1518a
NBodyLeapFrog.acceleration(x::SVector, y::SVector, m, G) = NBodyLeapFrog.acceleration(V3(x...), V3(y...), m, G)

# ╔═╡ e0db6922-2a0d-11eb-040f-8ba9cd6d6f5b
function solar_problem(r0, v0; dt, n, params=nothing)
	function dp(v::AbstractVector{T}, x, params, t) where T
		@show x
		res = zeros(T, length(x))
		for i=1:length(res)
			for j=1:length(res)
				if i!= j
					res[i] += NBodyLeapFrog.acceleration(x[i], x[j], params.m[j], params.G)
				end
			end
		end
		return res
	end
	function dq(v, x, params, t)
		@show v
		v
	end
	ODEProblem(DynamicalODEFunction{false}(dp,dq), ArrayPartition(v0,r0), (0.0, n*dt), params)
end

# ╔═╡ e7a72072-2a49-11eb-3bec-3766113cf185
ArrayPartition

# ╔═╡ bf0d3bc4-2a3f-11eb-2c09-c7e1e709235d
v, r, m = let
	n = 100
    planets = Bodies.chunit_day2year.(Bodies.set)
    getfield.(planets, :v),
    getfield.(planets, :r),
    getfield.(planets, :m)
end

# ╔═╡ a4179630-2a48-11eb-0738-53f92fca7ef6
prob = solar_problem(r, v, dt=0.01, n=4000, params=(m=m, G=NBodyLeapFrog.G_year_AU))

# ╔═╡ d6dc5db2-2a48-11eb-2e5a-2f3574592e7e
DiffEqBase.recursive_unitless_bottom_eltype(::Type{V3{T}}) where T = T#V3{T}

# ╔═╡ 0a5d0150-2a49-11eb-3a04-cbe528382566
function DiffEqBase.recursive_unitless_eltype(::Type{V3{T}}) where T
	V3{T}
end

# ╔═╡ b47297c8-2a48-11eb-0f80-b352b8100a64
DiffEqBase.recursive_unitless_bottom_eltype(r |> typeof)

# ╔═╡ 6b99eb4a-2a49-11eb-2a9a-b5cf52a26ded
begin
	Base.one(::Type{V3{T}}) where T = one(T)
	Base.one(::V3{T}) where T = one(T)
	Base.oneunit(x::V3) = one(x)
	Base.:(/)(x::V3, y::Real) = V3(x.x/y, x.y/y, x.z/y)
	Base.broadcastable(x::V3) = Ref(x)
	Base.isnan(x::V3) = isnan(x.x) || isnan(x.y) || isnan(x.z)
end

# ╔═╡ 2a4f0bce-2a45-11eb-19e8-41512b3e66be
 #solve(prob)

# ╔═╡ 462eb43c-2a01-11eb-1726-593e84125e48
alg = McAte2()

# ╔═╡ 2f7fb8e6-2a01-11eb-229c-210763e47bca
integrator = DiffEqBase.__init(prob, alg; dt=0.01, tstops=40.0)

# ╔═╡ 32f10746-2a01-11eb-01ec-6b90be753181
 solve!(integrator)

# ╔═╡ 372f56c6-2a01-11eb-0826-432eefb6198e
integrator.sol

# ╔═╡ Cell order:
# ╠═92ae949c-2a00-11eb-0f11-6570ad25144b
# ╠═897e68d2-2a3f-11eb-32e4-0b56dec432e8
# ╠═4aaf7264-2a45-11eb-1a8e-7976e5bbc634
# ╠═4f659e96-2a45-11eb-03ab-01dbf5e82eb4
# ╠═7e0214da-2a45-11eb-26a5-9bb0c6b1518a
# ╠═e0db6922-2a0d-11eb-040f-8ba9cd6d6f5b
# ╠═e7a72072-2a49-11eb-3bec-3766113cf185
# ╠═bf0d3bc4-2a3f-11eb-2c09-c7e1e709235d
# ╠═a4179630-2a48-11eb-0738-53f92fca7ef6
# ╠═d6dc5db2-2a48-11eb-2e5a-2f3574592e7e
# ╠═0a5d0150-2a49-11eb-3a04-cbe528382566
# ╠═b47297c8-2a48-11eb-0f80-b352b8100a64
# ╠═6b99eb4a-2a49-11eb-2a9a-b5cf52a26ded
# ╠═2a4f0bce-2a45-11eb-19e8-41512b3e66be
# ╠═462eb43c-2a01-11eb-1726-593e84125e48
# ╠═2f7fb8e6-2a01-11eb-229c-210763e47bca
# ╠═32f10746-2a01-11eb-01ec-6b90be753181
# ╠═372f56c6-2a01-11eb-0826-432eefb6198e
# ╠═4ba29d18-2a63-11eb-16cf-d13492677a59
