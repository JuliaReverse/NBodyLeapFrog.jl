### A Pluto.jl notebook ###
# v0.12.10

using Markdown
using InteractiveUtils

# ╔═╡ 83a2c2cc-296b-11eb-38bf-fd97aa3c6fda
using Revise, PlutoUI, NBodyLeapFrog, Plots

# ╔═╡ 814c93d0-296c-11eb-333d-216ce303ce34
plotly();

# ╔═╡ 8e419bc2-296b-11eb-0e0f-afdab78b1392
function showres(vs)
    a = map(x->x.x, vs)
    b = map(x->x.y, vs)
    plt = plot([], [])
    for i in 1:size(a, 1)
        plot!(plt,a[i,:],b[i,:]; aspect_ratio=:equal)
    end
    plt
end

# ╔═╡ de67315c-296b-11eb-0ae2-93757c2e2331
planets = Bodies.chunit_day2year.(Bodies.set);

# ╔═╡ bbe14adc-296b-11eb-1cb2-c149ac6e6321
r1, v1 = leapfrog(FR, planets; keep_history=true, n=4000, dt=0.01);

# ╔═╡ 6b217206-296c-11eb-0e56-a742db172508
r2, v2 = leapfrog(PEFRL, planets; keep_history=true, n=4000, dt=0.01);

# ╔═╡ ba521b86-296b-11eb-2bd4-ebf29f78b87a
showres(vcat(r1, r2))

# ╔═╡ Cell order:
# ╠═83a2c2cc-296b-11eb-38bf-fd97aa3c6fda
# ╠═814c93d0-296c-11eb-333d-216ce303ce34
# ╠═8e419bc2-296b-11eb-0e0f-afdab78b1392
# ╠═de67315c-296b-11eb-0ae2-93757c2e2331
# ╠═bbe14adc-296b-11eb-1cb2-c149ac6e6321
# ╠═6b217206-296c-11eb-0e56-a742db172508
# ╠═ba521b86-296b-11eb-2bd4-ebf29f78b87a
