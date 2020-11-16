using NBodyLeapFrog
using Test

@testset "NBodyLeapFrog.jl" begin
    # Write your tests here.
end

planets = Bodies.chunit_day2year.(Bodies.set)
nplanets = length(planets)
r, v, a = nOrbit(planets; n = 4000, dt = 0.01)

using Plots
function showres(vs)
    a = map(x->x.x, vs)
    b = map(x->x.y, vs)
    plt = plot([], [])
    for i in 1:nplanets
        plot!(plt,a[i,:],b[i,:])
    end
    plt
end

showres(r)

# python 5.6s
# before opt, 0.34s