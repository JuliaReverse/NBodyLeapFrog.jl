using NBodyLeapFrog
using Test

@testset "NBodyLeapFrog.jl" begin
    planets = Bodies.chunit_day2year.(Bodies.set)
    nplanets = length(planets)
    r, v, a = nOrbit(planets; n = 55, dt = 0.01)
    r2, v2, a2 = fast_nOrbit(planets; n = 55, dt = 0.01)
    @test all(a->Bodies.distance(a[1],a[2]) < 1e-8, zip(a[:,end], a2))
end