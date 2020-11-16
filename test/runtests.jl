using NBodyLeapFrog
using NBodyLeapFrog: acceleration, sqdistance, distance
using Test
using NiLang, NiLang.AD

@testset "functions" begin
    ra, rb = V3(rand(3)...), V3(rand(3)...)
    y1 = ra + rb
    @instr ra += rb
    @test y1 ≈ ra

    ra, rb = V3(rand(3)...), V3(rand(3)...)
    y1 = ra + rb * 3.3
    @instr ra += rb * 3.3
    @test y1 ≈ ra

    ra, rb = V3(rand(3)...), V3(rand(3)...)
    y1 = ra + 3.3 * rb
    @instr ra += 3.3 * rb
    @test y1 ≈ ra

    ra, rb = V3(rand(3)...), V3(rand(3)...)
    y1 = distance(ra, rb)
    y2 = 0.0
    @instr y2 += sqdistance(ra, rb)
    @test y1^2 ≈ y2

    ra, rb = V3(rand(3)...), V3(rand(3)...)
    y1 = acceleration(ra, rb, 7.4, 2.0)
    y2 = zero(y1)
    @instr y2 += acceleration(ra, rb, 7.4, 2.0)
    @test y1 ≈ y2
end

@testset "NBodyLeapFrog.jl" begin
    planets = Bodies.chunit_day2year.(Bodies.set)
    nplanets = length(planets)
    r, v, a = nOrbit(planets; n = 55, dt = 0.01)
    r2, v2, a2 = fast_nOrbit(planets; n = 55, dt = 0.01)
    @test all(a->a[1]≈ a[2], zip(a[:,end], a2))
    @test all(r->r[1]≈r[2], zip(r[:,end], r2[:,end]))
end

@testset "NBodyLeapFrog.jl" begin
    planets = Bodies.chunit_day2year.(Bodies.set)
    nplanets = length(planets)
    n = 55
    r, v = nOrbit(planets; n = n, dt = 0.01)

    v2 = zeros(V3{Float64}, (nplanets,2))
    r2 = zeros(V3{Float64}, (nplanets,2))
    @instr i_leapfrog(r2, v2, planets; n = n, dt = 0.01)
    @test check_inv(i_leapfrog, (r2, v2, planets); n = n, dt = 0.01)
    @test all(r->r[1]≈r[2], zip(r[:,end], r2[:,end]))

    @i function loss(y!, r2, v2, planets; n, dt)
        i_leapfrog(r2, v2, planets; n = n, dt = dt)
        y! += sqdistance(r2[2,mod1(n+1,2)], r2[3,mod1(n+1,2)])
    end
    _, _, _, _, gp = Grad(loss)(Val(1), 0.0, copy(r2), copy(v2), planets; n = n, dt = 0.01)
    for i=1:nplanets
        # gradients on r
        for j=1:3
            ps = copy(planets)
            l0 = loss(0.0, copy(r2), copy(v2), copy(ps); n=n, dt=0.01)[1]

            z = zeros(3)
            z[j] += 1e-5
            ps[i] = Body(ps[i].r+V3(z...), ps[i].v, ps[i].m)
            l1 = loss(0.0, copy(r2), copy(v2), copy(ps); n=n, dt=0.01)[1]
            ng = (l1-l0)/1e-5
            @test isapprox(ng, gp[i].r[j].g; atol=1e-5, rtol=1e-2)
        end

        # gradients on v
        for j=1:3
            ps = copy(planets)
            l0 = loss(0.0, copy(r2), copy(v2), copy(ps); n=n, dt=0.01)[1]

            z = zeros(3)
            z[j] += 1e-5
            ps[i] = Body(ps[i].r, ps[i].v+V3(z...), ps[i].m)
            l1 = loss(0.0, copy(r2), copy(v2), copy(ps); n=n, dt=0.01)[1]
            ng = (l1-l0)/1e-5
            @test isapprox(ng, gp[i].v[j].g; atol=1e-5, rtol=1e-2)
        end

        # gradients on m
        ps = copy(planets)
        l0 = loss(0.0, copy(r2), copy(v2), copy(ps); n=n, dt=0.01)[1]

        ps[i] = Body(ps[i].r, ps[i].v, ps[i].m+1e-5)
        l1 = loss(0.0, copy(r2), copy(v2), copy(ps); n=n, dt=0.01)[1]
        ng = (l1-l0)/1e-5
        @test isapprox(ng, gp[i].m.g; atol=1e-5, rtol=1e-2)
    end
end