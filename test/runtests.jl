using NBodyLeapFrog
using Test

@testset "NBodyLeapFrog.jl" begin
    # Write your tests here.
end

#=
a,b,c,va,vb,vc,aa,ab,ac,energy, ChangeInE,totalH, ChangeInH,hRatio,hTotalMag, rBary,vBary = nOrbit(x,y,z,x,y,z,vx,vy,vz,vx,vy,vz,G,m,m,N,planets,dt);

plt = plot([], [])
for i in 1:planets
    plot!(plt,a[i,:],b[i,:])
end
=#