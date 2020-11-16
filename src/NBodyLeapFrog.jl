module NBodyLeapFrog

using LinearAlgebra, NiLang

export Bodies, leapfrog, fast_leapfrog
export Body, V3

include("V3.jl")
include("starData.jl")
include("Core.jl")
include("reversible.jl")

end
