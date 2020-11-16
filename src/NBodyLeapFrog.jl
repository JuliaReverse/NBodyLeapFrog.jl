module NBodyLeapFrog

using LinearAlgebra, NiLang

export Bodies, nOrbit, fast_nOrbit
export Body, V3

include("V3.jl")
include("starData.jl")
include("Core.jl")
include("reversible.jl")

end
