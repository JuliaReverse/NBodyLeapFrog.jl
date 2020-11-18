# NBodyLeapFrog

[![Build Status](https://travis-ci.com/JuliaReverse/NBodyLeapFrog.jl.svg?branch=master)](https://travis-ci.com/JuliaReverse/NBodyLeapFrog.jl)
[![Coverage](https://codecov.io/gh/JuliaReverse/NBodyLeapFrog.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/JuliaReverse/NBodyLeapFrog.jl)


This is translated from the following repository: https://github.com/jdlamstein/N-body-Gravity-Simulator


## To start
```julia-pkg
pkg> dev https://github.com/JuliaReverse/NBodyLeapFrog.jl.git

pkg> add Plots Pluto PlutoUI NiLang
```

Then open a Julia REPL and type
```julia
julia> using Pluto

julia> Pluto.run("~/.julia/dev/NBodyLeapFrog/notebooks/orbitals.jl")
```

Note the above path might vary in different systems.

You will see

![gradient](grad.png)

Enjoy!
