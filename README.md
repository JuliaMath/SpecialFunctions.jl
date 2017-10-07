# SpecialFunctions.jl

Special mathematical functions in Julia, including Bessel, Hankel, Airy, error, Dawson, sine and cosine integrals,
eta, zeta, digamma, inverse digamma, trigamma, and polygamma functions.
Most of these functions were formerly part of Base.

Note: On Julia 0.7, this package downloads and/or builds
[openspecfun](https://github.com/JuliaLang/openspecfun), which is no longer built as part
of Julia.
Binaries are available for macOS, Windows, and Linux (glibc >= 2.6).
To force compilation of the library from source, set an environment variable called
`JULIA_SPECIALFUNCTIONS_BUILD_SOURCE` equal to `true` before running `Pkg.build`.

[![Travis](https://travis-ci.org/JuliaMath/SpecialFunctions.jl.svg?branch=master)](https://travis-ci.org/JuliaMath/SpecialFunctions.jl)
[![AppVeyor](https://ci.appveyor.com/api/projects/status/ccfgkm2cjcggu158/branch/master?svg=true)](https://ci.appveyor.com/project/ararslan/specialfunctions-jl/branch/master)
[![Coverage Status](https://coveralls.io/repos/github/JuliaMath/SpecialFunctions.jl/badge.svg?branch=master)](https://coveralls.io/github/JuliaMath/SpecialFunctions.jl?branch=master)

Latest release:
[![SpecialFunctions](http://pkg.julialang.org/badges/SpecialFunctions_0.6.svg)](http://pkg.julialang.org/?pkg=SpecialFunctions)

Documentation:
[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://JuliaMath.github.io/SpecialFunctions.jl/stable)
[![](https://img.shields.io/badge/docs-latest-blue.svg)](https://JuliaMath.github.io/SpecialFunctions.jl/latest)
