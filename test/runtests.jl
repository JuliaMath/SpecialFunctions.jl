# This file contains code that was formerly a part of Julia. License is MIT: http://julialang.org/license

using SpecialFunctions
using Test
using Base.MathConstants: γ

using SpecialFunctions: AmosException, f64

# useful test functions for relative error, which differ from isapprox
# in that relerrc separately looks at the real and imaginary parts
relerr(z, x) = z == x ? 0.0 : abs(z - x) / abs(x)
relerrc(z, x) = max(relerr(real(z),real(x)), relerr(imag(z),imag(x)))
≅(a,b) = relerrc(a,b) ≤ 1e-13

tests = [
    "bessel",
    "beta_inc",
    "betanc",
    "ellip",
    "erf",
    "gamma_inc",
    "gamma",
    "sincosint",
    "other_tests"
]

const testdir = dirname(@__FILE__)


for t in tests
    tp = joinpath(testdir, "$(t).jl")
    @testset "$(t)" begin
      include(tp)
    end
end
