# This file contains code that was formerly a part of Julia. License is MIT: http://julialang.org/license

using SpecialFunctions
using Test
using Base.MathConstants: γ

using SpecialFunctions: AmosException, f64

# useful test functions for relative error, which differ from isapprox
# in that relerrc separately looks at the real and imaginary parts
relerr(z, x) = z == x ? 0.0 : abs(z - x) / abs(x)
relerr(z::T, x::T) where {T <: Complex} = max(relerr(real(z),real(x)), relerr(imag(z),imag(x)))
checktol(err::Float16) = err ≤ 5e-2
checktol(err::Float32) = err ≤ 1e-6
checktol(err::Float64) = err ≤ 1e-13
≅(a,b) = checktol(relerr(a,b))


tests = [
    "bessel",
    "beta_inc",
    "betanc",
    "ellip",
    "erf",
    "expint",
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
