# Type stability of the exported functions.
#
# The result type has to carry the precision that promoting the arguments gives: an
# `Integer` order does not by itself force `Float64`, and a `Float16` or `Float32`
# argument must not be widened. `@inferred(f(x)) isa T` checks both properties at
# once, since a method can be perfectly inferrable and still return the wrong type.
#
# The entries marked `broken` are the cases that get this wrong today.

@testset "error functions ($T)" for T in (Float16, Float32, Float64)
    @test @inferred(erf(T(1))) isa T
    @test @inferred(erfc(T(1))) isa T
    @test @inferred(erfcx(T(1))) isa T
    @test @inferred(erfi(T(1))) isa T
    @test @inferred(dawson(T(1))) isa T
    @test @inferred(erfinv(T(1)/2)) isa T
    @test @inferred(erfcinv(T(1)/2)) isa T
    @test @inferred(logerf(T(1)/2, T(1))) isa T
    @test @inferred(logerfc(T(1))) isa T
    @test @inferred(logerfcx(T(1))) isa T
end

@testset "error functions ($C)" for C in (ComplexF32, ComplexF64)
    @test @inferred(erf(C(1, 1))) isa C
    @test @inferred(erfc(C(1, 1))) isa C
    @test @inferred(erfcx(C(1, 1))) isa C
    @test @inferred(erfi(C(1, 1))) isa C
    @test @inferred(dawson(C(1, 1))) isa C
    @test @inferred(faddeeva(C(1, 1))) isa C
end

@testset "gamma family ($T)" for T in (Float16, Float32, Float64)
    @test @inferred(gamma(T(3)/2)) isa T
    @test @inferred(loggamma(T(3)/2)) isa T
    @test @inferred(logabsgamma(T(3)/2)) isa Tuple{T,Int}
    @test @inferred(digamma(T(3)/2)) isa T
    @test @inferred(trigamma(T(3)/2)) isa T
    @test @inferred(invdigamma(T(1))) isa T
    @test @inferred(polygamma(3, T(3)/2)) isa T
    @test @inferred(zeta(T(3))) isa T
    @test @inferred(zeta(T(3), T(2))) isa T
    @test @inferred(eta(T(3))) isa T
    @test @inferred(beta(T(3), T(2))) isa T
    @test @inferred(logbeta(T(3), T(2))) isa T
    @test @inferred(logabsbeta(T(3), T(2))) isa Tuple{T,Int}
end

@testset "gamma family ($C)" for C in (ComplexF32, ComplexF64)
    @test @inferred(gamma(C(3, 1))) isa C
    @test @inferred(loggamma(C(3, 1))) isa C
    @test @inferred(digamma(C(3, 1))) isa C
    @test @inferred(trigamma(C(3, 1))) isa C
    @test @inferred(polygamma(3, C(3, 1))) isa C
    @test @inferred(zeta(C(3, 1))) isa C
    @test @inferred(zeta(C(3, 1), C(2, 1))) isa C
    @test @inferred(eta(C(3, 1))) isa C
    @test @inferred(beta(C(3, 1), C(2, 1))) isa C
    @test @inferred(logbeta(C(3, 1), C(2, 1))) isa C
end

# `gamma(a, x)` and `loggamma(a, x)` are computed via `expint` for non-positive-integer
# `a`, which widens a `Float16`/`Float32` argument to `Float64`
@testset "incomplete gamma ($T)" for T in (Float16, Float32, Float64)
    @test @inferred(gamma(T(2), T(3))) isa T broken = T !== Float64
    @test @inferred(loggamma(T(2), T(3))) isa T broken = T !== Float64
end

@testset "incomplete gamma and beta ($T)" for T in (Float16, Float32, Float64)
    @test @inferred(gamma_inc(T(2), T(3), 0)) isa Tuple{T,T}
    @test @inferred(gamma_inc_inv(T(2), T(1)/2, T(1)/2)) isa T
    @test @inferred(beta_inc(T(2), T(3), T(1)/2)) isa Tuple{T,T}
    @test @inferred(beta_inc_inv(T(2), T(3), T(1)/2)) isa Tuple{T,T}
end

@testset "Airy functions ($T)" for T in (Float16, Float32, Float64)
    @test @inferred(airyai(T(1))) isa T
    @test @inferred(airyaiprime(T(1))) isa T
    @test @inferred(airybi(T(1))) isa T
    @test @inferred(airybiprime(T(1))) isa T
    @test @inferred(airyaix(T(1))) isa T
    @test @inferred(airyaiprimex(T(1))) isa T
    @test @inferred(airybix(T(1))) isa T
    @test @inferred(airybiprimex(T(1))) isa T
end

@testset "Airy functions ($C)" for C in (ComplexF32, ComplexF64)
    @test @inferred(airyai(C(1, 1))) isa C
    @test @inferred(airyaiprime(C(1, 1))) isa C
    @test @inferred(airybi(C(1, 1))) isa C
    @test @inferred(airybiprime(C(1, 1))) isa C
    @test @inferred(airyaix(C(1, 1))) isa C
    @test @inferred(airybix(C(1, 1))) isa C
end

@testset "elliptic and trigonometric integrals ($T)" for T in (Float16, Float32, Float64)
    @test @inferred(ellipk(T(1)/2)) isa T
    @test @inferred(ellipe(T(1)/2)) isa T
    @test @inferred(sinint(T(1))) isa T
    @test @inferred(cosint(T(1))) isa T
end

@testset "exponential integrals ($T)" for T in (Float16, Float32, Float64)
    @test @inferred(expint(T(1))) isa T
    @test @inferred(expintx(T(1))) isa T
    @test @inferred(expinti(T(1))) isa T
end

@testset "exponential integrals ($C)" for C in (ComplexF32, ComplexF64)
    @test @inferred(expint(C(1, 1))) isa C
    @test @inferred(expintx(C(1, 1))) isa C
end

# `expint(ν, z)` widens `Float16`/`Float32` arguments to `Float64` on several paths
@testset "exponential integrals with an order ($T)" for T in (Float16, Float32, Float64)
    @test @inferred(expint(2, T(3))) isa T broken = T !== Float64
    @test @inferred(expint(T(3)/2, T(3))) isa T broken = T !== Float64
    @test @inferred(expintx(2, T(3))) isa T broken = T !== Float64
end

@testset "exponential integrals with an order ($C)" for C in (ComplexF32, ComplexF64)
    @test @inferred(expint(2, C(1, 1))) isa C broken = C !== ComplexF64
    @test @inferred(expintx(2, C(1, 1))) isa C broken = C !== ComplexF64
end

@testset "Bessel functions of order 0 and 1 ($T)" for T in (Float16, Float32, Float64)
    @test @inferred(besselj0(T(1))) isa T
    @test @inferred(besselj1(T(1))) isa T
    @test @inferred(bessely0(T(1))) isa T
    @test @inferred(bessely1(T(1))) isa T
    @test @inferred(jinc(T(1))) isa T
end

@testset "Bessel functions of general order ($T)" for T in (Float16, Float32, Float64)
    @test @inferred(besselj(T(3)/2, T(1))) isa T
    @test @inferred(bessely(T(3)/2, T(1))) isa T
    @test @inferred(besseli(2, T(1))) isa T
    @test @inferred(besselk(2, T(1))) isa T
    @test @inferred(besseljx(2, T(1))) isa T
    @test @inferred(besselyx(2, T(1))) isa T
    @test @inferred(besselix(2, T(1))) isa T
    @test @inferred(besselkx(2, T(1))) isa T
    @test @inferred(sphericalbesselj(2, T(1))) isa T
    @test @inferred(sphericalbessely(2, T(1))) isa T
    @test @inferred(hankelh1(2, T(1))) isa Complex{T}
    @test @inferred(hankelh2(2, T(1))) isa Complex{T}
end

@testset "Bessel functions of integer order ($T)" for T in (Float16, Float32, Float64)
    @test @inferred(besselj(2, T(1))) isa T
    @test @inferred(bessely(2, T(1))) isa T
end

@testset "Bessel functions of general order ($C)" for C in (ComplexF32, ComplexF64)
    @test @inferred(besselj0(C(1, 1))) isa C
    @test @inferred(besselj1(C(1, 1))) isa C
    @test @inferred(bessely0(C(1, 1))) isa C
    @test @inferred(bessely1(C(1, 1))) isa C
    @test @inferred(jinc(C(1, 1))) isa C
    @test @inferred(besselj(2, C(1, 1))) isa C
    @test @inferred(bessely(2, C(1, 1))) isa C
    @test @inferred(besseli(2, C(1, 1))) isa C
    @test @inferred(besselk(2, C(1, 1))) isa C
    @test @inferred(besseljx(2, C(1, 1))) isa C
    @test @inferred(besselyx(2, C(1, 1))) isa C
    @test @inferred(besselix(2, C(1, 1))) isa C
    @test @inferred(besselkx(2, C(1, 1))) isa C
    @test @inferred(hankelh1(2, C(1, 1))) isa C
    @test @inferred(hankelh2(2, C(1, 1))) isa C
    @test @inferred(hankelh1x(2, C(1, 1))) isa C
    @test @inferred(hankelh2x(2, C(1, 1))) isa C
    @test @inferred(sphericalbesselj(2, C(1, 1))) isa C
    @test @inferred(sphericalbessely(2, C(1, 1))) isa C
end

# `Integer` and `Rational` arguments promote to `Float64`
@testset "integer and rational arguments" begin
    @test @inferred(erf(1)) isa Float64
    @test @inferred(erf(1//2)) isa Float64
    @test @inferred(erfc(1)) isa Float64
    @test @inferred(erfcx(1)) isa Float64
    @test @inferred(erfinv(1//2)) isa Float64
    @test @inferred(erfcinv(1//2)) isa Float64
    @test @inferred(gamma(3//2)) isa Float64
    @test @inferred(logabsgamma(3//2)) isa Tuple{Float64,Int}
    @test @inferred(digamma(2)) isa Float64
    @test @inferred(trigamma(2)) isa Float64
    @test @inferred(invdigamma(1)) isa Float64
    @test @inferred(polygamma(3, 2)) isa Float64
    @test @inferred(zeta(3)) isa Float64
    @test @inferred(eta(3)) isa Float64
    @test @inferred(beta(3, 2)) isa Float64
    @test @inferred(logbeta(3, 2)) isa Float64
    @test @inferred(logfactorial(5)) isa Float64
    @test @inferred(gamma(2, 3)) isa Float64
    @test @inferred(loggamma(2, 3)) isa Float64
    @test @inferred(gamma(2, 3.0)) isa Float64
    @test @inferred(gamma(2.0, 3)) isa Float64
    @test @inferred(gamma(2, 3//1)) isa Float64
    @test @inferred(gamma_inc(2, 3, 0)) isa Tuple{Float64,Float64}
    @test @inferred(beta_inc(2, 3, 1//2)) isa Tuple{Float64,Float64}
    @test @inferred(besselj0(1)) isa Float64
    @test @inferred(besselj(2, 1)) isa Float64
    @test @inferred(besseli(2, 1)) isa Float64
    @test @inferred(sphericalbesselj(2, 1)) isa Float64
    @test @inferred(hankelh1(2, 1)) isa ComplexF64
    @test @inferred(jinc(1)) isa Float64
    @test @inferred(airyai(1)) isa Float64
    @test @inferred(ellipk(1//2)) isa Float64
    @test @inferred(sinint(1)) isa Float64
    @test @inferred(cosint(1)) isa Float64
    @test @inferred(expint(1)) isa Float64
    @test @inferred(expinti(2)) isa Float64
    @test @inferred(expint(2, 3)) isa Float64
    @test @inferred(expintx(2, 3)) isa Float64
    @test @inferred(expint(2, 3//2)) isa Float64
    # the second element is documented as the sign, so it should be an `Int` like
    # `logabsgamma`'s, but it comes back as a `Float64`
    @test @inferred(logabsbinomial(5, 2)) isa Tuple{Float64,Int} broken = true
end

@testset "BigFloat arguments" begin
    @test @inferred(erf(big(1))) isa BigFloat
    @test @inferred(gamma(big(3)/2)) isa BigFloat
    @test @inferred(zeta(big(3))) isa BigFloat
    @test @inferred(gamma(big(2), big(3))) isa BigFloat
    @test @inferred(loggamma(big(2), big(3))) isa BigFloat
    @test @inferred(besselj0(big(1))) isa BigFloat
    @test @inferred(besselj(2, big(1))) isa BigFloat
    @test @inferred(expinti(big(2))) isa BigFloat
    @test @inferred(expint(big(2), big(3)/2)) isa BigFloat
    @test @inferred(expintx(big(2), big(3)/2)) isa BigFloat
    # `expint(x::BigFloat)` is inferred as `Union{Float64,BigFloat}`
    @test @inferred(expint(big(1))) isa BigFloat broken = true
    # `airyai` lacks the `f(x::Real) = f(float(x))` promotion that `besselj0` has, so a
    # `BigInt` never reaches its `BigFloat` method
    @test @inferred(airyai(big(1))) isa BigFloat broken = true
end
