### domain errors
using SpecialFunctions: inve # FIXME temporary until the fate of inve is decided
#using IrrationalConstants

@test_throws DomainError lambertw(-2.0, 0)
@test_throws DomainError lambertw(-2.0, -1)
@test_throws DomainError lambertw(-2.0, 1)
@test isnan(@inferred(lambertw(NaN)))

## math constant e
@test_throws DomainError lambertw(MathConstants.e, 1)
@test_throws DomainError lambertw(MathConstants.e, -1)
@test_throws DomainError lambertw(.3, 2)

## integer arguments return floating point types
@test @inferred(lambertw(0)) isa AbstractFloat
@test @inferred(lambertw(0)) == 0

### math constant, MathConstants.e e

# could return math const e, but this would break type stability
@test @inferred(lambertw(1)) isa AbstractFloat
@test @inferred(lambertw(1)) == float(LambertW.Omega)
@test @inferred(lambertw(big(1))) == big(LambertW.Omega)
@test @inferred(lambertw(MathConstants.e, 0)) == 1

## value at branch point where real branches meet
@test lambertw(-inve, 0) == lambertw(-inve, -1) == -1
@test typeof(lambertw(-inve, 0)) == typeof(lambertw(-inve, -1)) <: AbstractFloat

## convert irrationals to float

@test @inferred(lambertw(pi)) ≈ 1.0736581947961492
@test @inferred(lambertw(pi, 0)) ≈ 1.0736581947961492

### infinite args or return values

@test @inferred(lambertw(0, -1)) == @inferred(lambertw(0.0, -1)) == -Inf
@test @inferred(lambertw(Inf, 0)) == Inf
@test @inferred(lambertw(complex(Inf, 1), 0)) == complex(Inf, 1)
@test @inferred(lambertw(complex(Inf, 0), 1)) == complex(Inf, 2pi)
@test @inferred(lambertw(complex(-Inf, 0), 1)) == complex(Inf, 3pi)
@test @inferred(lambertw(complex(0.0, 0.0), -1)) == complex(-Inf, 0.0)

## default branch is k = 0
@test @inferred(lambertw(1.0)) == @inferred(lambertw(1.0, 0))

## BigInt args return BigFloats
@test @inferred(lambertw(BigInt(0))) isa BigFloat
@test @inferred(lambertw(BigInt(3))) isa BigFloat

## Any Integer type allowed for second argument
@test lambertw(-0.2, -1) == lambertw(-0.2, BigInt(-1))

## BigInt for second arg does not promote the type
@test typeof(lambertw(-0.2, -1)) == typeof(lambertw(-0.2, BigInt(-1)))

for (z, k, res) in [(0, 0 , 0), (complex(0, 0), 0 , 0),
                    (complex(0.0, 0), 0 , 0), (complex(1.0, 0), 0, 0.567143290409783873)]
    if Int != Int32
        @test lambertw(z, k) ≈ res
        @test lambertw(z) ≈ res
    else
        @test lambertw(z, k) ≈ res rtol=1e-14
        @test lambertw(z) ≈ res rtol=1e-14
    end
end

@testset "complex z=$z, k=$k" for (z, k) in
            ((complex(1, 1), 2), (complex(1, 1), 0), (complex(.6, .6), 0),
             (complex(.6, -.6), 0))
    w = lambertw(z, k)
    @test w*exp(w) ≈ z atol=1e-15
end

@test lambertw(complex(-3.0, -4.0), 0) ≈ Complex(1.075073066569255, -1.3251023817343588) atol=1e-14
@test lambertw(complex(-3.0, -4.0), 1) ≈ Complex(0.5887666813694675, 2.7118802109452247) atol=1e-14
@test lambertw(complex(.3, .3)) ≈ Complex(0.26763519642648767, 0.1837481231767825)

# bug fix
# The routine will start at -1/e + eps * im, rather than -1/e + 0im,
# otherwise root finding will fail
@test lambertw(-inve + 0im, -1) ≈ -1 atol=1e-7

# lambertw for BigFloat is more precise than Float64. Note
# that 70 digits in test is about 35 digits in W
@testset "lambertw() for BigFloat z=$z" for z in
            [BigFloat(1), BigFloat(2), complex(BigFloat(1), BigFloat(1))]
    W = lambertw(z)
    @test z ≈ W*exp(W) atol=BigFloat(10)^(-70)
end

@testset "LambertW.Omega" begin
    @test isapprox(LambertW.Ω * exp(LambertW.Ω), 1)
    @test LambertW.Omega === LambertW.Ω

    # lower than default precision
    setprecision(BigFloat, 196) do
        o = big(LambertW.Ω)
        @test precision(o) == 196
        @test isapprox(o * exp(o), 1, atol=eps(BigFloat))

        oalias = big(LambertW.Omega)
        @test o == oalias
    end

    # higher than default precision
    setprecision(BigFloat, 2048) do
        o = big(LambertW.Ω)
        @test precision(o) == 2048
        @test isapprox(o * exp(o), 1, atol=eps(BigFloat))

        oalias = big(LambertW.Omega)
        @test o == oalias
    end
end

###  expansion about branch point
@testset "lambertwbp()" begin
    # not a domain error, but not implemented
    @test_throws ArgumentError lambertwbp(1, 1)
    @test_throws ArgumentError lambertwbp(inve + 1e-5, 2)
    @test_throws DomainError lambertwbp(inve + 1e-5, 0)
    @test_throws DomainError lambertwbp(inve + 1e-5, -1)

    # Expansions about branch point converges almost to machine precision
    # except near the radius of convergence.
    # Complex args are not tested here.

    @testset "double-precision expansion near branch point using BigFloats" begin
        setprecision(2048) do
            z = BigFloat(10)^(-12)
            for _ in 1:300
                @test lambertwbp(Float64(z)) ≈ 1 + lambertw(z - big(inve)) atol=5e-16
                @test lambertwbp(Float64(z), -1) ≈ 1 + lambertw(z - big(inve), -1) atol=1e-15

                z *= 1.1
                if z > 0.23 break end
            end
        end
    end

    # test the expansion about branch point for k=-1,
    # by comparing to exact BigFloat calculation.
    @test @inferred(lambertwbp(1e-20, -1)) ≈ 1 + lambertw(-big(inve) + BigFloat(10)^(-20), -1) atol=1e-16
    @test @inferred(lambertwbp(Complex(.01, .01), -1)) ≈ Complex(-0.27550382080412062443536, -0.12778889284946406573511) atol=1e-16
end
