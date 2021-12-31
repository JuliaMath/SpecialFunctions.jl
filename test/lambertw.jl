### domain errors

@test_throws DomainError lambertw(-2.0, 0)
@test_throws DomainError lambertw(-2.0, -1)
@test_throws DomainError lambertw(-2.0, 1)
@test isnan(@inferred(lambertw(NaN)))

## math constant e
@test_throws DomainError lambertw(MathConstants.e, 1)
@test_throws DomainError lambertw(MathConstants.e, -1)

## integer arguments return floating point types
@test @inferred(lambertw(0)) isa AbstractFloat
@test @inferred(lambertw(0)) == 0

### math constant, MathConstants.e e

# could return math const e, but this would break type stability
@test @inferred(lambertw(1)) isa AbstractFloat
@test @inferred(lambertw(MathConstants.e, 0)) == 1

## value at branch point where real branches meet
@test lambertw(-inv(MathConstants.e), 0) == lambertw(-inv(MathConstants.e), -1) == -1
@test typeof(lambertw(-inv(MathConstants.e), 0)) == typeof(lambertw(-inv(MathConstants.e), -1)) <: AbstractFloat

## convert irrationals to float

@test isapprox(@inferred(lambertw(pi)), 1.0736581947961492)
@test isapprox(@inferred(lambertw(pi, 0)), 1.0736581947961492)

### infinite args or return values

@test lambertw(0, -1) == lambertw(0.0, -1) == -Inf
@test lambertw(Inf, 0) == Inf
@test @inferred(lambertw(complex(Inf, 1), 0)) == complex(Inf, 1)
@test lambertw(complex(Inf, 0), 1) == complex(Inf, 2pi)
@test lambertw(complex(-Inf, 0), 1) == complex(Inf, 3pi)
@test @inferred(lambertw(complex(0.0, 0.0), -1)) == complex(-Inf, 0.0)

## default branch is  k = 0
@test lambertw(1.0) == lambertw(1.0, 0)

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
        @test isapprox(lambertw(z, k), res)
        @test isapprox(lambertw(z), res)
    else
        @test isapprox(lambertw(z, k), res; rtol = 1e-14)
        @test isapprox(lambertw(z), res; rtol = 1e-14)
    end
end

for (z, k) in ((complex(1, 1), 2), (complex(1, 1), 0), (complex(.6, .6), 0),
              (complex(.6, -.6), 0))
    let w
        @test (w = lambertw(z, k) ; true)
        @test abs(w*exp(w) - z) < 1e-15
    end
end

@test abs(lambertw(complex(-3.0, -4.0), 0) - Complex(1.075073066569255, -1.3251023817343588)) < 1e-14
@test abs(lambertw(complex(-3.0, -4.0), 1) - Complex(0.5887666813694675, 2.7118802109452247)) < 1e-14
@test (lambertw(complex(.3, .3), 0); true)

# bug fix
# The routine will start at -1/e + eps * im, rather than -1/e + 0im,
# otherwise root finding will fail
if Int != Int32
    @test abs(lambertw(-inv(MathConstants.e) + 0im, -1)) == 1
else
    @test abs(lambertw(-inv(MathConstants.e) + 0im, -1) + 1) < 1e-7
end
# lambertw for BigFloat is more precise than Float64. Note
# that 70 digits in test is about 35 digits in W
let W
    for z in [BigFloat(1), BigFloat(2), complex(BigFloat(1), BigFloat(1))]
        @test (W = lambertw(z); true)
        @test abs(z - W * exp(W)) < BigFloat(1)^(-70)
    end
end

###  ω constant

## get ω from recursion and compare to value from lambertw
let sp = precision(BigFloat)
    setprecision(512)
    @test lambertw(big(1)) == big(SpecialFunctions.omega)
    setprecision(sp)
end

@test lambertw(1) == float(SpecialFunctions.omega)
@test convert(Float16, SpecialFunctions.omega) == convert(Float16, 0.5674)
@test convert(Float32, SpecialFunctions.omega) == 0.56714326f0
@test lambertw(BigInt(1)) == big(SpecialFunctions.omega)

###  expansion about branch point

# not a domain error, but not implemented
@test_throws ArgumentError lambertwbp(1, 1)

@test_throws DomainError lambertw(.3, 2)

# Expansions about branch point converges almost to machine precision
# except near the radius of convergence.
# Complex args are not tested here.

if Int != Int32

# Test double-precision expansion near branch point using BigFloats
let sp = precision(BigFloat), z = BigFloat(1)/10^12, wo, xdiff
    setprecision(2048)
    for i in 1:300
        innerarg = z - inv(big(MathConstants.e))

        # branch k = 0
        @test (wo = lambertwbp(Float64(z)); xdiff = abs(-1 + wo - lambertw(innerarg)); true)
        if xdiff > 5e-16
            println(Float64(z), " ", Float64(xdiff))
        end
        @test xdiff < 5e-16

        #  branch k = -1
        @test (wo = lambertwbp(Float64(z), -1); xdiff = abs(-1 + wo - lambertw(innerarg, -1)); true)
        if xdiff > 5e-16
            println(Float64(z), " ", Float64(xdiff))
        end
        @test xdiff < 5e-16
        z  *= 1.1
        if z > 0.23 break end

    end
    setprecision(sp)
end

# test the expansion about branch point for k=-1,
# by comparing to exact BigFloat calculation.
@test @inferred(lambertwbp(1e-20, -1)) - 1 - lambertw(-inv(big(MathConstants.e)) + BigFloat(10)^BigFloat(-20), -1) < 1e-16

@test abs(@inferred(lambertwbp(Complex(.01, .01), -1)) - Complex(-0.2755038208041206, -0.1277888928494641)) < 1e-14

end  # if Int != Int32
