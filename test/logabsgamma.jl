# From openlibm/test/libm-test-ulps.h, openlibm/test/libm-test.c

# lgamma_test block
# for (T, lgamma) in ((Float64, _lgamma_r), (Float32, _lgammaf_r))
for T in (Float64, Float32)
    @testset "lgamma_test, $T" begin
        @test logabsgamma(T(Inf))[1] === T(Inf)
        @test logabsgamma(T(0))[1] === T(Inf)
        @test logabsgamma(T(NaN))[1] === T(NaN)

        @test logabsgamma(T(-3))[1] === T(Inf)
        @test logabsgamma(T(-Inf))[1] === T(Inf)

        # logabsgamma(1) == 0, lgamma (1) sets signgam to 1
        y, signgam = logabsgamma(T(1))
        @test y === T(0.0)
        @test signgam == 1

        # logabsgamma(3) == log(2), lgamma (3) sets signgam to 1
        y, signgam = logabsgamma(T(3))
        @test y === log(T(2.0))
        @test signgam == 1

        # logabsgamma(0.5) == log(sqrt(pi)), logabsgamma(0.5) sets signgam to 1
        y, signgam = logabsgamma(T(0.5))
        @test y === T(0.5log(π))
        @test signgam == 1

        # logabsgamma(-0.5) == log(2sqrt(pi)), logabsgamma(-0.5) sets signgam to -1
        y, signgam = logabsgamma(T(-0.5))
        @test y === T(0.5log(4π))
        @test signgam == -1

        # In the two "broken" tests, an exact match not possible, even
        # in Float64, thus, we check for as close a tolerance as
        # possible.

        # logabsgamma(0.7) == 0.26086724653166651439, logabsgamma(0.7) sets signgam to 1
        y, signgam = logabsgamma(T(0.7))
        # @test_broken y === 0.26086724653166651439
        if T === Float64
            @test y ≈ 0.26086724653166651439 atol=6e-17
        else
            @test y ≈ 0.26086724653166651439 atol=3e-8
        end
        @test signgam == 1

        # logabsgamma(1.2) == -0.853740900033158497197e-1, logabsgamma(1.2) sets signgam to 1
        y, signgam = logabsgamma(T(1.2))
        # @test_broken y === -0.853740900033158497197e-1
        if T === Float64
            @test y ≈ -0.853740900033158497197e-1 atol=2e-17
        else
            @test y ≈ -0.853740900033158497197e-1 atol=2e-8
        end
        @test signgam == 1
    end
end

# Comparison against the "gold standard" (in context of SpecialFunctions.jl)
using OpenLibm_jll

function openlibm_logabsgamma(x::Float64)
    signp = Ref{Int32}()
    y = ccall((:lgamma_r,libopenlibm),  Float64, (Float64, Ptr{Int32}), x, signp)
    return y, Int(signp[])
end
function openlibm_logabsgamma(x::Float32)
    signp = Ref{Int32}()
    y = ccall((:lgammaf_r,libopenlibm),  Float32, (Float32, Ptr{Int32}), x, signp)
    return y, Int(signp[])
end

meetstol(x, atol) = isapprox(openlibm_logabsgamma(x)[1], logabsgamma(x)[1], atol=atol)

@testset "logabsgamma validation against OpenLibm, Float64" begin
    @test all(x -> meetstol(x, 1e-13), -50:1e-4:50)
    @test all(x -> meetstol(x, 1e-12), 50:1e-4:500)
    @test all(x -> meetstol(x, 1e-11), 500:1e-3:1000)
    @test all(x -> meetstol(x, 1e-10), 1000:1e-1:10000)
    @test meetstol(-2.1^-71, 1e-15)
    @test meetstol(2.1^-71, 1e-15)
    @test meetstol(2.0^60, 1e-15)
end

@testset "logabsgamma validation against OpenLibm, Float32" begin
    @test all(x -> meetstol(x, 1f-5), -0.0f0:1f-4:25.0f0)
    @test all(x -> meetstol(x, 1f-4), -25.0f0:1f-4:0.0f0)
    @test all(x -> meetstol(x, 1f-4), 25.0f0:1f-4:100.0f0)
    @test all(x -> meetstol(x, 1f-3), 100.0f0:1f-4:1000.0f0)
    @test all(x -> meetstol(x, 1f-1), 1000.0f0:1f-1:10000.0f0)
    @test meetstol(-2.1f0^-71, 1f-6)
    @test meetstol(2.1f0^-71, 1f-6)
    @test meetstol(2.0f0^60, 1e-15)
end

#### A few critical values, but first, a few useful functions
ixword(x::Float64) = (reinterpret(UInt64, x) >>> 32 & Int32) & 0x7fffffff
ixword(x::Float32) = reinterpret(Int32, x) & 0x7fffffff

function ulp(x::Float32)
    z′ = logabsgamma(x)[1]
    z = logabsgamma(Float64(x))[1] # Test against our Float64 implementation
    isinf(z′) && isinf(oftype(x, z)) && return 0.0
    iszero(z′) && iszero(z) && return 0.0
    e = exponent(z′)
    abs(z′ - z) * 2.0^(precision(x) - 1 - e)
end
function ulp(x)
    z′ = logabsgamma(x)[1]
    z = logabsgamma(big(x))[1] # Dispatches to MPFR
    isinf(z′) && isinf(oftype(x, z)) && return big(0.0)
    iszero(z′) && iszero(z) && return big(0.0)
    e = exponent(z′)
    abs(z′ - z) * 2.0^(precision(x) - 1 - e)
end


# interval: 0 < x < 2
# First f64 which would not satisfy ixword(x) ≤ 0x40000000:
x = 2.000001907348633
# That is, all values between 2.0 and 2.000001907348633 would inappropriately fall
# into the trap of ≤ 0x40000000, thereby producing completely incorrect values.
# Hence, we must carefully test to ensure that this small region is computed
# appropriately.

@test ulp(x) < 0.15114421873672
@test ulp(2.0) == 0.0
r = 2.0:1e-9:prevfloat(x)
@test all(x -> ulp(x) < 1.0, r)

# first f32 greater than 2.0f0:
x = 2.0000002f0
@test ulp(x) < 0.197537131607533
# unlike f64, the representable distance between values is too small to matter
# -- i.e. prevfloat(2.0000002f0) == 2.0f0; we check behavior anyway
@test ulp(2.0f0) == 0.0


# interval: 2 < x < 8
# In fact, the first f64 which does not fall into the x < 8.0 branch is:
x = 8.000007629394531
@test ulp(prevfloat(x)) < 0.625524448102362
@test ulp(x) < 0.390118045607547
# This is overkill, but should reveal any toying around with this
# sensitive region.
r = 8.0:1e-9:prevfloat(x)
# Threshold of 1.5 can be met on everything except the combination of Julia v1.3, macOS-x64
# It would be preferable to reduce to 1.5 if/when the above combination is amenable.
@test all(x -> ulp(x) < 2.0, r)

# first f32 which does not fall into x < 8.0 branch is:
x = 8.000001f0
@test ulp(x) < 0.614982068538667
# But, unlike f64, the representable distance between values is too small.
# (i.e. prevfloat(8.000001f0) == 8.0f0)
# We still check appropriate behavior at 8.0f0
@test ulp(8.0f0) < 0.4006594736129046
