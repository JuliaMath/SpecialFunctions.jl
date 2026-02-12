using DelimitedFiles

# using Nemo
# """
# Compute Wright's generalized Bessel function as a series using Nemo
# for arbitrary precision.
# """
# function wrightbessel_nemo(a, b, x; n=1000)
#     ctx = RealField()
#     x′ = ctx(x)
#
#     series_term(k) =  x′^k / factorial(k, ctx) * rgamma(ctx(k * a + b))
#
#     res = zero(x)
#     k = 0
#     for _ in 1:n
#         res += series_term(k)
#         k += 1
#     end
#
#     return res
# end

@testset "`wrightbessel` with x=0 (rgamma)" begin
    for a in (0.0, 1e-6, 0.1, 0.5, 1.0, 10.0)
        for b in (0.0, 1e-6, 0.1, 0.5, 1.0, 10.0)
            @test wrightbessel(a, b, 0.0) == SpecialFunctions.rgamma(b)
        end
    end
end

# Test relation of `wrightbessel` and modified Bessel function of the first kind.
# besseli(z) = (1/2*z)^v * Φ(1, v+1; 1/4*z^2).
# See https://dlmf.nist.gov/10.46.E2
@testset "`wrightbessel` with a=1 (besseli)" begin
    for b in (0.0, 1e-6, 0.1, 0.5, 1.0, 10.0)
        for x in (0.0, 1e-6, 0.1, 0.5, 1.0)
            if x != 0
                v = b - 1
                wb = wrightbessel(1.0, v + 1, x^2 / 4)
                @test (x / 2)^v * wb ≈ besseli(v, x)
            end
        end
    end
end

# Test functional relation of `wrightbessel`.
#
# Φ(a, b-1, z) = a*z*Φ(a, b+a, z) + (b-1)*Φ(a, b, z)
#
# Note that d/dx Φ(a, b, x) = Φ(a, b-1, x)
# See Eq. (22) of
# B. Stankovic, On the Function of E. M. Wright,
# Publ. de l'Institut Mathématique, Beograd,
# Nouvelle Série. 10 (1970), 113-124.
@testset "`wrightbessel` functional relation" begin
    for a in (0.0, 1e-6, 0.1, 0.5, 1.0, 10.0)
        for b in (1.0, 1.0 + 1e-3, 2.0, 5.0, 10.0)
            for x in (0.0, 1e-6, 0.1, 0.5, 1.0, 5.0, 10.0, 100.0)
                lhs = wrightbessel(a, b - 1, x)
                rhs = a * x * wrightbessel(a, b + a, x) + (b - 1) * wrightbessel(a, b, x)
                @test lhs ≈ rhs
            end
        end
    end
end

# grid of rows (a, b, x, expected, accuracy) that do not reach 1e-11 accuracy
const low_accuracy = [
    (0.1, 100.0, 709.7827128933841, 8.026353022981087e+34, 2e-8),
    (0.5, 10.0, 709.7827128933841, 2.680788404494657e+48, 9e-8),
    (0.5, 10.0, 1000.0, 2.005901980702872e+64, 1e-8),
    (0.5, 100.0, 1000.0, 3.4112367580445246e-117, 6e-8),
    (1.0, 20.0, 100000.0, 1.7717158630699857e+225, 3e-11),
    (1.0, 100.0, 100000.0, 1.0269334596230763e+22, NaN),
    (1.0000000000000222, 20.0, 100000.0, 1.7717158630001672e+225, 3e-11),
    (1.0000000000000222, 100.0, 100000.0, 1.0269334595866202e+22, NaN),
    (1.5, 0.0, 500.0, 15648961196.432373, 3e-11),
    (1.5, 2.220446049250313e-14, 500.0, 15648961196.431465, 3e-11),
    (1.5, 1e-10, 500.0, 15648961192.344728, 3e-11),
    (1.5, 1e-05, 500.0, 15648552437.334162, 3e-11),
    (1.5, 0.1, 500.0, 12049870581.10317, 2e-11),
    (1.5, 20.0, 100000.0, 7.81930438331405e+43, 3e-9),
    (1.5, 100.0, 100000.0, 9.653370857459075e-130, NaN)
]

# Test cases of test_data that do not reach relative accuracy of 1e-11
@testset "`wrightbessel` data grid failures" begin
    for (a, b, x, expected, accuracy) in low_accuracy
        res = wrightbessel(a, b, x)
        #@test Float64(wrightbessel_nemo(a, b, x)) ≈ expected
        @test_broken res ≈ expected rtol=1e-11
        # Here we test for reduced accuracy or even NaN
        if isnan(accuracy)
            @test isnan(wrightbessel(a, b, x))
        else
            @test wrightbessel(a, b, x) ≈ expected rtol=accuracy
        end
    end
end

# File contains precomputed values which can be reproduced using Nemo
data = readdlm(joinpath(@__DIR__, "data/wrightbessel.txt"))
for i in 1:size(data, 1)
    a, b, x, expected = data[i, :]
    #@test Float64(wrightbessel_nemo(a, b, x)) ≈ expected
    @test wrightbessel(a, b, x) ≈ expected rtol=1e-11
end

@testset "rgamma_zero and exp_inf constants" begin
    @test SpecialFunctions.rgamma(SpecialFunctions.rgamma_zero) == 0
    @test SpecialFunctions.rgamma(prevfloat(SpecialFunctions.rgamma_zero)) > 0

    @test isinf(exp(SpecialFunctions.exp_inf))
    @test isfinite(exp(prevfloat(SpecialFunctions.exp_inf)))

    @test isnan(wrightbessel(SpecialFunctions.rgamma_zero+0.5, 1.0, 1.0))
    @test isnan(wrightbessel(1.0, SpecialFunctions.rgamma_zero+0.5, 1.0))
end

@testset "wrightbessel errors" begin
    @test_throws ArgumentError wrightbessel(NaN, 1.0, 1.0)
    @test_throws ArgumentError wrightbessel(1.0, NaN, 1.0)
    @test_throws ArgumentError wrightbessel(1.0, 1.0, NaN)

    @test_throws ArgumentError wrightbessel(Inf, 1.0, 1.0)
    @test_throws ArgumentError wrightbessel(1.0, Inf, 1.0)
    @test wrightbessel(1.0, 1.0, Inf) === Inf

    @test_throws ArgumentError wrightbessel(-1.0, 1.0, 1.0)
    @test_throws ArgumentError wrightbessel(1.0, -1.0, 1.0)
    @test_throws ArgumentError wrightbessel(1.0, 1.0, -1.0)
end