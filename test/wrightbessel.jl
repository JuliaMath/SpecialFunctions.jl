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

@testset "`wrightbessel` systematic tests" begin
    # File contains precomputed values which can be reproduced using Nemo
    data = readdlm(joinpath(@__DIR__, "data/wrightbessel.txt"))
    for i in 1:size(data, 1)
        a, b, x, expected = data[i, :]
        #@test Float64(wrightbessel_nemo(a, b, x)) ≈ expected
        @test wrightbessel(a, b, x) ≈ expected rtol=1e-11
        @test exp(logwrightbessel(a, b, x)) ≈ expected rtol=1e-11
        @test logwrightbessel(a, b, x) ≈ log(expected) rtol=1e-11 atol=1e-15
    end
end

# Test cases that do not reach relative accuracy of 1e-11
@testset "`wrightbessel` low accuracy cases" begin
    low_accuracy = [
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

    for (a, b, x, expected, accuracy) in low_accuracy
        res = wrightbessel(a, b, x)
        #@test Float64(wrightbessel_nemo(a, b, x)) ≈ expected
        @test_broken res ≈ expected rtol=1e-11
        # Here we test for reduced accuracy or even NaN
        if isnan(accuracy)
            @test isnan(wrightbessel(a, b, x))
            @test isnan(logwrightbessel(a, b, x))
        else
            @test wrightbessel(a, b, x) ≈ expected rtol=accuracy
            @test logwrightbessel(a, b, x) ≈ log(expected) rtol=accuracy
        end
    end
end

@testset "`logwrightbessel` equals log of `wrightbessel`" begin
    for a in (0.0, 0.1, 0.5, 1.5, 5.0, 10.0)
        for b in (1.0, 2.0)
            for x in (1e-3, 1.0, 1.5, 5.0, 10.0)
                @test logwrightbessel(a, b, x) ≈ log(wrightbessel(a, b, x)) rtol=1e-8
            end
        end
    end
end

@testset "`logwrightbessel` extreme cases" begin
    log_test_cases = [
        (0.0, 0.0, 0.0, -Inf, 1e-11),
        (0.0, 0.0, 1.0, -Inf, 1e-11),
        (0.0, 1.0, 1.23, 1.23, 1e-11),
        (0.0, 1.0, 1e50, 1e50, 1e-11),
        (1e-5, 0.0, 700.0, 695.0421608273609, 1e-11),
        (1e-5, 0.0, 1e3, 995.40052566540066, 1e-11),
        (1e-5, 100.0, 1e3, 640.8197935670078, 1e-11),
        (1e-3, 0.0, 1e4, 9987.2229532297262, 1e-11),
        (1e-3, 0.0, 1e5, 99641.920687169507, 1e-11),
        (1e-3, 0.0, 1e6, 994118.55560054416, 1e-11),
        (1e-3, 10.0, 1e5, 99595.47710802537, 1e-11),
        (1e-3, 50.0, 1e5, 99401.240922855647, 1e-3),
        (1e-3, 100.0, 1e5, 99143.465191656527, NaN),
        (0.5, 0.0, 1e5, 4074.1112442197941, 1e-11),
        (0.5, 0.0, 1e7, 87724.552120038896, 1e-11),
        (0.5, 100.0, 1e5, 3350.3928746306163, NaN),
        (0.5, 100.0, 1e7, 86696.109975301719, 1e-11),
        (1.0, 0.0, 1e5, 634.06765787997266, 1e-11),
        (1.0, 0.0, 1e8, 20003.339639312035, 1e-11),
        (1.5, 0.0, 1e5, 197.01777556071194, 1e-11),
        (1.5, 0.0, 1e8, 3108.987414395706, 1e-11),
        (1.5, 100.0, 1e8, 2354.8915946283275, NaN),
        (5.0, 0.0, 1e5, 9.8980480013203547, 1e-11),
        (5.0, 0.0, 1e8, 33.642337258687465, 1e-11),
        (5.0, 0.0, 1e12, 157.53704288117429, 1e-11),
        (5.0, 100.0, 1e5, -359.13419630792148, 1e-11),
        (5.0, 100.0, 1e12, -337.07722086995229, 1e-4),
        (5.0, 100.0, 1e20, 2588.2471229986845, 2e-6),
        (100.0, 0.0, 1e5, -347.62127990460517, 1e-11),
        (100.0, 0.0, 1e20, -313.08250350969449, 1e-11),
        (100.0, 100.0, 1e5, -359.1342053695754, 1e-11),
        (100.0, 100.0, 1e20, -359.1342053695754, 1e-11),
    ]

    for (a, b, x, phi, accuracy) in log_test_cases
        if isnan(accuracy)
            @test isnan(logwrightbessel(a, b, x))
        else
            @test logwrightbessel(a, b, x) ≈ phi rtol=accuracy
        end
    end
end

@testset "rgamma_zero and exp_inf constants" begin
    @test SpecialFunctions.rgamma(SpecialFunctions.rgamma_zero) == 0
    @test SpecialFunctions.rgamma(prevfloat(SpecialFunctions.rgamma_zero)) > 0

    @test isinf(exp(SpecialFunctions.exp_inf))
    @test isfinite(exp(prevfloat(SpecialFunctions.exp_inf)))

    @test isnan(wrightbessel(SpecialFunctions.rgamma_zero+0.5, 1.0, 1.0))
    @test isnan(wrightbessel(1.0, SpecialFunctions.rgamma_zero+0.5, 1.0))

    @test isnan(logwrightbessel(SpecialFunctions.rgamma_zero+0.5, 1.0, 1.0))
    @test isnan(logwrightbessel(1.0, SpecialFunctions.rgamma_zero+0.5, 1.0))
end

@testset "$f errors" for f in (wrightbessel, logwrightbessel)
    @test_throws ArgumentError f(NaN, 1.0, 1.0)
    @test_throws ArgumentError f(1.0, NaN, 1.0)
    @test_throws ArgumentError f(1.0, 1.0, NaN)

    @test_throws ArgumentError f(Inf, 1.0, 1.0)
    @test_throws ArgumentError f(1.0, Inf, 1.0)
    @test f(1.0, 1.0, Inf) === Inf

    @test_throws ArgumentError f(-1.0, 1.0, 1.0)
    @test_throws ArgumentError f(1.0, -1.0, 1.0)
    @test_throws ArgumentError f(1.0, 1.0, -1.0)
end