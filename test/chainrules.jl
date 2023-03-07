@testset "chainrules" begin
    Random.seed!(1)

    @testset "general: single input" begin
        for x in (1.0, -1.0, 0.0, 0.5, 10.0, -17.1, 1.5 + 0.7im)
            test_scalar(erf, x)
            test_scalar(erfc, x)
            test_scalar(erfcx, x)
            test_scalar(erfi, x)

            test_scalar(airyai, x)
            test_scalar(airyaiprime, x)
            test_scalar(airybi, x)
            test_scalar(airybiprime, x)

            test_scalar(dawson, x)

            if x isa Real
                test_scalar(logerfc, x)
                test_scalar(logerfcx, x)

                test_scalar(invdigamma, x)
            end

            if x isa Real && 0 < x < 1
                test_scalar(erfinv, x)
                test_scalar(erfcinv, x)
            end

            if x isa Real && x > 0 || x isa Complex
                test_scalar(gamma, x)
                test_scalar(digamma, x)
                test_scalar(trigamma, x)

                if x isa Real
                    test_scalar(airyaix, x)
                    test_scalar(airyaiprimex, x)
                end
            end

            if x isa Real && x < 1
                test_scalar(ellipk, x)
                test_scalar(ellipe, x)
            end
        end
    end

    @testset "Bessel functions" begin
        for x in (1.5, 2.5, 10.5, -0.6, -2.6, -3.3, 1.6 + 1.6im, 1.6 - 1.6im, -4.6 + 1.6im)
            test_scalar(besselj0, x)
            test_scalar(besselj1, x)

            isreal(x) && x < 0 && continue

            test_scalar(bessely0, x)
            test_scalar(bessely1, x)

            for nu in (-1.5, 2.2, 4.0)
                test_frule(besseli, nu, x)
                test_rrule(besseli, nu, x)
                test_frule(besselix, nu, x) # derivative is `NotImplemented`
                test_frule(besselix, nu ⊢ NoTangent(), x) # derivative is a number
                test_rrule(besselix, nu, x)

                test_frule(besselj, nu, x)
                test_rrule(besselj, nu, x)
                test_frule(besseljx, nu, x) # derivative is `NotImplemented`
                test_frule(besseljx, nu ⊢ NoTangent(), x) # derivative is a number
                test_rrule(besseljx, nu, x)

                test_frule(besselk, nu, x)
                test_rrule(besselk, nu, x)
                test_frule(besselkx, nu, x)
                test_rrule(besselkx, nu, x)

                test_frule(bessely, nu, x)
                test_rrule(bessely, nu, x)
                test_frule(besselyx, nu, x) # derivative is `NotImplemented`
                test_frule(besselyx, nu ⊢ NoTangent(), x) # derivative is a number
                test_rrule(besselyx, nu, x)

                test_frule(hankelh1, nu, x)
                test_rrule(hankelh1, nu, x)
                test_frule(hankelh1x, nu, x)
                test_rrule(hankelh1x, nu, x)

                test_frule(hankelh2, nu, x)
                test_rrule(hankelh2, nu, x)
                test_frule(hankelh2x, nu, x)
                test_rrule(hankelh2x, nu, x)
            end
        end
    end

    @testset "erf, beta, and logbeta" begin
        test_points = (1.5, 2.5, 10.5, 1.6 + 1.6im, 1.6 - 1.6im, 4.6 + 1.6im)
        for x in test_points, y in test_points
            test_frule(beta, x, y)
            test_rrule(beta, x, y)

            test_frule(logbeta, x, y)
            test_rrule(logbeta, x, y)

            if x isa Real && y isa Real
                test_frule(erf, x, y)
                test_rrule(erf, x, y)
            end
        end
    end

    @testset "log gamma and co" begin
        # It is important that we have negative numbers with both odd and even integer parts
        test_points = (1.5, 2.5, 10.5, -0.6, -2.6, -3.3, 1.6 + 1.6im, 1.6 - 1.6im, -4.6 + 1.6im)
        for x in test_points
            for m in (0, 1, 2, 3)
                test_frule(polygamma, m, x; rtol=1e-8)
                test_rrule(polygamma, m, x; rtol=1e-8)
            end

            isreal(x) && x < 0 && continue
            test_scalar(loggamma, x)
            for a in test_points
                test_frule(gamma, a, x; rtol=1e-8)
                test_rrule(gamma, a, x; rtol=1e-8)

                test_frule(loggamma, a, x)
                test_rrule(loggamma, a, x)
            end

            isreal(x) || continue
            test_frule(logabsgamma, x)
            test_rrule(logabsgamma, x; output_tangent=(randn(), randn()))
            for a in test_points
                isreal(a) && a > 0 || continue
                test_frule(gamma_inc, a, x, 0)
                test_rrule(gamma_inc, a, x, 0; output_tangent=(randn(), randn()))
            end
        end
    end

    @testset "exponential integrals" begin
        for x in (1.5, 2.5, 10.5, 1.6 + 1.6im, 1.6 - 1.6im, -4.6 + 1.6im)
            test_scalar(expint, x)
            test_scalar(expintx, x)

            for nu in (-1.5, 2.2, 4.0)
                test_frule(expint, nu, x)
                test_rrule(expint, nu, x)

                test_frule(expintx, nu, x)
                test_rrule(expintx, nu, x)
            end

            isreal(x) || continue
            test_scalar(expinti, x)
            test_scalar(sinint, x)
            test_scalar(cosint, x)
        end
    end

    # https://github.com/JuliaMath/SpecialFunctions.jl/issues/307
    @testset "promotions" begin
        # one argument
        for f in (erf, erfc, logerfc, erfcinv, erfcx, logerfcx, erfi, erfinv, sinint)
            _, ẏ = frule((NoTangent(), 1f0), f, 1f0)
            @test ẏ isa Float32
            _, back = rrule(f, 1f0)
            _, x̄ = back(1f0)
            @test x̄ isa Float32
        end

        # two arguments
        _, ẏ = frule((NoTangent(), 1f0, 1f0), erf, 1f0, 1f0)
        @test ẏ isa Float32
        _, back = rrule(erf, 1f0, 1f0)
        _, x̄ = back(1f0)
        @test x̄ isa Float32
    end
end
