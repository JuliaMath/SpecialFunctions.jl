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

    @testset "incomplete gamma functions: large arguments" begin
        for (a, x) in ((150.0, 150.0), (1e6, 1e6))
            (p, q), (∂p, ∂q) = frule((NoTangent(), NoTangent(), 1.0, NoTangent()), gamma_inc, a, x, 0)
            @test (p, q) == gamma_inc(a, x, 0)
            @test ∂p == -∂q == gamma_inc_dx(a, x)
            Ω, pullback = rrule(gamma_inc, a, x, 0)
            @test Ω == (p, q)
            @test pullback(Tangent{typeof(Ω)}(2.0, 1.0))[3] == ∂p
        end
        test_frule(gamma_inc, 150.0, 150.0, 0)
        test_rrule(gamma_inc, 150.0, 150.0, 0)

        # Γ(n, x) = (n - 1)! e^(-x) Σ_{k < n} x^k / k!
        function ∂loggamma(n, x)
            x = big(x)
            return -(x^(n - 1) / factorial(big(n - 1))) / sum(x^k / factorial(big(k)) for k in 0:(n - 1))
        end
        for (n, x) in ((2, 1e10), (30, 1e6), (5, 0.1))
            Ω, ∂Ω = frule((NoTangent(), NoTangent(), 1.0), loggamma, n, x)
            @test Ω == loggamma(n, x)
            @test ∂Ω ≈ ∂loggamma(n, x) rtol = 4 * eps()
            Ω, pullback = rrule(loggamma, n, x)
            @test Ω == loggamma(n, x)
            @test pullback(1.0)[3] == ∂Ω
        end
    end

    @testset "gamma_inc_dx" begin
        # includes the poles of `digamma(a)` at `a = 0, -2`
        for (a, x) in ((2.5, 1.0), (0.3, 2.0), (-2.5, 1.5), (-2.0, 1.5), (0.0, 1.5), (30.0, 25.0), (150.0, 10.0))
            test_frule(gamma_inc_dx, a, x)
            test_rrule(gamma_inc_dx, a, x)
        end

        function partials(a, x)
            Ω, pullback = rrule(gamma_inc_dx, a, x)
            _, ∂a, ∂x = pullback(one(Ω))
            return Ω, ∂a, ∂x
        end
        # f = x^(a - 1) e^(-x) / Γ(a), ∂f/∂a = f * (log(x) - ψ(a)), ∂f/∂x = f * (a - 1 - x) / x
        function partials_ref(a, x)
            setprecision(BigFloat, 512) do
                a, x = big(a), big(x)
                # 1 / Γ(a) = sin(πa) Γ(1 - a) / π
                e = exp((a - 1) * log(x) - x)
                f = a > 0 ? e / gamma(a) : e * sinpi(a) * gamma(1 - a) / big(π)
                ∂a = isinteger(a) && a <= 0 ? cospi(a) * e * gamma(1 - a) : f * (log(x) - digamma(a))
                return f, ∂a, f * (a - 1 - x) / x
            end
        end
        for T in (Float32, Float64), (a, x) in (
                (2.5, 1.0), (12.0, 11.5), (150.0, 150.0), (150.0, 149.5), (1e6, 1e6), (1e6, 999_000.0), (1e10, 1e10),
                (0.0, 1.5), (-2.0, 1.5),
            )
            Ω, ∂a, ∂x = partials(T(a), T(x))
            _, ∂a_ref, ∂x_ref = partials_ref(T(a), T(x))
            @test Ω === gamma_inc_dx(T(a), T(x))
            @test ∂a isa T
            @test ∂a ≈ T(∂a_ref) rtol = 8 * eps(T)
            @test ∂x isa T
            @test ∂x ≈ T(∂x_ref) rtol = 8 * eps(T)
        end
        for (a, x) in ((2.5, 1.0), (0.0, 1.5), (-2.0, 1.5), (1e6, 1e6))
            Ω, ∂a, ∂x = partials(big(a), big(x))
            _, ∂a_ref, ∂x_ref = partials_ref(a, x)
            @test ∂a ≈ ∂a_ref rtol = 8 * eps(BigFloat)
            @test ∂x ≈ ∂x_ref rtol = 8 * eps(BigFloat)
        end
        @test all(p -> p isa Float16, partials(Float16(2.5), Float16(1)))

        # `f` underflows
        @test partials(3.0, 1e-160)[3] ≈ 1e-160 rtol = 4 * eps()
        for T in (Float32, Float64)
            @test partials(T(2), T(0)) === (T(0), T(0), T(1))
            @test partials(T(2.5), T(Inf)) === (T(0), T(0), T(0))
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
