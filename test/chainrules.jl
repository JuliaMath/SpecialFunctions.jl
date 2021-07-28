@testset "chainrules" begin
    Random.seed!(1)

    @testset "general: single input" begin
        for x in (1.0, -1.0, 0.0, 0.5, 10.0, -17.1, 1.5 + 0.7im)
            test_scalar(erf, x)
            test_scalar(erfc, x)
            test_scalar(erfi, x)

            test_scalar(airyai, x)
            test_scalar(airyaiprime, x)
            test_scalar(airybi, x)
            test_scalar(airybiprime, x)

            test_scalar(erfcx, x)
            test_scalar(dawson, x)

            if x isa Real
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

                test_frule(besselj, nu, x)
                test_rrule(besselj, nu, x)

                test_frule(besselk, nu, x)
                test_rrule(besselk, nu, x)

                test_frule(bessely, nu, x)
                test_rrule(bessely, nu, x)

                # use complex numbers in `rrule` for FiniteDifferences
                test_frule(hankelh1, nu, x)
                test_rrule(hankelh1, nu, complex(x))

                # use complex numbers in `rrule` for FiniteDifferences
                test_frule(hankelh2, nu, x)
                test_rrule(hankelh2, nu, complex(x))
            end
        end
    end

    @testset "incomplete beta" begin
        test_points = (1e-20, 0.5, 0.8, 0.9, 0.99, 1.5, 1.7, 10.5, 100.5)
        for a in test_points, b in test_points, x in 0.0:0.01:1.0
            test_frule(beta_inc, a, b, x)
            test_rrule(beta_inc, a, b, x)
        end
    end

    @testset "beta and logbeta" begin
        test_points = (1.5, 2.5, 10.5, 1.6 + 1.6im, 1.6 - 1.6im, 4.6 + 1.6im)
        for _x in test_points, _y in test_points
            # ensure all complex if any complex for FiniteDifferences
            x, y = promote(_x, _y)
            test_frule(beta, x, y)
            test_rrule(beta, x, y)

            test_frule(logbeta, x, y)
            test_rrule(logbeta, x, y)
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
                # ensure all complex if any complex for FiniteDifferences
                _a, _x = promote(a, x)
                test_frule(gamma, _a, _x; rtol=1e-8)
                test_rrule(gamma, _a, _x; rtol=1e-8)

                test_frule(loggamma, _a, _x)
                test_rrule(loggamma, _a, _x)
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
                # ensure all complex if any complex for FiniteDifferences
                _x, _nu = promote(x, nu)

                test_frule(expint, _nu, _x)
                test_rrule(expint, _nu, _x)

                test_frule(expintx, _nu, _x)
                test_rrule(expintx, _nu, _x)
            end

            isreal(x) || continue
            test_scalar(expinti, x)
            test_scalar(sinint, x)
            test_scalar(cosint, x)
        end
    end
end
