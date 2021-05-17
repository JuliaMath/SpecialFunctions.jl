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
        for x in (1.5, 2.5, 10.5, -0.6, -2.6, -3.3, 1.6 + 1.6im, 1.6 - 1.6im, -4.6 + 1.6im)
            for m in (0, 1, 2, 3)
                test_frule(polygamma, m, x)
                test_rrule(polygamma, m, x)
            end

            isreal(x) && x < 0 && continue
            test_scalar(loggamma, x)

            isreal(x) || continue
            test_frule(logabsgamma, x)
            test_rrule(logabsgamma, x; output_tangent=(randn(), randn()))
        end
    end
end
