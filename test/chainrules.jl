@testset "chainrules" begin
    Random.seed!(1)

    @testset "general" begin
        for x in (1.0, -1.0, 0.0, 0.5, 10.0, -17.1, 1.5 + 0.7im)
            test_scalar(erf, x)
            test_scalar(erfc, x)
            test_scalar(erfi, x)

            test_scalar(airyai, x)
            test_scalar(airyaiprime, x)
            test_scalar(airybi, x)
            test_scalar(airybiprime, x)

            test_scalar(besselj0, x)
            test_scalar(besselj1, x)

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
                test_scalar(bessely0, x)
                test_scalar(bessely1, x)
                test_scalar(gamma, x)
                test_scalar(digamma, x)
                test_scalar(trigamma, x)
            end
        end

        @testset "beta and logbeta" begin
            test_points = (1.5, 2.5, 10.5, 1.6 + 1.6im, 1.6 - 1.6im, 4.6 + 1.6im)
            for _x in test_points, _y in test_points
                # ensure all complex if any complex for FiniteDifferences
                x, y = promote(_x, _y)
                T = typeof(x)

                Δx, x̄ = randn(T, 2)
                Δy, ȳ = randn(T, 2)
                Δz = randn(T)

                frule_test(beta, (x, Δx), (y, Δy))
                rrule_test(beta, Δz, (x, x̄), (y, ȳ))

                frule_test(logbeta, (x, Δx), (y, Δy))
                rrule_test(logbeta, Δz, (x, x̄), (y, ȳ))
            end
        end

        @testset "log gamma and co" begin
            # It is important that we have negative numbers with both odd and even integer parts
            for x in (1.5, 2.5, 10.5, -0.6, -2.6, -3.3, 1.6 + 1.6im, 1.6 - 1.6im, -4.6 + 1.6im)
                isreal(x) && x < 0 && continue
                test_scalar(loggamma, x)

                isreal(x) || continue

                Δx, x̄ = randn(2)
                Δz = (randn(), randn())

                frule_test(logabsgamma, (x, Δx))
                rrule_test(logabsgamma, Δz, (x, x̄))
            end
        end
    end
end
