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

    @testset "beta_inc and beta_inc_inv" begin
        @testset "beta_inc and beta_inc_inv minimal (no-FD identities)" begin
            a = 1.2
            b = 2.3
            x = 0.4
            # Direct derivative checks without FD: ∂I/∂x equals beta pdf
            pdf = x^(a - 1) * (1 - x)^(b - 1) / beta(a, b)
            _, Δx = frule((NoTangent(), 0.0, 0.0, 1.0), beta_inc, a, b, x)
            @test isapprox(Δx[1], pdf; rtol=1e-12, atol=1e-12)

            # Symmetry check: ∂I/∂a(a,b,x) = -∂I/∂b(b,a,1-x)
            _, Δa = frule((NoTangent(), 1.0, 0.0, 0.0), beta_inc, a, b, x)
            _, Δb_sw = frule((NoTangent(), 0.0, 1.0, 0.0), beta_inc, b, a, 1 - x)
            @test isapprox(Δa[1], -Δb_sw[1]; rtol=1e-10, atol=1e-12)

            # Composition identity f(g(p)) = p: forward-mode differential equals 1 for dp, 0 for da,db
            p = first(beta_inc(a, b, x))
            x_inv, _ = beta_inc_inv(a, b, p)
            # Check primal composition
            p_roundtrip = first(beta_inc(a, b, x_inv))
            @test isapprox(p_roundtrip, p; rtol=1e-12, atol=1e-12)
            # Forward through g then f: dp
            _, Δx_inv_dp = frule((NoTangent(), 0.0, 0.0, 1.0), beta_inc_inv, a, b, p)
            _, Δp_from_dp = frule((NoTangent(), 0.0, 0.0, Δx_inv_dp[1]), beta_inc, a, b, x_inv)
            @test isapprox(Δp_from_dp[1], 1.0; rtol=1e-9, atol=1e-12)
            # Forward da
            _, Δx_inv_da = frule((NoTangent(), 1.0, 0.0, 0.0), beta_inc_inv, a, b, p)
            _, Δp_from_da = frule((NoTangent(), 1.0, 0.0, Δx_inv_da[1]), beta_inc, a, b, x_inv)
            @test isapprox(Δp_from_da[1], 0.0; rtol=1e-9, atol=1e-12)
            # Forward db
            _, Δx_inv_db = frule((NoTangent(), 0.0, 1.0, 0.0), beta_inc_inv, a, b, p)
            _, Δp_from_db = frule((NoTangent(), 0.0, 1.0, Δx_inv_db[1]), beta_inc, a, b, x_inv)
            @test isapprox(Δp_from_db[1], 0.0; rtol=1e-9, atol=1e-12)

            # Reverse-mode chain for composition: pullback through f then g
            # Pullback of f at (a,b,x_inv)
            _, pb_f = rrule(beta_inc, a, b, x_inv)
            _, āf, b̄f, x̄f = pb_f((1.0, 0.0))
            # Pullback of g at (a,b,p) with cotangent x̄f for x
            _, pb_g = rrule(beta_inc_inv, a, b, p)
            _, āg, b̄g, p̄g = pb_g((x̄f, 0.0))
            ā_total = āf + āg
            b̄_total = b̄f + b̄g
            p̄_total = p̄g
            @test isapprox(ā_total, 0.0; rtol=1e-10, atol=1e-12)
            @test isapprox(b̄_total, 0.0; rtol=1e-10, atol=1e-12)
            @test isapprox(p̄_total, 1.0; rtol=1e-9, atol=1e-12)
        end

        @testset "incomplete beta: basic test_frule/test_rrule" begin
            # Use a small, representative set of interior points (avoid endpoints for FD)
            test_points = (0.2, 0.5, 0.8)
            ab = (0.7, 2.5)

            # 3-argument beta_inc(a,b,x)
            for a in ab, b in ab, x in test_points
                0.0 < x < 1.0 || continue
                test_frule(beta_inc, a, b, x)
                test_rrule(beta_inc, a, b, x)
            end

            # Inverse beta: beta_inc_inv(a,b,p)
            for a in ab, b in ab, p in test_points
                0.0 < p < 1.0 || continue
                test_frule(beta_inc_inv, a, b, p)
                test_rrule(beta_inc_inv, a, b, p)
            end

            # Float32 promotion sanity (lightweight)
            a32 = 1.5f0; b32 = 2.25f0; x32 = 0.3f0
            # Finite-difference checks for Float32 are noisier; use looser tolerances
            test_frule(beta_inc, a32, b32, x32; rtol=5e-4, atol=1e-6)
            test_rrule(beta_inc, a32, b32, x32; rtol=5e-4, atol=1e-6)
            p32 = first(beta_inc(a32, b32, x32))
            test_frule(beta_inc_inv, a32, b32, p32; rtol=5e-4, atol=1e-6)
            test_rrule(beta_inc_inv, a32, b32, p32; rtol=5e-4, atol=1e-6)
        end

        @testset "4-arg beta_inc identities (y = 1 - x)" begin
            test_points = (0.2, 0.5, 0.8)
            ab = (0.7, 2.5)

            for a in ab, b in ab, x in test_points
                0.0 < x < 1.0 || continue
                y = 1 - x
                # Primal consistency: 4-arg matches 3-arg when y = 1 - x
                p3, q3 = beta_inc(a, b, x)
                p4, q4 = beta_inc(a, b, x, y)
                @test isapprox(p4, p3; rtol=1e-12, atol=1e-12)
                @test isapprox(q4, q3; rtol=1e-12, atol=1e-12)

                # Analytical pdf
                pdf = x^(a - 1) * (1 - x)^(b - 1) / beta(a, b)

                # Constrained x-variation: dx = 1, dy = -1 => dp = 2 * pdf, dq = -dp
                _, Δxy = frule((NoTangent(), 0.0, 0.0, 1.0, -1.0), beta_inc, a, b, x, y)
                @test isapprox(Δxy[1], 2 * pdf; rtol=1e-11, atol=1e-12)
                @test isapprox(Δxy[2], -Δxy[1]; rtol=1e-11, atol=1e-12)

                # Parameter derivatives should match 3-arg ones
                _, Δa3 = frule((NoTangent(), 1.0, 0.0, 0.0), beta_inc, a, b, x)
                _, Δb3 = frule((NoTangent(), 0.0, 1.0, 0.0), beta_inc, a, b, x)
                _, Δa4 = frule((NoTangent(), 1.0, 0.0, 0.0, 0.0), beta_inc, a, b, x, y)
                _, Δb4 = frule((NoTangent(), 0.0, 1.0, 0.0, 0.0), beta_inc, a, b, x, y)
                @test isapprox(Δa4[1], Δa3[1]; rtol=1e-11, atol=1e-12)
                @test isapprox(Δb4[1], Δb3[1]; rtol=1e-11, atol=1e-12)

                # Reverse-mode: compare pullbacks for 3-arg vs constrained 4-arg
                _, pb3 = rrule(beta_inc, a, b, x)
                _, ā3, b̄3, x̄3 = pb3((1.0, 0.0))
                _, pb4 = rrule(beta_inc, a, b, x, y)
                _, ā4, b̄4, x̄4, ȳ4 = pb4((1.0, 0.0))
                @test isapprox(ā4, ā3; rtol=1e-11, atol=1e-12)
                @test isapprox(b̄4, b̄3; rtol=1e-11, atol=1e-12)
                # Unconstrained pullbacks should satisfy x̄4 ≈ x̄3 and ȳ4 ≈ -x̄3
                @test isapprox(x̄4, x̄3; rtol=1e-11, atol=1e-12)
                @test isapprox(ȳ4, -x̄3; rtol=1e-11, atol=1e-12)
                # Effective pullback along the constraint y = 1 - x equals 2*x̄3
                x̄_eff = x̄4 - ȳ4
                @test isapprox(x̄_eff, 2 * x̄3; rtol=1e-11, atol=1e-12)
            end
        end

    end
end