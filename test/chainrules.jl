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
            @test isapprox(Δx[1], pdf; rtol=1e-12)

            # Symmetry check: ∂I/∂a(a,b,x) = -∂I/∂b(b,a,1-x)
            _, Δa = frule((NoTangent(), 1.0, 0.0, 0.0), beta_inc, a, b, x)
            _, Δb_sw = frule((NoTangent(), 0.0, 1.0, 0.0), beta_inc, b, a, 1 - x)
            @test isapprox(Δa[1], -Δb_sw[1]; rtol=1e-10)

            # Composition identity f(g(p)) = p: forward-mode differential equals 1 for dp, 0 for da,db
            p = first(beta_inc(a, b, x))
            x_inv, _ = beta_inc_inv(a, b, p)
            # Check primal composition
            p_roundtrip = first(beta_inc(a, b, x_inv))
            @test isapprox(p_roundtrip, p; rtol=1e-12)
            # Forward through g then f: dp
            _, Δx_inv_dp = frule((NoTangent(), 0.0, 0.0, 1.0), beta_inc_inv, a, b, p)
            _, Δp_from_dp = frule((NoTangent(), 0.0, 0.0, Δx_inv_dp[1]), beta_inc, a, b, x_inv)
            @test isapprox(Δp_from_dp[1], 1.0; rtol=1e-9)
            # Forward da
            _, Δx_inv_da = frule((NoTangent(), 1.0, 0.0, 0.0), beta_inc_inv, a, b, p)
            _, Δp_from_da = frule((NoTangent(), 1.0, 0.0, Δx_inv_da[1]), beta_inc, a, b, x_inv)
            @test isapprox(Δp_from_da[1], 0.0, rtol=1e-9, atol=1e-15)
            # Forward db
            _, Δx_inv_db = frule((NoTangent(), 0.0, 1.0, 0.0), beta_inc_inv, a, b, p)
            _, Δp_from_db = frule((NoTangent(), 0.0, 1.0, Δx_inv_db[1]), beta_inc, a, b, x_inv)
            @test isapprox(Δp_from_db[1], 0.0, rtol=1e-9, atol=1e-15)

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
            @test iszero(ā_total)
            @test iszero(b̄_total)
            @test isapprox(p̄_total, 1.0; rtol=1e-9)
        end

        @testset "incomplete beta: basic test_frule/test_rrule" begin
            # Use an expanded set of interior points (avoid endpoints for FD) to exercise many branches:
            # Rationale for x values:
            # - Include values around 0.1, 0.3, 0.5, 0.7, 0.9 to trigger different code paths.
            # - Include 0.14 and 0.28 to straddle the bx ≤ 0.7 power-series threshold for b ≈ 5 and 2.5.
            # - Include values near 0.5 (0.49, 0.51) to probe near-symmetry and tail swaps.
            # - Include additional midpoints to increase chance that x ≈ a/(a+b) for some (a,b), which makes λ ≈ 0
            #   in the large-parameter regime (key for choosing symmetric asymptotics when min(a,b) > 100).
            # - Add a few more around 0.6–0.8 to exercise continued fraction vs. asymptotics for large (a,b).
            test_points = (
                0.05, 0.08, 0.10, 0.12, 0.14, 0.18, 0.20, 0.22, 0.26,
                0.28, 0.30, 0.32, 0.35, 0.38, 0.40, 0.42, 0.45,
                0.49, 0.50, 0.51, 0.55, 0.58, 0.60, 0.62, 0.65,
                0.68, 0.70, 0.72, 0.76, 0.80, 0.85, 0.90
            )
            # Rationale for a,b values:
            # - <1: 0.4, 0.6 to stress small-parameter power series branches.
            # - Near 1: 0.9, 1.1 to test branch boundaries and continuity across a≈1, b≈1.
            # - Moderate: 2.5, 5.0 where multiple algorithm choices engage based on x and bx.
            # - Large (≥15, ≥40) to drive large-parameter regimes: 16.0, 45.0.
            # - Very large (≫100): 100.5, 150.0 to ensure symmetric vs asymmetric asymptotics are exercised when λ
            #   is small/large, and continued fractions are robust for large shapes.
            ab = (0.4, 0.6, 0.9, 1.1, 2.5, 5.0, 16.0, 45.0, 100.5, 150.0)

            # 3-argument beta_inc(a,b,x)
            for a in ab, b in ab, x in test_points
                0.0 < x < 1.0 || continue
                @testset "a=$a b=$b x=$x" begin
                    test_frule(beta_inc, a, b, x)
                    test_rrule(beta_inc, a, b, x)
                end
            end

            # Inverse beta: beta_inc_inv(a,b,p)
            for a in ab, b in ab, p in test_points
                0.0 < p < 1.0 || continue
                @testset "a=$a b=$b p=$p" begin
                    test_frule(beta_inc_inv, a, b, p)
                    test_rrule(beta_inc_inv, a, b, p)
                end
            end

            # Float32 promotion sanity (lightweight)
            a32 = 1.5f0; b32 = 2.25f0; x32 = 0.3f0
            # The Float32 primal limits the accuracy of the finite-difference reference;
            # fixed tangents keep these checks reproducible across RNG changes.
            test_frule(
                beta_inc, a32 ⊢ 0.75f0, b32 ⊢ -0.5f0, x32 ⊢ 0.25f0; rtol=5e-4,
            )
            test_rrule(
                beta_inc, a32 ⊢ 0.75f0, b32 ⊢ -0.5f0, x32 ⊢ 0.25f0;
                output_tangent=(0.75f0, -0.25f0), rtol=5e-4,
            )
            p32 = first(beta_inc(a32, b32, x32))
            # The inverse checks need a looser tolerance for the same reason; in
            # particular, their Float32 finite-difference rrule reference fails at 1e-4.
            test_frule(
                beta_inc_inv, a32 ⊢ 0.75f0, b32 ⊢ -0.5f0, p32 ⊢ 0.25f0; rtol=5e-4,
            )
            test_rrule(
                beta_inc_inv, a32 ⊢ 0.75f0, b32 ⊢ -0.5f0, p32 ⊢ 0.25f0;
                output_tangent=(0.75f0, -0.25f0), rtol=5e-4,
            )

            @testset "_beta_inc_grad sensitive cases" begin
                ext = Base.get_extension(SpecialFunctions, :SpecialFunctionsChainRulesCoreExt)
                @test ext !== nothing

                for (a, b, x) in (
                    (2.0, 1.0, 0.25),       # internal q == 1 without tail swap
                    (1.0, 2.0, 0.75),       # internal q == 1 after tail swap
                    (1e-6, 2.0, 0.01),      # tiny/imbalanced left-tail shape
                    (2.0, 1e-6, 0.99),      # tiny/imbalanced right-tail shape
                    (1000.0, 1000.0, 0.5),  # large central shape
                    (1e6, 1e6, 0.5),        # requires more than the historical 200 approximants
                    (1e8, 1e8, 0.5),        # requires thousands of approximants but converges below 10,000
                )
                    @test all(isfinite, ext._beta_inc_grad(a, b, x))
                end

                grad = @test_logs (:warn, r"_beta_inc_grad reached maxapp approximants before convergence") ext._beta_inc_grad(1.2, 2.3, 0.4; maxapp=2)
                @test all(isnan, grad)

                @test @inferred(ext._beta_inc_grad(1.2, 2.3, 0.4)) isa NTuple{3,Float64}
                beta_grad = ext._beta_inc_grad
                allocation_count(f) = @allocated f(1.2, 2.3, 0.4)
                beta_grad(1.2, 2.3, 0.4) # compile before allocation check
                @test allocation_count(beta_grad) == 0

                @test ext._beta_inc_grad(2.0, 3.0, 0.0) == (0.0, 0.0, 0.0)
                @test ext._beta_inc_grad(1.0, 3.0, 0.0) == (0.0, 0.0, 3.0)
                @test ext._beta_inc_grad(0.5, 3.0, 0.0) == (0.0, 0.0, Inf)
                @test ext._beta_inc_grad(3.0, 2.0, 1.0) == (0.0, 0.0, 0.0)
                @test ext._beta_inc_grad(3.0, 1.0, 1.0) == (0.0, 0.0, 3.0)
                @test ext._beta_inc_grad(3.0, 0.5, 1.0) == (0.0, 0.0, Inf)
                @test isequal(ext._beta_inc_grad(0.0, 2.0, 0.4), (NaN, NaN, 0.0))
                @test isequal(ext._beta_inc_grad(2.0, 0.0, 0.4), (NaN, NaN, 0.0))

                for args in ((NaN, 2.0, 0.5), (2.0, NaN, 0.5), (2.0, 3.0, NaN),
                             (Inf, 2.0, 0.5), (2.0, Inf, 0.5), (2.0, 3.0, Inf))
                    @test_logs @test all(isnan, ext._beta_inc_grad(args...))
                end

                # 600-bit references for the primal's cancellation-resistant beta density.
                for (a, b, x, expected) in (
                    (1e6, 1e8, 0.01, 6.530042978671665e-18),
                    (1e8, 1e8, 0.5, 11283.791656850386),
                )
                    @test ext._beta_inc_grad(a, b, x)[3] ≈ expected rtol=1e-12
                end

                # Ten minimum approximants prevent the eps floor from accepting a
                # visibly immature parameter derivative when its magnitude is tiny.
                tiny_grad = ext._beta_inc_grad(3e5, 1e6, 0.2345)
                mature_tiny_grad = ext._beta_inc_grad(3e5, 1e6, 0.2345; minapp=40)
                @test all(isapprox.(tiny_grad[1:2], mature_tiny_grad[1:2]; rtol=1e-12))

                @test ext._beta_inc_grad(Float16(1.2), Float16(2.3), Float16(0.4)) isa NTuple{3,Float16}
                _, delta16 = frule(
                    (NoTangent(), Float16(0), Float16(0), Float16(1)),
                    beta_inc, Float16(1.2), Float16(2.3), Float16(0.4),
                )
                @test delta16[1] isa Float16
                @test delta16[1] ≈ Float16(ext._beta_inc_grad(1.2, 2.3, 0.4)[3]) rtol=2e-3

                _, delta_int = frule((NoTangent(), 0.0, 0.0, 1.0), beta_inc_inv, 1, 2, 0.5)
                @test all(isfinite, delta_int)
                _, pullback_int = rrule(beta_inc_inv, 1, 2, 0.5)
                @test all(isfinite, pullback_int((1.0, 0.0))[2:4])
            end
        end

        @testset "4-arg beta_inc identities (y = 1 - x)" begin
            # All four arguments participate in promotion of the derivative calculation.
            mixed_args = (1.5f0, 2.25f0, 0.3f0, 1 - Float64(0.3f0))
            mixed_output = beta_inc(mixed_args...)
            mixed_derivatives = ChainRulesCore.derivatives_given_output(
                mixed_output, beta_inc, mixed_args...,
            )
            @test eltype(first(mixed_derivatives)) === Float64

            # Preserve the additional endpoint precision carried by the explicit y.
            endpoint_args = (2.0, 0.5, 1.0, 1e-16)
            endpoint_output = beta_inc(endpoint_args...)
            endpoint_derivatives = ChainRulesCore.derivatives_given_output(
                endpoint_output, beta_inc, endpoint_args...,
            )
            endpoint_dx = first(endpoint_derivatives)[3] - first(endpoint_derivatives)[4]
            expected_dx = SpecialFunctions.beta_integrand(
                endpoint_args..., -log(endpoint_args[3]) - log(endpoint_args[4]),
            )
            @test endpoint_dx ≈ expected_dx rtol=1e-12

            # Exercise more regimes while keeping y = 1 - x constraint.
            # Same rationale as above for x and (a,b) coverage.
            test_points = (
                0.05, 0.10, 0.12, 0.14, 0.20, 0.28, 0.35, 0.40, 0.49, 0.50, 0.51, 0.60, 0.65, 0.70, 0.72, 0.80, 0.90
            )
            ab = (0.4, 0.6, 0.9, 1.1, 2.5, 5.0, 16.0, 45.0, 100.5, 150.0)

            for a in ab, b in ab, x in test_points
                0.0 < x < 1.0 || continue
                y = 1 - x
                # Primal consistency: 4-arg matches 3-arg when y = 1 - x
                p3, q3 = beta_inc(a, b, x)
                p4, q4 = beta_inc(a, b, x, y)
                @test isapprox(p4, p3; rtol=1e-12)
                @test isapprox(q4, q3; rtol=1e-12)

                # Analytical pdf
                pdf = x^(a - 1) * (1 - x)^(b - 1) / beta(a, b)

                # Constrained x-variation: the symmetric split recovers one pdf.
                _, Δxy = frule((NoTangent(), 0.0, 0.0, 1.0, -1.0), beta_inc, a, b, x, y)
                @test isapprox(Δxy[1], pdf; rtol=1e-11)
                @test isapprox(Δxy[2], -Δxy[1]; rtol=1e-11)

                # Parameter derivatives should match 3-arg ones
                _, Δa3 = frule((NoTangent(), 1.0, 0.0, 0.0), beta_inc, a, b, x)
                _, Δb3 = frule((NoTangent(), 0.0, 1.0, 0.0), beta_inc, a, b, x)
                _, Δa4 = frule((NoTangent(), 1.0, 0.0, 0.0, 0.0), beta_inc, a, b, x, y)
                _, Δb4 = frule((NoTangent(), 0.0, 1.0, 0.0, 0.0), beta_inc, a, b, x, y)
                @test isapprox(Δa4[1], Δa3[1]; rtol=1e-11)
                @test isapprox(Δb4[1], Δb3[1]; rtol=1e-11)

                # Reverse-mode: compare pullbacks for 3-arg vs constrained 4-arg
                _, pb3 = rrule(beta_inc, a, b, x)
                _, ā3, b̄3, x̄3 = pb3((1.0, 0.0))
                _, pb4 = rrule(beta_inc, a, b, x, y)
                _, ā4, b̄4, x̄4, ȳ4 = pb4((1.0, 0.0))
                @test isapprox(ā4, ā3; rtol=1e-11)
                @test isapprox(b̄4, b̄3; rtol=1e-11)
                # The symmetric split assigns half of the constrained derivative to each input.
                @test isapprox(x̄4, x̄3 / 2; rtol=1e-11)
                @test isapprox(ȳ4, -x̄3 / 2; rtol=1e-11)
                # Effective pullback along the constraint y = 1 - x equals x̄3.
                x̄_eff = x̄4 - ȳ4
                @test isapprox(x̄_eff, x̄3; rtol=1e-11)
            end
        end

        @testset "4-arg beta_inc_inv identities (q = 1 - p)" begin
            a, b, p = 1.2, 2.3, 0.4
            x = first(beta_inc_inv(a, b, p))
            pdf = x^(a - 1) * (1 - x)^(b - 1) / beta(a, b)
            _, delta = frule((NoTangent(), 0.0, 0.0, 1.0, -1.0), beta_inc_inv, a, b, p, 1 - p)
            @test delta[1] ≈ inv(pdf) rtol=1e-12
            @test delta[2] ≈ -inv(pdf) rtol=1e-12
        end

    end
end
