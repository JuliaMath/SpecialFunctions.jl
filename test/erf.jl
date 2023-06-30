@testset "error functions" begin
    @testset "real argument" begin
        for T in (Float16, Float32, Float64)
            @test @inferred(erf(T(1))) isa T
            @test erf(T(1))   ≈ T(0.84270079294971486934) rtol=2*eps(T)

            @test @inferred(erfc(T(1))) isa T
            @test erfc(T(1))  ≈ T(0.15729920705028513066) rtol=2*eps(T)

            @test @inferred(erfcx(T(1))) isa T
            @test erfcx(T(1)) ≈ T(0.42758357615580700442) rtol=2*eps(T)

            @test @inferred(logerfc(T(1))) isa T
            @test logerfc(T(-100)) ≈ T(0.6931471805599453) rtol=2*eps(T)
            @test logerfc(T(1000)) ≈ T(-1.0000074801207219e6) rtol=2*eps(T)
            @test logerfc(T(10000)) ≈ T(log(erfc(BigFloat(10000, precision=100)))) rtol=2*eps(T)

            @test @inferred(logerfcx(T(1))) isa T
            @test logerfcx(T(1)) ≈ T(-0.849605509933248248576017509499) rtol=2eps(T)
            @test logerfcx(T(-1)) ≈ T(1.61123231767807049464268192445) rtol=2eps(T)
            @test logerfcx(T(-100)) ≈ T(10000.6931471805599453094172321) rtol=2eps(T)
            @test logerfcx(T(100)) ≈ T(-5.17758512266433257046678208395) rtol=2eps(T)
            @test logerfcx(T(-1000)) ≈ T(1.00000069314718055994530941723e6) rtol=2eps(T)
            @test logerfcx(T(1000)) ≈ T(-7.48012072190621214066734919080) rtol=2eps(T)

            @test @inferred(erfi(T(1))) isa T
            @test erfi(T(1)) ≈ T(1.6504257587975428760) rtol=2*eps(T)

            @test @inferred(erfinv(T(1))) isa T
            @test erfinv(T(0.84270079294971486934)) ≈ 1 rtol=2*eps(T)

            @test @inferred(erfcinv(T(1))) isa T
            @test erfcinv(T(0.15729920705028513066)) ≈ 1 rtol=2*eps(T)

            @test @inferred(dawson(T(1))) isa T
            @test dawson(T(1)) ≈ T(0.53807950691276841914) rtol=2*eps(T)

            @test @inferred(faddeeva(T(1))) isa Complex{T}
            @test faddeeva(T(1)) ≈ 0.36787944117144233402+0.60715770584139372446im rtol=2*eps(T)
        end

        @test logerfc(1000) ≈ -1.0000074801207219e6 rtol=2*eps(Float32)
        @test erfinv(Integer(0)) == 0 == erfinv(0//1)
        @test erfcinv(Integer(1)) == 0 == erfcinv(1//1)
        @test faddeeva(0) == faddeeva(0//1) == 1
    end

    @testset "complex arguments" begin
        @test erf(ComplexF16(1+2im)) ≈ -0.53664356577856503399-5.0491437034470346695im
        @test erf(ComplexF32(1+2im)) ≈ -0.53664356577856503399-5.0491437034470346695im
        @test erf(ComplexF64(1+2im)) ≈ -0.53664356577856503399-5.0491437034470346695im

        @test erfc(ComplexF16(1+2im)) ≈ 1.5366435657785650340+5.0491437034470346695im
        @test erfc(ComplexF32(1+2im)) ≈ 1.5366435657785650340+5.0491437034470346695im
        @test erfc(ComplexF64(1+2im)) ≈ 1.5366435657785650340+5.0491437034470346695im

        @test erfcx(ComplexF16(1+2im)) ≈ 0.14023958136627794370-0.22221344017989910261im
        @test erfcx(ComplexF32(1+2im)) ≈ 0.14023958136627794370-0.22221344017989910261im
        @test erfcx(ComplexF64(1+2im)) ≈ 0.14023958136627794370-0.22221344017989910261im

        @test erfi(ComplexF16(1+2im)) ≈ -0.011259006028815025076+1.0036063427256517509im
        @test erfi(ComplexF32(1+2im)) ≈ -0.011259006028815025076+1.0036063427256517509im
        @test erfi(ComplexF64(1+2im)) ≈ -0.011259006028815025076+1.0036063427256517509im

        @test_throws MethodError erfinv(Complex(1))

        @test_throws MethodError erfcinv(Complex(1))

        @test dawson(ComplexF16(1+2im)) ≈ -13.388927316482919244-11.828715103889593303im
        @test dawson(ComplexF32(1+2im)) ≈ -13.388927316482919244-11.828715103889593303im
        @test dawson(ComplexF64(1+2im)) ≈ -13.388927316482919244-11.828715103889593303im

        @test faddeeva(ComplexF16(1+2im)) ≈ 0.21849261527489066692+0.09299780939260188228im
        @test faddeeva(ComplexF32(1+2im)) ≈ 0.21849261527489066692+0.09299780939260188228im
        @test faddeeva(ComplexF64(1+2im)) ≈ 0.21849261527489066692+0.09299780939260188228im
    end

    @testset "BigFloat arguments" begin
        @test erf(BigFloat(1))  ≈ 0.84270079294971486934 rtol=2*eps()

        @test erfc(BigFloat(1)) ≈ 0.15729920705028513066 rtol=2*eps()

        @test erfcx(BigFloat(1))        ≈ 0.42758357615580700442    rtol=2*eps()
        @test erfcx(BigFloat(2e9))      ≈ 2.820947917738782e-10     rtol=2*eps()
        @test erfcx(BigFloat(1.8))      ≈ erfcx(1.8)                rtol=4*eps()
        @test erfcx(BigFloat(1.8e8))    ≈ erfcx(1.8e8)              rtol=4*eps()
        @test erfcx(BigFloat(1.8e88))   ≈ erfcx(1.8e88)             rtol=4*eps()
        @test isnan(erfcx(BigFloat(NaN)))

        @test logerfc(BigFloat(1000, precision=100)) ≈ -1.0000074801207219e6 rtol=2*eps(Float64)
        @test isnan(logerfc(BigFloat(NaN)))

        @test_throws MethodError erfi(big(1.0))

        @test_throws MethodError dawson(BigFloat(1))

        @test_throws MethodError faddeeva(BigFloat(1))

        for y in (big"1e-1000", big"1e-60", big"0.1", big"0.5", big"1.0", 1+big"1e-50", big"1.2", 2-big"1e-50")
            @test erfc(erfcinv(y)) ≈ y
        end
        for y in (big"1e-1000", big"1e-60", big"0.1", big"0.5", 1-big"1e-50")
            @test erf(erfinv(y)) ≈ y
            @test erf(erfinv(-y)) ≈ -y
        end
        @test erfcinv(big(0)) == -erfcinv(big(2)) == erfinv(big(1)) == -erfinv(big(-1)) == Inf
        for x in (1+big"1e-3", -1-big"1e-3")
            @test_throws DomainError erfinv(x)
            @test_throws DomainError erfcinv(1-x)
        end
    end

    @testset "Other float types" begin
        struct NotAFloat <: AbstractFloat end

        @test_throws MethodError erf(NotAFloat())
        @test_throws MethodError erfc(NotAFloat())
        @test_throws MethodError erfcx(NotAFloat())
        @test_throws MethodError erfi(NotAFloat())
        @test_throws MethodError erfinv(NotAFloat())
        @test_throws MethodError erfcinv(NotAFloat())
        @test_throws MethodError dawson(NotAFloat())
        @test_throws MethodError faddeeva(NotAFloat())
    end

    @testset "inverse" begin
        for elty in [Float32,Float64]
            for x in exp10.(range(-200, stop=-0.01, length=50))
                @test isapprox(erf(erfinv(x)), x, atol=1e-12*x)
                @test isapprox(erf(erfinv(-x)), -x, atol=1e-12*x)
                @test isapprox(erfc(erfcinv(2*x)), 2*x, atol=1e-12*x)
                if x > 1e-20
                    xf = Float32(x)
                    @test isapprox(erf(erfinv(xf)), xf, atol=1e-5*xf)
                    @test isapprox(erf(erfinv(-xf)), -xf, atol=1e-5*xf)
                    @test isapprox(erfc(erfcinv(2xf)), 2xf, atol=1e-5*xf)
                end
            end
            @test erfinv(one(elty))  == Inf
            @test erfinv(-one(elty)) == -Inf
            @test_throws DomainError erfinv(convert(elty,2.0))

            @test erfcinv(zero(elty)) == Inf
            @test_throws DomainError erfcinv(-one(elty))
        end

        @test erfinv(one(Int))  == erfinv(1.0)
        @test erfcinv(one(Int)) == erfcinv(1.0)
    end

    @testset "two args version" begin
        @test erf(10, 11) ≈ 2.08848758232167861905709161e-45
        @test erf(11, 10) ≈ -2.08848758232167861905709161e-45
        @test erf(-11, -10) ≈ 2.08848758232167861905709161e-45
        @test erf(-1, 1) ≈ 1.68540158589942973868244127017
        @test erf(1e-30, 2e-30) ≈ 1.12837916709551257389615890312e-30
        @test isnan(erf(NaN, 0))
        @test isnan(erf(0, NaN))
        @test isnan(erf(NaN, NaN))

        @test logerf(-5, 35) ≈ 0.693147180559176579520017808560
        @test logerf(30, 35) ≈ -903.974117110643878079600243618
        @test logerf(-35, -30) ≈ -903.974117110643878079600243618
        @test logerf(10, Inf) ≈ -102.879889024844888574804787140
        @test logerf(-Inf, Inf) ≈ 0.693147180559945309417232121458
        @test logerf(Inf, Inf) == -Inf
        @test logerf(-Inf, -Inf) == -Inf
        @test logerf(-Inf, Inf) ≈ log(2)
        @test logerf(-1e-6, 1e-6) ≈ -13.0015811397694169056785314498
        @test isnan(logerf(NaN, 0))
        @test isnan(logerf(0, NaN))
        @test isnan(logerf(NaN, NaN))
        @test logerf(-1e-30, 1e-30) ≈ -68.2636233716261799887769930733
        @test logerf(1e-30, 2e-30) ≈ -68.9567705521861252981942251947
        @test logerf(-2e-30, -1e-30) ≈ -68.9567705521861252981942251947
        @test_throws DomainError logerf(2, 1)
    end
end
