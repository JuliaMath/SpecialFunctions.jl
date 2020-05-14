@testset "error functions" begin
    @testset "real argument" begin
        @test erf(Float16(1))   ≈ 0.84270079294971486934 rtol=2*eps(Float16)
        @test erf(Float32(1))   ≈ 0.84270079294971486934 rtol=2*eps(Float32)
        @test erf(Float64(1))   ≈ 0.84270079294971486934 rtol=2*eps(Float64)

        @test erfc(Float16(1))  ≈ 0.15729920705028513066 rtol=2*eps(Float16)
        @test erfc(Float32(1))  ≈ 0.15729920705028513066 rtol=2*eps(Float32)
        @test erfc(Float64(1))  ≈ 0.15729920705028513066 rtol=2*eps(Float64)

        @test_throws MethodError erfcx(Float16(1))
        @test erfcx(Float32(1)) ≈ 0.42758357615580700442 rtol=2*eps(Float32)
        @test erfcx(Float64(1)) ≈ 0.42758357615580700442 rtol=2*eps(Float64)

        @test_throws MethodError logerfc(Float16(1))
        @test_throws MethodError logerfc(Float16(-1))
        @test logerfc(Float32(-100)) ≈ 0.6931471805599453 rtol=2*eps(Float32)
        @test logerfc(Float64(-100)) ≈ 0.6931471805599453 rtol=2*eps(Float64)
        @test logerfc(Float32(1000)) ≈ -1.0000074801207219e6 rtol=2*eps(Float32)
        @test logerfc(Float64(1000)) ≈ -1.0000074801207219e6 rtol=2*eps(Float64)
        @test logerfc(1000) ≈ -1.0000074801207219e6 rtol=2*eps(Float32)
        @test logerfc(Float32(10000)) ≈ log(erfc(BigFloat(10000, 100))) rtol=2*eps(Float32)
        @test logerfc(Float64(10000)) ≈ log(erfc(BigFloat(10000, 100))) rtol=2*eps(Float64)

        @test_throws MethodError logerfcx(Float16(1))
        @test_throws MethodError logerfcx(Float16(-1))
        @test iszero(logerfcx(0))
        @test logerfcx(Float32(1)) ≈ -0.849605509933248248576017509499 rtol=2eps(Float32)
        @test logerfcx(Float64(1)) ≈ -0.849605509933248248576017509499 rtol=2eps(Float32)
        @test logerfcx(Float32(-1)) ≈ 1.61123231767807049464268192445 rtol=2eps(Float32)
        @test logerfcx(Float64(-1)) ≈ 1.61123231767807049464268192445 rtol=2eps(Float32)
        @test logerfcx(Float32(-100)) ≈ 10000.6931471805599453094172321 rtol=2eps(Float32)
        @test logerfcx(Float64(-100)) ≈ 10000.6931471805599453094172321 rtol=2eps(Float64)
        @test logerfcx(Float32(100)) ≈ -5.17758512266433257046678208395 rtol=2eps(Float32)
        @test logerfcx(Float64(100)) ≈ -5.17758512266433257046678208395 rtol=2eps(Float64)
        @test logerfcx(Float32(-1000)) ≈ 1.00000069314718055994530941723e6 rtol=2eps(Float32)
        @test logerfcx(Float64(-1000)) ≈ 1.00000069314718055994530941723e6 rtol=2eps(Float64)
        @test logerfcx(Float32(1000)) ≈ -7.48012072190621214066734919080 rtol=2eps(Float32)
        @test logerfcx(Float64(1000)) ≈ -7.48012072190621214066734919080 rtol=2eps(Float64)

        @test_throws MethodError erfi(Float16(1))
        @test erfi(Float32(1)) ≈ 1.6504257587975428760 rtol=2*eps(Float32)
        @test erfi(Float64(1)) ≈ 1.6504257587975428760 rtol=2*eps(Float64)

        @test erfinv(Integer(0)) == 0
        @test_throws MethodError erfinv(Float16(1))
        @test erfinv(Float32(0.84270079294971486934)) ≈ 1 rtol=2*eps(Float32)
        @test erfinv(Float64(0.84270079294971486934)) ≈ 1 rtol=2*eps(Float64)

        @test erfcinv(Integer(1)) == 0
        @test_throws MethodError erfcinv(Float16(1))
        @test erfcinv(Float32(0.15729920705028513066)) ≈ 1 rtol=2*eps(Float32)
        @test erfcinv(Float64(0.15729920705028513066)) ≈ 1 rtol=2*eps(Float64)

        @test_throws MethodError dawson(Float16(1))
        @test dawson(Float32(1)) ≈ 0.53807950691276841914 rtol=2*eps(Float32)
        @test dawson(Float64(1)) ≈ 0.53807950691276841914 rtol=2*eps(Float64)
    end

    @testset "complex arguments" begin
        @test erf(ComplexF16(1+2im)) ≈ -0.53664356577856503399-5.0491437034470346695im
        @test erf(ComplexF32(1+2im)) ≈ -0.53664356577856503399-5.0491437034470346695im
        @test erf(ComplexF64(1+2im)) ≈ -0.53664356577856503399-5.0491437034470346695im

        @test erfc(ComplexF16(1+2im)) ≈ 1.5366435657785650340+5.0491437034470346695im
        @test erfc(ComplexF32(1+2im)) ≈ 1.5366435657785650340+5.0491437034470346695im
        @test erfc(ComplexF64(1+2im)) ≈ 1.5366435657785650340+5.0491437034470346695im

        @test_throws MethodError erfcx(ComplexF16(1))
        @test erfcx(ComplexF32(1+2im)) ≈ 0.14023958136627794370-0.22221344017989910261im
        @test erfcx(ComplexF64(1+2im)) ≈ 0.14023958136627794370-0.22221344017989910261im

        @test_throws MethodError erfi(ComplexF16(1))
        @test erfi(ComplexF32(1+2im)) ≈ -0.011259006028815025076+1.0036063427256517509im
        @test erfi(ComplexF64(1+2im)) ≈ -0.011259006028815025076+1.0036063427256517509im

        @test_throws MethodError erfinv(Complex(1))

        @test_throws MethodError erfcinv(Complex(1))

        @test_throws MethodError dawson(ComplexF16(1))
        @test dawson(ComplexF32(1+2im)) ≈ -13.388927316482919244-11.828715103889593303im
        @test dawson(ComplexF64(1+2im)) ≈ -13.388927316482919244-11.828715103889593303im
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

        @test logerfc(BigFloat(1000, 100)) ≈ -1.0000074801207219e6 rtol=2*eps(Float64)
        @test isnan(logerfc(BigFloat(NaN)))

        @test_throws MethodError erfi(big(1.0))

        @test_throws MethodError erfinv(BigFloat(1))

        @test_throws MethodError erfcinv(BigFloat(1))

        @test_throws MethodError dawson(BigFloat(1))
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
    end
end
