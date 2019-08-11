@testset "error functions" begin
    @testset "real argument" begin
        @test erf(Float16(1))   ≈ 0.84270079294971486934 rtol=2*eps(Float16)
        @test erf(Float32(1))   ≈ 0.84270079294971486934 rtol=2*eps(Float32)
        @test erf(Float64(1))   ≈ 0.84270079294971486934 rtol=2*eps(Float64)

        @test erfc(Float16(1))  ≈ 0.15729920705028513066 rtol=2*eps(Float16)
        @test erfc(Float32(1))  ≈ 0.15729920705028513066 rtol=2*eps(Float32)
        @test erfc(Float64(1))  ≈ 0.15729920705028513066 rtol=2*eps(Float64)

        @test_throws MethodError erfcx(Float16(1))
        @test erfcx(Float32(1)) ≈ 0.42758357615580700442    rtol=2*eps(Float32)
        @test erfcx(Float64(1)) ≈ 0.42758357615580700442    rtol=2*eps(Float64)

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
end
