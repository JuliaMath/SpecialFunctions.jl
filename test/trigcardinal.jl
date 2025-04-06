@testset "Cardinal trigonometric" begin
    @testset "sincu (unnormalized sinc)" begin
        a = 1.0
        @test sincu(a) ≈ sin(a)/a
        @test sincu(-a) ≈ sin(-a)/-a

        b = complex(2.0, 3.0)
        @test sincu(b) ≈ sin(b)/b
        @test sincu(-b) ≈ sin(-b)/-b

        # Right below threshold
        c = sqrt(sqrt(eps(Float64)))
        @test sincu(c) == 1 - c^2/6
        d = complex(0,c)
        @test sincu(d) == 1 - d^2/6

        # Limits/nans
        @test sincu(0) == 1.0
        @test sincu(Inf) == 0.0
        @test sincu(-Inf) == 0.0
        @test isnan(sincu(NaN))
    end

    @testset "sinhcu (unnormalized sinhc)" begin
        a = 1.0
        @test sinhcu(a) ≈ sinh(a)/a
        @test sinhcu(-a) ≈ sinh(-a)/-a

        b = complex(2.0, 3.0)
        @test sinhcu(b) ≈ sinh(b)/b
        @test sinhcu(-b) ≈ sinh(-b)/-b

        # 0th order approximation
        c = sqrt(eps(Float64))-eps(Float64)
        @test sinhcu(c) == 1.0
        d = complex(0,c)
        @test sinhcu(d) == 1.0

        # 2nd order approximation
        e = sqrt(eps(Float64))
        @test sinhcu(e) == 1 + e^2/6
        f = complex(0,e)
        @test sinhcu(f) == 1 + f^2/6

        # 4th order approximation
        g = sqrt(sqrt(eps(Float64)))
        @test sinhcu(g) == 1 + g^2/6 + g^4/120
        h = complex(0,g)
        @test sinhcu(h) == 1 + h^2/6 + h^4/120

        # Limits/nans
        @test sinhcu(0.0) == 1.0
        @test sinhcu(Inf*im) == 0.0
        @test sinhcu(-Inf*im) == 0.0
        @test isnan(sinhcu(NaN))
  end
end
