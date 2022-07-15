using QuadGK
@testset "fresnel" begin
    # Generate random complex number
    z = randn(ComplexF64)
    # Test by comparing to numeric solution
    @test fresnelcos(z) ≈ quadgk(t->cos(π*t^2/2),0,z)[1]
    @test fresnelsin(z) ≈ quadgk(t->sin(π*t^2/2),0,z)[1]
    # Test just for code coverage 😄
    @test (fresnelcos(z),fresnelsin(z)) == fresnelsincos(z)

    # Generate random real number
    z = randn(Float64)
    # Test by comparing to numeric solution
    @test fresnelcos(z) ≈ quadgk(t->cos(π*t^2/2),0,z)[1]
    @test fresnelsin(z) ≈ quadgk(t->sin(π*t^2/2),0,z)[1]
    # Test just for code coverage 😄
    @test (fresnelcos(z),fresnelsin(z)) == fresnelsincos(z)

    # Precise values come from WolframAlpha calculator
    # One could add more decimals and more tests if needed

    @test fresnelsin(1.)      ≈ 0.4382591473903
    @test fresnelcos(1.)      ≈ 0.7798934003768
    @test fresnelsin(sqrt(2)*im) ≈ -0.7139722140219*im
    @test fresnelcos(sqrt(2)*im) ≈ 0.5288915951112*im
end