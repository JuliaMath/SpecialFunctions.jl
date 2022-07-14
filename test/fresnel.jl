@testset "fresnel" begin
    # Precise values come from WolframAlpha calculator
    # One could add more decimals and more tests if needed

    @test fresnels(1.)      ≈ 0.4382591473903
    @test fresnelc(1.)      ≈ 0.7798934003768
    @test fresnels(sqrt(2)) ≈ 0.7139722140219
    @test fresnelc(sqrt(2)) ≈ 0.5288915951112
    
    z = rand(ComplexF64)
    @test fresnel(z) == (fresnelc(z),fresnels(z))
end