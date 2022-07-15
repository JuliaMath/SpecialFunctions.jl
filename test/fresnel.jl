@testset "fresnel" begin
    # Precise values come from WolframAlpha calculator
    # One could add more decimals and more tests if needed

    @test fresnelsin(1.)      ≈ 0.4382591473903
    @test fresnelcos(1.)      ≈ 0.7798934003768
    @test fresnelsin(sqrt(2)) ≈ 0.7139722140219
    @test fresnelcos(sqrt(2)) ≈ 0.5288915951112
    
    z = rand(ComplexF64)
    @test fresnelsincos(z) == (fresnelcos(z),fresnelsin(z))
end