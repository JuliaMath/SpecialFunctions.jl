@testset "Owen's T function" begin
    @test all(owens_t.(-3:0.1:3, 0) .== 0.)
    @test all(isapprox.(owens_t.(0, -3:0.1:3), 1 ./ (2 .* pi) .* atan.(-3:0.1:3), rtol=1e-6))
    @test all(isapprox.([owens_t(-h, a) for h in -3:0.1:3, a in -3:0.1:3], [owens_t(h, a) for h in -3:0.1:3, a in -3:0.1:3], rtol=1e-6))
    @test all(isapprox.([owens_t(h, -a) for h in -3:0.1:3, a in -3:0.1:3], [-owens_t(h, a) for h in -3:0.1:3, a in -3:0.1:3], rtol=1e-6))
end