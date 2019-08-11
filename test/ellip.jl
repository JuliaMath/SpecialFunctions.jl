#Computed using Wolframalpha EllipticK and EllipticE functions.
@testset "elliptic integrals" begin
    @testset "complete of 1st kind" begin
        @test ellipk(-0.5)  ≈ 1.41573720842595619           rtol=2*eps()
        @test ellipk(0)     ≈ 1.570796326794896619231322    rtol=2*eps()
        @test ellipk(0.01)  ≈ 1.57474556151735595           rtol=2*eps()
        @test ellipk(0.17)  ≈ 1.6448064907988806            rtol=2*eps()
        @test ellipk(0.25)  ≈ 1.685750354812596             rtol=2*eps()
        @test ellipk(0.33)  ≈ 1.73186477825209792           rtol=2*eps()
        @test ellipk(0.45)  ≈ 1.81388393681698264           rtol=2*eps()
        @test ellipk(0.5)   ≈ 1.854074677301371918          rtol=2*eps()
        @test ellipk(0.69)  ≈ 2.0608816467301313            rtol=2*eps()
        @test ellipk(0.75)  ≈ 2.15651564749964323           rtol=2*eps()
        @test ellipk(0.84)  ≈ 2.3592635547450067            rtol=2*eps()
        @test ellipk(0.88)  ≈ 2.492635323239715564          rtol=2*eps()
        @test ellipk(0.92)  ≈ 2.683551406315229344          rtol=2*eps()
        @test ellipk(1.0) == Inf
        @test ellipk(Float16(0.92)) ≈ 2.683551406315229344  rtol=2*eps(Float16)
        @test ellipk(Float32(0.92)) ≈ 2.683551406315229344  rtol=2*eps(Float32)
        @test_throws MethodError ellipk(BigFloat(0.5))
        @test_throws DomainError ellipk(1.1)
    end

    @testset "complete of 2nd kind" begin
        @test ellipe(-0.1)  ≈ 1.6093590249375295            rtol=2*eps()
        @test ellipe(0)     ≈ 1.570796326794896619231322    rtol=2*eps()
        @test ellipe(0.01)  ≈ 1.5668619420216682            rtol=2*eps()
        @test ellipe(0.15)  ≈ 1.5101218320928197            rtol=2*eps()
        @test ellipe(0.21)  ≈ 1.4847605813318776            rtol=2*eps()
        @test ellipe(0.3)   ≈ 1.4453630644126652            rtol=2*eps()
        @test ellipe(0.42)  ≈ 1.3898829914929717            rtol=2*eps()
        @test ellipe(0.5)   ≈ 1.3506438810476755            rtol=2*eps()
        @test ellipe(0.66)  ≈ 1.2650125751607508            rtol=2*eps()
        @test ellipe(0.76)  ≈ 1.2047136418292115            rtol=2*eps()
        @test ellipe(0.8)   ≈ 1.17848992432783852           rtol=2*eps()
        @test ellipe(0.865) ≈ 1.1322436887003925            rtol=2*eps()
        @test ellipe(0.99)  ≈ 1.0159935450252239            rtol=2*eps()
        @test ellipe(1.0)   ≈ 1.00                          rtol=2*eps()
        @test ellipe(Float16(0.865)) ≈ 1.1322436887003925   rtol=2*eps(Float16)
        @test ellipe(Float32(0.865)) ≈ 1.1322436887003925   rtol=2*eps(Float32)
        @test_throws MethodError ellipe(BigFloat(-1))
        @test_throws DomainError ellipe(1.2)
    end
end
