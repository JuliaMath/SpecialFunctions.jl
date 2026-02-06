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
        @test isnan(ellipk(NaN))
        @test ellipk(-Inf) == 0.0
        @test ellipk(-1e30) ≈ 0.0                          atol=1e-13
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
        @test isnan(ellipe(NaN))
        @test ellipe(-Inf) == Inf
        @test ellipe(-1e16) > ellipe(-1e15)
        @test_throws MethodError ellipe(BigFloat(-1))
        @test_throws DomainError ellipe(1.2)
    end
end

@testset "elliptic edge cases" begin
    # Boundary value tests at each polynomial table transition (idx = 0,1,2,...,9)
    # These test the boundaries where the binary tree branching changes paths
    @testset "ellipk boundary values" begin
        # Table boundaries: 0.0, 0.1, 0.2, ..., 0.9, 1.0
        # Reference values from current implementation (regression test baseline)
        @test ellipk(0.0)  ≈ 1.5707963267948966    rtol=2*eps()  # π/2
        @test ellipk(0.1)  ≈ 1.6124413487202192    rtol=2*eps()
        @test ellipk(0.2)  ≈ 1.6596235986105277    rtol=2*eps()
        @test ellipk(0.3)  ≈ 1.7138894481787912    rtol=2*eps()
        @test ellipk(0.4)  ≈ 1.7775193714912532    rtol=2*eps()
        @test ellipk(0.5)  ≈ 1.8540746773013721    rtol=2*eps()
        @test ellipk(0.6)  ≈ 1.9495677498060260    rtol=2*eps()
        @test ellipk(0.7)  ≈ 2.0753631352924691    rtol=2*eps()
        @test ellipk(0.8)  ≈ 2.2572053268208530    rtol=2*eps()
        @test ellipk(0.85) ≈ 2.3890164863255676    rtol=2*eps()  # Table 10-11 boundary
        @test ellipk(0.9)  ≈ 2.5780921133481733    rtol=2*eps()
        @test ellipk(0.95) ≈ 2.9083372484445515    rtol=2*eps()
        @test ellipk(0.99) ≈ 3.6956373629898747    rtol=2*eps()

        # Near-boundary tests using prevfloat/nextfloat
        @test ellipk(prevfloat(0.1)) ≈ ellipk(0.1) rtol=1e-10
        @test ellipk(nextfloat(0.1)) ≈ ellipk(0.1) rtol=1e-10
        @test ellipk(prevfloat(0.5)) ≈ ellipk(0.5) rtol=1e-10
        @test ellipk(nextfloat(0.5)) ≈ ellipk(0.5) rtol=1e-10
        @test ellipk(prevfloat(0.9)) ≈ ellipk(0.9) rtol=1e-10
        @test ellipk(nextfloat(0.9)) ≈ ellipk(0.9) rtol=1e-10
    end

    @testset "ellipe boundary values" begin
        # Table boundaries: 0.0, 0.1, 0.2, ..., 0.9, 1.0
        # Reference values from current implementation (regression test baseline)
        @test ellipe(0.0)  ≈ 1.5707963267948966    rtol=2*eps()  # π/2
        @test ellipe(0.1)  ≈ 1.5307576368977631    rtol=2*eps()
        @test ellipe(0.2)  ≈ 1.4890350580958527    rtol=2*eps()
        @test ellipe(0.3)  ≈ 1.4453630644126654    rtol=2*eps()
        @test ellipe(0.4)  ≈ 1.3993921388974322    rtol=2*eps()
        @test ellipe(0.5)  ≈ 1.3506438810476753    rtol=2*eps()
        @test ellipe(0.6)  ≈ 1.2984280350469133    rtol=2*eps()
        @test ellipe(0.7)  ≈ 1.2416705679458224    rtol=2*eps()
        @test ellipe(0.8)  ≈ 1.1784899243278386    rtol=2*eps()
        @test ellipe(0.85) ≈ 1.1433957918831656    rtol=2*eps()  # Table 10-11 boundary
        @test ellipe(0.9)  ≈ 1.1047747327040722    rtol=2*eps()
        @test ellipe(0.95) ≈ 1.0604737277662784    rtol=2*eps()
        @test ellipe(0.99) ≈ 1.0159935450252240    rtol=2*eps()

        # Near-boundary tests
        @test ellipe(prevfloat(0.1)) ≈ ellipe(0.1) rtol=1e-10
        @test ellipe(nextfloat(0.1)) ≈ ellipe(0.1) rtol=1e-10
        @test ellipe(prevfloat(0.5)) ≈ ellipe(0.5) rtol=1e-10
        @test ellipe(nextfloat(0.5)) ≈ ellipe(0.5) rtol=1e-10
        @test ellipe(prevfloat(0.9)) ≈ ellipe(0.9) rtol=1e-10
        @test ellipe(nextfloat(0.9)) ≈ ellipe(0.9) rtol=1e-10
    end

    @testset "negative m values" begin
        # Negative m uses transformation: x = m/(m-1), K(m) = K(x)/sqrt(1-m)
        # Reference values from current implementation (regression test baseline)
        @test ellipk(-0.1)  ≈ 1.5335928197134567   rtol=2*eps()
        @test ellipk(-0.5)  ≈ 1.4157372084259563   rtol=2*eps()
        @test ellipk(-1.0)  ≈ 1.3110287771460600   rtol=2*eps()
        @test ellipk(-2.0)  ≈ 1.1714200841467699   rtol=2*eps()
        @test ellipk(-10.0) ≈ 0.7908718902387385   rtol=2*eps()

        @test ellipe(-0.1)  ≈ 1.6093590249375296   rtol=2*eps()
        @test ellipe(-0.5)  ≈ 1.7517712756948178   rtol=2*eps()
        @test ellipe(-1.0)  ≈ 1.9100988945138559   rtol=2*eps()
        @test ellipe(-2.0)  ≈ 2.1844381427462012   rtol=2*eps()
        @test ellipe(-10.0) ≈ 3.6391380384177681   rtol=2*eps()
    end

    @testset "special edge cases" begin
        # Values very close to 1.0
        @test ellipk(prevfloat(1.0)) > 10.0  # Should be very large
        @test ellipe(prevfloat(1.0)) ≈ 1.0 rtol=1e-6  # Should be close to 1

        # Very small positive values
        @test ellipk(1e-15) ≈ π/2 rtol=1e-10
        @test ellipe(1e-15) ≈ π/2 rtol=1e-10

        # Consistency: E(m) ≤ π/2 for m ∈ [0,1]
        for m in 0.0:0.1:1.0
            @test ellipe(m) <= π/2 + eps()
        end

        # Consistency: K(m) ≥ π/2 for m ∈ [0,1)
        for m in 0.0:0.1:0.9
            @test ellipk(m) >= π/2 - eps()
        end
    end
end
