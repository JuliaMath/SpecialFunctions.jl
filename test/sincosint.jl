# Computed via wolframalpha.com: SinIntegral[SetPrecision[Table[x,{x, 1,20,1}],20]] and CosIntegral[SetPrecision[Table[x,{x, 1,20,1}],20]]

@testset "trigonometric integrals" begin
    sinintvals = [  0.9460830703671830149, 1.6054129768026934849, 1.848652527999468256,
                    1.75820313894905306, 1.54993124494467414, 1.4246875512805065,
                    1.4545966142480936, 1.5741868217069421, 1.665040075829602,
                    1.658347594218874, 1.578306806945727416, 1.504971241526373371,
                    1.499361722862824564, 1.556211050077665054, 1.618194443708368739,
                    1.631302268270032886, 1.590136415870701122, 1.536608096861185462,
                    1.518630031769363932, 1.548241701043439840]
    cosintvals = [0.3374039229009681347, 0.4229808287748649957, 0.119629786008000328,
                    -0.14098169788693041, -0.19002974965664388, -0.06805724389324713,
                    0.07669527848218452, 0.122433882532010, 0.0553475313331336,
                    -0.045456433004455, -0.08956313549547997948, -0.04978000688411367560,
                    0.02676412556403455504, 0.06939635592758454727, 0.04627867767436043960,
                    -0.01420019012019002240, -0.05524268226081385053, -0.04347510299950100478,
                    0.00515037100842612857, 0.04441982084535331654]
    for x in 1:20
        @test sinint(x)     ≅ sinintvals[x]
        @test sinint(-x)    ≅ -sinintvals[x]
        @test cosint(x)     ≅ cosintvals[x]
    end

    @test sinint(Float16(1)) == Float16(sinint(1))
    @test cosint(Float16(1)) == Float16(cosint(1))
    @test sinint(Float32(1)) == Float32(sinint(1))
    @test cosint(Float32(1)) == Float32(cosint(1))

    @test sinint(1//2) == sinint(0.5)
    @test cosint(1//2) == cosint(0.5)

    @test sinint(1e300)     ≅ π/2
    @test cosint(1e300)     ≅ -8.17881912115908554103E-301
    @test sinint(1e-300)    ≅ 1.0E-300
    @test cosint(1e-300)    ≅ -690.1983122333121

    @test sinint(Inf) == π/2
    @test cosint(Inf) == 0.0
    @test isnan(sinint(NaN))
    @test isnan(cosint(NaN))

    @test_throws ErrorException sinint(big(1.0))
    @test_throws ErrorException cosint(big(1.0))
    @test_throws DomainError cosint(-1.0)
    @test_throws DomainError cosint(Float32(-1.0))
    @test_throws DomainError cosint(Float16(-1.0))
    @test_throws DomainError cosint(-1//2)
    @test_throws MethodError sinint(Complex(1))
    @test_throws MethodError cosint(Complex(1))
end
