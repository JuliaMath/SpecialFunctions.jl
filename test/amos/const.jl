
@testset "AMOS.I1_MACH,D1_MACH" begin
    # I1_MACH
    @test size(AMOS.I1_MACH) == (16,)
    @test AMOS.I1_MACH[6] == sizeof(Int32)
    @test AMOS.I1_MACH[9] == typemax(Int32) 
    
    # D1_MACH
    @test size(AMOS.D1_MACH) == (5,)
    @test AMOS.D1_MACH[3] == eps(Float64)/2
    @test AMOS.D1_MACH[4] == eps(Float64)
    @test AMOS.D1_MACH[5] == log10(2)
end

@testset "AMOS.GAMMALN_GLN" begin
    @test size(AMOS.GAMMALN_GLN) == (100,)
    for idx in 1:100
        @test AMOS.GAMMALN_GLN[idx] == log(gamma(idx))
    end
    
    @test size(AMOS.GAMMALN_CF) == (22,)
    # TODO:  test AMOS.GAMMALN_CF
end
