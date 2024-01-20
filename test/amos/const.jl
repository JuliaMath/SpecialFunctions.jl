
@testset "AMOS.I1_MACH,D1_MACH" begin
    # I1_MACH
    @test AMOS.I1_MACH[6] == sizeof(Int32)
    @test AMOS.I1_MACH[9] == typemax(Int32) 
    
    # D1_MACH
    @test AMOS.D1_MACH[3] == eps(Float64)/2
    @test AMOS.D1_MACH[4] == eps(Float64)
    @test AMOS.D1_MACH[5] == log10(2)
end
