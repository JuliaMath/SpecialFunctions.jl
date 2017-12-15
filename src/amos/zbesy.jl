function ZBESY(ZR::Float64,ZI::Float64,FNU::Float64,KODE::Int32,N::Int32,CYR::AbstractArray{Float64},CYI::AbstractArray{Float64},NZ::Int32,CWRKR::AbstractArray{Float64},CWRKI::AbstractArray{Float64},IERR::Int32)
    AA::Float64 = zero(Float64)
    ASCLE::Float64 = zero(Float64)
    ATOL::Float64 = zero(Float64)
    BB::Float64 = zero(Float64)
    C1I::Float64 = zero(Float64)
    C1R::Float64 = zero(Float64)
    C2I::Float64 = zero(Float64)
    C2R::Float64 = zero(Float64)
    ELIM::Float64 = zero(Float64)
    EXI::Float64 = zero(Float64)
    EXR::Float64 = zero(Float64)
    EY::Float64 = zero(Float64)
    HCII::Float64 = zero(Float64)
    I::Int32 = zero(Int32)
    K::Int32 = zero(Int32)
    K1::Int32 = zero(Int32)
    K2::Int32 = zero(Int32)
    NZ1::Int32 = zero(Int32)
    NZ2::Int32 = zero(Int32)
    RTOL::Float64 = zero(Float64)
    STI::Float64 = zero(Float64)
    STR::Float64 = zero(Float64)
    TAY::Float64 = zero(Float64)
    TOL::Float64 = zero(Float64)
    IERR = Int32(0)
    NZ = Int32(0)
    if ZR == 0.0 && ZI == 0.0
        IERR = Int32(1)
    end
    if FNU < 0.0
        IERR = Int32(1)
    end
    if KODE < Int32(1) || KODE > Int32(2)
        IERR = Int32(1)
    end
    if N < Int32(1)
        IERR = Int32(1)
    end
    if IERR != Int32(0)
        return (NZ,IERR)
    end
    HCII = 0.5
    (NZ1,IERR) = ZBESH(ZR,ZI,FNU,KODE,Int32(1),N,CYR,CYI,NZ1,IERR)
    if IERR != Int32(0) && IERR != Int32(3)
        @goto line170
    end
    (NZ2,IERR) = ZBESH(ZR,ZI,FNU,KODE,Int32(2),N,CWRKR,CWRKI,NZ2,IERR)
    if IERR != Int32(0) && IERR != Int32(3)
        @goto line170
    end
    NZ = min(NZ1,NZ2)
    if KODE == Int32(2)
        @goto line60
    end
    for I = Int32(1):N
        STR = CWRKR[I] - CYR[I]
        STI = CWRKI[I] - CYI[I]
        CYR[I] = -STI * HCII
        CYI[I] = STR * HCII
        @label line50
    end
    return (NZ,IERR)
    @label line60
    TOL = max(D1MACH4,1.0e-18)
    K1 = I1MACH15
    K2 = I1MACH16
    K = min(abs(K1),abs(K2))
    R1M5 = D1MACH5
    ELIM = 2.303 * (DBLE(FLOAT(K)) * R1M5 - 3.0)
    EXR = cos(ZR)
    EXI = sin(ZR)
    EY = 0.0
    TAY = abs(ZI + ZI)
    if TAY < ELIM
        EY = exp(-TAY)
    end
    if ZI < 0.0
        @goto line90
    end
    C1R = EXR * EY
    C1I = EXI * EY
    C2R = EXR
    C2I = -EXI
    @label line70
    NZ = Int32(0)
    RTOL = 1.0 / TOL
    ASCLE = D1MACH1 * RTOL * 1000.0
    for I = Int32(1):N
        AA = CWRKR[I]
        BB = CWRKI[I]
        ATOL = 1.0
        if max(abs(AA),abs(BB)) > ASCLE
            @goto line75
        end
        AA = AA * RTOL
        BB = BB * RTOL
        ATOL = TOL
        @label line75
        STR = (AA * C2R - BB * C2I) * ATOL
        STI = (AA * C2I + BB * C2R) * ATOL
        AA = CYR[I]
        BB = CYI[I]
        ATOL = 1.0
        if max(abs(AA),abs(BB)) > ASCLE
            @goto line85
        end
        AA = AA * RTOL
        BB = BB * RTOL
        ATOL = TOL
        @label line85
        STR = STR - (AA * C1R - BB * C1I) * ATOL
        STI = STI - (AA * C1I + BB * C1R) * ATOL
        CYR[I] = -STI * HCII
        CYI[I] = STR * HCII
        if STR == 0.0 && STI == 0.0 && EY == 0.0
            NZ = NZ + Int32(1)
        end
        @label line80
    end
    return (NZ,IERR)
    @label line90
    C1R = EXR
    C1I = EXI
    C2R = EXR * EY
    C2I = -EXI * EY
    @goto line70
    @label line170
    NZ = Int32(0)
    return (NZ,IERR)
end
