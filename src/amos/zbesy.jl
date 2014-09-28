function ZBESY(ZR::Float64,ZI::Float64,FNU::Float64,KODE::Integer,N::Integer,CYR::AbstractArray{Float64},CYI::AbstractArray{Float64},NZ::Integer,CWRKR::AbstractArray{Float64},CWRKI::AbstractArray{Float64},IERR::Integer)
    AA::Float64 = 0
    ASCLE::Float64 = 0
    ATOL::Float64 = 0
    BB::Float64 = 0
    C1I::Float64 = 0
    C1R::Float64 = 0
    C2I::Float64 = 0
    C2R::Float64 = 0
    ELIM::Float64 = 0
    EXI::Float64 = 0
    EXR::Float64 = 0
    EY::Float64 = 0
    HCII::Float64 = 0
    I::Int32 = 0
    K::Int32 = 0
    K1::Int32 = 0
    K2::Int32 = 0
    NZ1::Int32 = 0
    NZ2::Int32 = 0
    RTOL::Float64 = 0
    STI::Float64 = 0
    STR::Float64 = 0
    TAY::Float64 = 0
    TOL::Float64 = 0
    IERR = 0
    NZ = 0
    if ZR == 0.0 && ZI == 0.0
        IERR = 1
    end
    if FNU < 0.0
        IERR = 1
    end
    if KODE < 1 || KODE > 2
        IERR = 1
    end
    if N < 1
        IERR = 1
    end
    if IERR != 0
        return (NZ,IERR)
    end
    HCII = 0.5
    (NZ1,IERR) = ZBESH(ZR,ZI,FNU,KODE,1,N,CYR,CYI,NZ1,IERR)
    if IERR != 0 && IERR != 3
        @goto line170
    end
    (NZ2,IERR) = ZBESH(ZR,ZI,FNU,KODE,2,N,CWRKR,CWRKI,NZ2,IERR)
    if IERR != 0 && IERR != 3
        @goto line170
    end
    NZ = MIN0(NZ1,NZ2)
    if KODE == 2
        @goto line60
    end
    for I = 1:N
        STR = CWRKR[I] - CYR[I]
        STI = CWRKI[I] - CYI[I]
        CYR[I] = -STI * HCII
        CYI[I] = STR * HCII
        @label line50
    end
    return (NZ,IERR)
    @label line60
    TOL = DMAX1(D1MACH4,1.0e-18)
    K1 = I1MACH15
    K2 = I1MACH16
    K = MIN0(IABS(K1),IABS(K2))
    R1M5 = D1MACH5
    ELIM = 2.303 * (DBLE(FLOAT(K)) * R1M5 - 3.0)
    EXR = DCOS(ZR)
    EXI = DSIN(ZR)
    EY = 0.0
    TAY = DABS(ZI + ZI)
    if TAY < ELIM
        EY = DEXP(-TAY)
    end
    if ZI < 0.0
        @goto line90
    end
    C1R = EXR * EY
    C1I = EXI * EY
    C2R = EXR
    C2I = -EXI
    @label line70
    NZ = 0
    RTOL = 1.0 / TOL
    ASCLE = D1MACH1 * RTOL * 1000.0
    for I = 1:N
        AA = CWRKR[I]
        BB = CWRKI[I]
        ATOL = 1.0
        if DMAX1(DABS(AA),DABS(BB)) > ASCLE
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
        if DMAX1(DABS(AA),DABS(BB)) > ASCLE
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
            NZ = NZ + 1
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
    NZ = 0
    return (NZ,IERR)
end
