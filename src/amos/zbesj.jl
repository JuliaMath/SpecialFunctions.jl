function ZBESJ(ZR::Float64,ZI::Float64,FNU::Float64,KODE::Int32,N::Int32,CYR::AbstractArray{Float64},CYI::AbstractArray{Float64},NZ::Int32,IERR::Int32)
    AA::Float64 = zero(Float64)
    ALIM::Float64 = zero(Float64)
    ARG::Float64 = zero(Float64)
    ASCLE::Float64 = zero(Float64)
    ATOL::Float64 = zero(Float64)
    AZ::Float64 = zero(Float64)
    BB::Float64 = zero(Float64)
    CII::Float64 = zero(Float64)
    CSGNI::Float64 = zero(Float64)
    CSGNR::Float64 = zero(Float64)
    DIG::Float64 = zero(Float64)
    ELIM::Float64 = zero(Float64)
    FN::Float64 = zero(Float64)
    FNUL::Float64 = zero(Float64)
    HPI::Float64 = zero(Float64)
    I::Int32 = zero(Int32)
    INU::Int32 = zero(Int32)
    INUH::Int32 = zero(Int32)
    IR::Int32 = zero(Int32)
    K::Int32 = zero(Int32)
    K1::Int32 = zero(Int32)
    K2::Int32 = zero(Int32)
    NL::Int32 = zero(Int32)
    R1M5::Float64 = zero(Float64)
    RL::Float64 = zero(Float64)
    RTOL::Float64 = zero(Float64)
    STI::Float64 = zero(Float64)
    STR::Float64 = zero(Float64)
    TOL::Float64 = zero(Float64)
    ZNI::Float64 = zero(Float64)
    ZNR::Float64 = zero(Float64)
    begin 
        HPI = 1.5707963267948966
    end
    IERR = int32(0)
    NZ = int32(0)
    if FNU < 0.0
        IERR = int32(1)
    end
    if KODE < int32(1) || KODE > int32(2)
        IERR = int32(1)
    end
    if N < int32(1)
        IERR = int32(1)
    end
    if IERR != int32(0)
        return (NZ,IERR)
    end
    TOL = DMAX1(D1MACH4,1.0e-18)
    K1 = I1MACH15
    K2 = I1MACH16
    R1M5 = D1MACH5
    K = MIN0(IABS(K1),IABS(K2))
    ELIM = 2.303 * (DBLE(FLOAT(K)) * R1M5 - 3.0)
    K1 = I1MACH14 - int32(1)
    AA = R1M5 * DBLE(FLOAT(K1))
    DIG = DMIN1(AA,18.0)
    AA = AA * 2.303
    ALIM = ELIM + DMAX1(-AA,-41.45)
    RL = 1.2DIG + 3.0
    FNUL = 10.0 + 6.0 * (DIG - 3.0)
    AZ = ZABS(COMPLEX(ZR,ZI))
    FN = FNU + DBLE(FLOAT(N - int32(1)))
    AA = 0.5 / TOL
    BB = DBLE(FLOAT(I1MACH9)) * 0.5
    AA = DMIN1(AA,BB)
    if AZ > AA
        @goto line260
    end
    if FN > AA
        @goto line260
    end
    AA = DSQRT(AA)
    if AZ > AA
        IERR = int32(3)
    end
    if FN > AA
        IERR = int32(3)
    end
    CII = 1.0
    INU = INT(SNGL(FNU))
    INUH = div(INU,int32(2))
    IR = INU - int32(2) * INUH
    ARG = (FNU - DBLE(FLOAT(INU - IR))) * HPI
    CSGNR = DCOS(ARG)
    CSGNI = DSIN(ARG)
    if MOD(INUH,int32(2)) == int32(0)
        @goto line40
    end
    CSGNR = -CSGNR
    CSGNI = -CSGNI
    @label line40
    ZNR = ZI
    ZNI = -ZR
    if ZI >= 0.0
        @goto line50
    end
    ZNR = -ZNR
    ZNI = -ZNI
    CSGNI = -CSGNI
    CII = -CII
    @label line50
    (NZ,) = ZBINU(ZNR,ZNI,FNU,KODE,N,CYR,CYI,NZ,RL,FNUL,TOL,ELIM,ALIM)
    if NZ < int32(0)
        @goto line130
    end
    NL = N - NZ
    if NL == int32(0)
        return (NZ,IERR)
    end
    RTOL = 1.0 / TOL
    ASCLE = D1MACH1 * RTOL * 1000.0
    for I = int32(1):NL
        AA = CYR[I]
        BB = CYI[I]
        ATOL = 1.0
        if DMAX1(DABS(AA),DABS(BB)) > ASCLE
            @goto line55
        end
        AA = AA * RTOL
        BB = BB * RTOL
        ATOL = TOL
        @label line55
        STR = AA * CSGNR - BB * CSGNI
        STI = AA * CSGNI + BB * CSGNR
        CYR[I] = STR * ATOL
        CYI[I] = STI * ATOL
        STR = -CSGNI * CII
        CSGNI = CSGNR * CII
        CSGNR = STR
        @label line60
    end
    return (NZ,IERR)
    @label line130
    if NZ == int32(-2)
        @goto line140
    end
    NZ = int32(0)
    IERR = int32(2)
    return (NZ,IERR)
    @label line140
    NZ = int32(0)
    IERR = int32(5)
    return (NZ,IERR)
    @label line260
    NZ = int32(0)
    IERR = int32(4)
    return (NZ,IERR)
end
