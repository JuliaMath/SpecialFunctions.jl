function ZBESJ(ZR::Float64,ZI::Float64,FNU::Float64,KODE::Integer,N::Integer,CYR::AbstractArray{Float64},CYI::AbstractArray{Float64},NZ::Integer,IERR::Integer)
    AA::Float64 = 0
    ALIM::Float64 = 0
    ARG::Float64 = 0
    ASCLE::Float64 = 0
    ATOL::Float64 = 0
    AZ::Float64 = 0
    BB::Float64 = 0
    CII::Float64 = 0
    CSGNI::Float64 = 0
    CSGNR::Float64 = 0
    DIG::Float64 = 0
    ELIM::Float64 = 0
    FN::Float64 = 0
    FNUL::Float64 = 0
    HPI::Float64 = 0
    I::Int32 = 0
    INU::Int32 = 0
    INUH::Int32 = 0
    IR::Int32 = 0
    K::Int32 = 0
    K1::Int32 = 0
    K2::Int32 = 0
    NL::Int32 = 0
    R1M5::Float64 = 0
    RL::Float64 = 0
    RTOL::Float64 = 0
    STI::Float64 = 0
    STR::Float64 = 0
    TOL::Float64 = 0
    ZNI::Float64 = 0
    ZNR::Float64 = 0
    begin 
        HPI = 1.5707963267948966
    end
    IERR = 0
    NZ = 0
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
    TOL = DMAX1(D1MACH4,1.0e-18)
    K1 = I1MACH15
    K2 = I1MACH16
    R1M5 = D1MACH5
    K = MIN0(IABS(K1),IABS(K2))
    ELIM = 2.303 * (DBLE(FLOAT(K)) * R1M5 - 3.0)
    K1 = I1MACH14 - 1
    AA = R1M5 * DBLE(FLOAT(K1))
    DIG = DMIN1(AA,18.0)
    AA = AA * 2.303
    ALIM = ELIM + DMAX1(-AA,-41.45)
    RL = 1.2DIG + 3.0
    FNUL = 10.0 + 6.0 * (DIG - 3.0)
    AZ = ZABS(COMPLEX(ZR,ZI))
    FN = FNU + DBLE(FLOAT(N - 1))
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
        IERR = 3
    end
    if FN > AA
        IERR = 3
    end
    CII = 1.0
    INU = INT(SNGL(FNU))
    INUH = div(INU,2)
    IR = INU - 2INUH
    ARG = (FNU - DBLE(FLOAT(INU - IR))) * HPI
    CSGNR = DCOS(ARG)
    CSGNI = DSIN(ARG)
    if MOD(INUH,2) == 0
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
    if NZ < 0
        @goto line130
    end
    NL = N - NZ
    if NL == 0
        return (NZ,IERR)
    end
    RTOL = 1.0 / TOL
    ASCLE = D1MACH1 * RTOL * 1000.0
    for I = 1:NL
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
    if NZ == -2
        @goto line140
    end
    NZ = 0
    IERR = 2
    return (NZ,IERR)
    @label line140
    NZ = 0
    IERR = 5
    return (NZ,IERR)
    @label line260
    NZ = 0
    IERR = 4
    return (NZ,IERR)
end
