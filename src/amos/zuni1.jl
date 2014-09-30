const _ZUNI1_CYR = Array(Float64,2)
const _ZUNI1_CYI = Array(Float64,2)
const _ZUNI1_CWRKR = Array(Float64,16)
const _ZUNI1_CWRKI = Array(Float64,16)
const _ZUNI1_CSSR = Array(Float64,3)
const _ZUNI1_CSRR = Array(Float64,3)
const _ZUNI1_BRY = Array(Float64,3)
function ZUNI1(ZR::Float64,ZI::Float64,FNU::Float64,KODE::Int32,N::Int32,YR::AbstractArray{Float64},YI::AbstractArray{Float64},NZ::Int32,NLAST::Int32,FNUL::Float64,TOL::Float64,ELIM::Float64,ALIM::Float64)
    APHI::Float64 = zero(Float64)
    ASCLE::Float64 = zero(Float64)
    const BRY = _ZUNI1_BRY
    C1R::Float64 = zero(Float64)
    C2I::Float64 = zero(Float64)
    C2M::Float64 = zero(Float64)
    C2R::Float64 = zero(Float64)
    CONER::Float64 = zero(Float64)
    CRSC::Float64 = zero(Float64)
    CSCL::Float64 = zero(Float64)
    const CSRR = _ZUNI1_CSRR
    const CSSR = _ZUNI1_CSSR
    const CWRKI = _ZUNI1_CWRKI
    const CWRKR = _ZUNI1_CWRKR
    const CYI = _ZUNI1_CYI
    const CYR = _ZUNI1_CYR
    FN::Float64 = zero(Float64)
    I::Int32 = zero(Int32)
    IFLAG::Int32 = zero(Int32)
    INIT::Int32 = zero(Int32)
    K::Int32 = zero(Int32)
    M::Int32 = zero(Int32)
    ND::Int32 = zero(Int32)
    NN::Int32 = zero(Int32)
    NUF::Int32 = zero(Int32)
    NW::Int32 = zero(Int32)
    PHII::Float64 = zero(Float64)
    PHIR::Float64 = zero(Float64)
    RAST::Float64 = zero(Float64)
    RS1::Float64 = zero(Float64)
    RZI::Float64 = zero(Float64)
    RZR::Float64 = zero(Float64)
    S1I::Float64 = zero(Float64)
    S1R::Float64 = zero(Float64)
    S2I::Float64 = zero(Float64)
    S2R::Float64 = zero(Float64)
    STI::Float64 = zero(Float64)
    STR::Float64 = zero(Float64)
    SUMI::Float64 = zero(Float64)
    SUMR::Float64 = zero(Float64)
    ZEROI::Float64 = zero(Float64)
    ZEROR::Float64 = zero(Float64)
    ZETA1I::Float64 = zero(Float64)
    ZETA1R::Float64 = zero(Float64)
    ZETA2I::Float64 = zero(Float64)
    ZETA2R::Float64 = zero(Float64)
    begin 
        ZEROR = 0.0
        ZEROI = 0.0
        CONER = 1.0
    end
    NZ = int32(0)
    ND = N
    NLAST = int32(0)
    CSCL = 1.0 / TOL
    CRSC = TOL
    CSSR[int32(1)] = CSCL
    CSSR[int32(2)] = CONER
    CSSR[int32(3)] = CRSC
    CSRR[int32(1)] = CRSC
    CSRR[int32(2)] = CONER
    CSRR[int32(3)] = CSCL
    BRY[int32(1)] = (1000.0D1MACH1) / TOL
    FN = DMAX1(FNU,1.0)
    INIT = int32(0)
    (INIT,PHIR,PHII,ZETA1R,ZETA1I,ZETA2R,ZETA2I,SUMR,SUMI) = ZUNIK(ZR,ZI,FN,int32(1),int32(1),TOL,INIT,PHIR,PHII,ZETA1R,ZETA1I,ZETA2R,ZETA2I,SUMR,SUMI,CWRKR,CWRKI)
    if KODE == int32(1)
        @goto line10
    end
    STR = ZR + ZETA2R
    STI = ZI + ZETA2I
    RAST = FN / ZABS(COMPLEX(STR,STI))
    STR = STR * RAST * RAST
    STI = -STI * RAST * RAST
    S1R = -ZETA1R + STR
    S1I = -ZETA1I + STI
    @goto line20
    @label line10
    S1R = -ZETA1R + ZETA2R
    S1I = -ZETA1I + ZETA2I
    @label line20
    RS1 = S1R
    if DABS(RS1) > ELIM
        @goto line130
    end
    @label line30
    NN = MIN0(int32(2),ND)
    for I = int32(1):NN
        FN = FNU + DBLE(FLOAT(ND - I))
        INIT = int32(0)
        (INIT,PHIR,PHII,ZETA1R,ZETA1I,ZETA2R,ZETA2I,SUMR,SUMI) = ZUNIK(ZR,ZI,FN,int32(1),int32(0),TOL,INIT,PHIR,PHII,ZETA1R,ZETA1I,ZETA2R,ZETA2I,SUMR,SUMI,CWRKR,CWRKI)
        if KODE == int32(1)
            @goto line40
        end
        STR = ZR + ZETA2R
        STI = ZI + ZETA2I
        RAST = FN / ZABS(COMPLEX(STR,STI))
        STR = STR * RAST * RAST
        STI = -STI * RAST * RAST
        S1R = -ZETA1R + STR
        S1I = -ZETA1I + STI + ZI
        @goto line50
        @label line40
        S1R = -ZETA1R + ZETA2R
        S1I = -ZETA1I + ZETA2I
        @label line50
        RS1 = S1R
        if DABS(RS1) > ELIM
            @goto line110
        end
        if I == int32(1)
            IFLAG = int32(2)
        end
        if DABS(RS1) < ALIM
            @goto line60
        end
        APHI = ZABS(COMPLEX(PHIR,PHII))
        RS1 = RS1 + DLOG(APHI)
        if DABS(RS1) > ELIM
            @goto line110
        end
        if I == int32(1)
            IFLAG = int32(1)
        end
        if RS1 < 0.0
            @goto line60
        end
        if I == int32(1)
            IFLAG = int32(3)
        end
        @label line60
        S2R = PHIR * SUMR - PHII * SUMI
        S2I = PHIR * SUMI + PHII * SUMR
        STR = DEXP(S1R) * CSSR[IFLAG]
        S1R = STR * DCOS(S1I)
        S1I = STR * DSIN(S1I)
        STR = S2R * S1R - S2I * S1I
        S2I = S2R * S1I + S2I * S1R
        S2R = STR
        if IFLAG != int32(1)
            @goto line70
        end
        (NW,) = ZUCHK(S2R,S2I,NW,BRY[int32(1)],TOL)
        if NW != int32(0)
            @goto line110
        end
        @label line70
        CYR[I] = S2R
        CYI[I] = S2I
        M = (ND - I) + int32(1)
        YR[M] = S2R * CSRR[IFLAG]
        YI[M] = S2I * CSRR[IFLAG]
        @label line80
    end
    if ND <= int32(2)
        @goto line100
    end
    RAST = 1.0 / ZABS(COMPLEX(ZR,ZI))
    STR = ZR * RAST
    STI = -ZI * RAST
    RZR = (STR + STR) * RAST
    RZI = (STI + STI) * RAST
    BRY[int32(2)] = 1.0 / BRY[int32(1)]
    BRY[int32(3)] = D1MACH2
    S1R = CYR[int32(1)]
    S1I = CYI[int32(1)]
    S2R = CYR[int32(2)]
    S2I = CYI[int32(2)]
    C1R = CSRR[IFLAG]
    ASCLE = BRY[IFLAG]
    K = ND - int32(2)
    FN = DBLE(FLOAT(K))
    for I = int32(3):ND
        C2R = S2R
        C2I = S2I
        S2R = S1R + (FNU + FN) * (RZR * C2R - RZI * C2I)
        S2I = S1I + (FNU + FN) * (RZR * C2I + RZI * C2R)
        S1R = C2R
        S1I = C2I
        C2R = S2R * C1R
        C2I = S2I * C1R
        YR[K] = C2R
        YI[K] = C2I
        K = K - int32(1)
        FN = FN - 1.0
        if IFLAG >= int32(3)
            @goto line90
        end
        STR = DABS(C2R)
        STI = DABS(C2I)
        C2M = DMAX1(STR,STI)
        if C2M <= ASCLE
            @goto line90
        end
        IFLAG = IFLAG + int32(1)
        ASCLE = BRY[IFLAG]
        S1R = S1R * C1R
        S1I = S1I * C1R
        S2R = C2R
        S2I = C2I
        S1R = S1R * CSSR[IFLAG]
        S1I = S1I * CSSR[IFLAG]
        S2R = S2R * CSSR[IFLAG]
        S2I = S2I * CSSR[IFLAG]
        C1R = CSRR[IFLAG]
        @label line90
    end
    @label line100
    return (NZ,NLAST)
    @label line110
    if RS1 > 0.0
        @goto line120
    end
    YR[ND] = ZEROR
    YI[ND] = ZEROI
    NZ = NZ + int32(1)
    ND = ND - int32(1)
    if ND == int32(0)
        @goto line100
    end
    (NUF,) = ZUOIK(ZR,ZI,FNU,KODE,int32(1),ND,YR,YI,NUF,TOL,ELIM,ALIM)
    if NUF < int32(0)
        @goto line120
    end
    ND = ND - NUF
    NZ = NZ + NUF
    if ND == int32(0)
        @goto line100
    end
    FN = FNU + DBLE(FLOAT(ND - int32(1)))
    if FN >= FNUL
        @goto line30
    end
    NLAST = ND
    return (NZ,NLAST)
    @label line120
    NZ = int32(-1)
    return (NZ,NLAST)
    @label line130
    if RS1 > 0.0
        @goto line120
    end
    NZ = N
    for I = int32(1):N
        YR[I] = ZEROR
        YI[I] = ZEROI
        @label line140
    end
    return (NZ,NLAST)
end
