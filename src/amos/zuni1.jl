function ZUNI1(ZR::Float64,ZI::Float64,FNU::Float64,KODE::Integer,N::Integer,YR::AbstractArray{Float64},YI::AbstractArray{Float64},NZ::Integer,NLAST::Integer,FNUL::Float64,TOL::Float64,ELIM::Float64,ALIM::Float64)
    APHI::Float64 = 0
    ASCLE::Float64 = 0
    BRY = Array(Float64,3)
    C1R::Float64 = 0
    C2I::Float64 = 0
    C2M::Float64 = 0
    C2R::Float64 = 0
    CONER::Float64 = 0
    CRSC::Float64 = 0
    CSCL::Float64 = 0
    CSRR = Array(Float64,3)
    CSSR = Array(Float64,3)
    CWRKI = Array(Float64,16)
    CWRKR = Array(Float64,16)
    CYI = Array(Float64,2)
    CYR = Array(Float64,2)
    FN::Float64 = 0
    I::Int32 = 0
    IFLAG::Int32 = 0
    INIT::Int32 = 0
    K::Int32 = 0
    M::Int32 = 0
    ND::Int32 = 0
    NN::Int32 = 0
    NUF::Int32 = 0
    NW::Int32 = 0
    PHII::Float64 = 0
    PHIR::Float64 = 0
    RAST::Float64 = 0
    RS1::Float64 = 0
    RZI::Float64 = 0
    RZR::Float64 = 0
    S1I::Float64 = 0
    S1R::Float64 = 0
    S2I::Float64 = 0
    S2R::Float64 = 0
    STI::Float64 = 0
    STR::Float64 = 0
    SUMI::Float64 = 0
    SUMR::Float64 = 0
    ZEROI::Float64 = 0
    ZEROR::Float64 = 0
    ZETA1I::Float64 = 0
    ZETA1R::Float64 = 0
    ZETA2I::Float64 = 0
    ZETA2R::Float64 = 0
    begin 
        ZEROR = 0.0
        ZEROI = 0.0
        CONER = 1.0
    end
    NZ = 0
    ND = N
    NLAST = 0
    CSCL = 1.0 / TOL
    CRSC = TOL
    CSSR[1] = CSCL
    CSSR[2] = CONER
    CSSR[3] = CRSC
    CSRR[1] = CRSC
    CSRR[2] = CONER
    CSRR[3] = CSCL
    BRY[1] = (1000.0D1MACH1) / TOL
    FN = DMAX1(FNU,1.0)
    INIT = 0
    (INIT,PHIR,PHII,ZETA1R,ZETA1I,ZETA2R,ZETA2I,SUMR,SUMI) = ZUNIK(ZR,ZI,FN,1,1,TOL,INIT,PHIR,PHII,ZETA1R,ZETA1I,ZETA2R,ZETA2I,SUMR,SUMI,CWRKR,CWRKI)
    if KODE == 1
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
    NN = MIN0(2,ND)
    for I = 1:NN
        FN = FNU + DBLE(FLOAT(ND - I))
        INIT = 0
        (INIT,PHIR,PHII,ZETA1R,ZETA1I,ZETA2R,ZETA2I,SUMR,SUMI) = ZUNIK(ZR,ZI,FN,1,0,TOL,INIT,PHIR,PHII,ZETA1R,ZETA1I,ZETA2R,ZETA2I,SUMR,SUMI,CWRKR,CWRKI)
        if KODE == 1
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
        if I == 1
            IFLAG = 2
        end
        if DABS(RS1) < ALIM
            @goto line60
        end
        APHI = ZABS(COMPLEX(PHIR,PHII))
        RS1 = RS1 + DLOG(APHI)
        if DABS(RS1) > ELIM
            @goto line110
        end
        if I == 1
            IFLAG = 1
        end
        if RS1 < 0.0
            @goto line60
        end
        if I == 1
            IFLAG = 3
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
        if IFLAG != 1
            @goto line70
        end
        (NW,) = ZUCHK(S2R,S2I,NW,BRY[1],TOL)
        if NW != 0
            @goto line110
        end
        @label line70
        CYR[I] = S2R
        CYI[I] = S2I
        M = (ND - I) + 1
        YR[M] = S2R * CSRR[IFLAG]
        YI[M] = S2I * CSRR[IFLAG]
        @label line80
    end
    if ND <= 2
        @goto line100
    end
    RAST = 1.0 / ZABS(COMPLEX(ZR,ZI))
    STR = ZR * RAST
    STI = -ZI * RAST
    RZR = (STR + STR) * RAST
    RZI = (STI + STI) * RAST
    BRY[2] = 1.0 / BRY[1]
    BRY[3] = D1MACH2
    S1R = CYR[1]
    S1I = CYI[1]
    S2R = CYR[2]
    S2I = CYI[2]
    C1R = CSRR[IFLAG]
    ASCLE = BRY[IFLAG]
    K = ND - 2
    FN = DBLE(FLOAT(K))
    for I = 3:ND
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
        K = K - 1
        FN = FN - 1.0
        if IFLAG >= 3
            @goto line90
        end
        STR = DABS(C2R)
        STI = DABS(C2I)
        C2M = DMAX1(STR,STI)
        if C2M <= ASCLE
            @goto line90
        end
        IFLAG = IFLAG + 1
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
    NZ = NZ + 1
    ND = ND - 1
    if ND == 0
        @goto line100
    end
    (NUF,) = ZUOIK(ZR,ZI,FNU,KODE,1,ND,YR,YI,NUF,TOL,ELIM,ALIM)
    if NUF < 0
        @goto line120
    end
    ND = ND - NUF
    NZ = NZ + NUF
    if ND == 0
        @goto line100
    end
    FN = FNU + DBLE(FLOAT(ND - 1))
    if FN >= FNUL
        @goto line30
    end
    NLAST = ND
    return (NZ,NLAST)
    @label line120
    NZ = -1
    return (NZ,NLAST)
    @label line130
    if RS1 > 0.0
        @goto line120
    end
    NZ = N
    for I = 1:N
        YR[I] = ZEROR
        YI[I] = ZEROI
        @label line140
    end
    return (NZ,NLAST)
end
