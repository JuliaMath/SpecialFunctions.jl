function ZACON(ZR::Float64,ZI::Float64,FNU::Float64,KODE::Integer,MR::Integer,N::Integer,YR::AbstractArray{Float64},YI::AbstractArray{Float64},NZ::Integer,RL::Float64,FNUL::Float64,TOL::Float64,ELIM::Float64,ALIM::Float64)
    ARG::Float64 = 0
    AS2::Float64 = 0
    ASCLE::Float64 = 0
    AZN::Float64 = 0
    BRY = Array(Float64,3)
    BSCLE::Float64 = 0
    C1I::Float64 = 0
    C1M::Float64 = 0
    C1R::Float64 = 0
    C2I::Float64 = 0
    C2R::Float64 = 0
    CKI::Float64 = 0
    CKR::Float64 = 0
    CONER::Float64 = 0
    CPN::Float64 = 0
    CSCL::Float64 = 0
    CSCR::Float64 = 0
    CSGNI::Float64 = 0
    CSGNR::Float64 = 0
    CSPNI::Float64 = 0
    CSPNR::Float64 = 0
    CSR::Float64 = 0
    CSRR = Array(Float64,3)
    CSSR = Array(Float64,3)
    CYI = Array(Float64,2)
    CYR = Array(Float64,2)
    FMR::Float64 = 0
    FN::Float64 = 0
    I::Int32 = 0
    INU::Int32 = 0
    IUF::Int32 = 0
    KFLAG::Int32 = 0
    NN::Int32 = 0
    NW::Int32 = 0
    PI::Float64 = 0
    PTI::Float64 = 0
    PTR::Float64 = 0
    RAZN::Float64 = 0
    RZI::Float64 = 0
    RZR::Float64 = 0
    S1I::Float64 = 0
    S1R::Float64 = 0
    S2I::Float64 = 0
    S2R::Float64 = 0
    SC1I::Float64 = 0
    SC1R::Float64 = 0
    SC2I::Float64 = 0
    SC2R::Float64 = 0
    SGN::Float64 = 0
    SPN::Float64 = 0
    STI::Float64 = 0
    STR::Float64 = 0
    YY::Float64 = 0
    ZEROR::Float64 = 0
    ZNI::Float64 = 0
    ZNR::Float64 = 0
    begin 
        PI = 3.141592653589793
    end
    begin 
        ZEROR = 0.0
        CONER = 1.0
    end
    NZ = 0
    ZNR = -ZR
    ZNI = -ZI
    NN = N
    (NW,) = ZBINU(ZNR,ZNI,FNU,KODE,NN,YR,YI,NW,RL,FNUL,TOL,ELIM,ALIM)
    if NW < 0
        @goto line90
    end
    NN = MIN0(2,N)
    (NW,) = ZBKNU(ZNR,ZNI,FNU,KODE,NN,CYR,CYI,NW,TOL,ELIM,ALIM)
    if NW != 0
        @goto line90
    end
    S1R = CYR[1]
    S1I = CYI[1]
    FMR = DBLE(FLOAT(MR))
    SGN = -(DSIGN(PI,FMR))
    CSGNR = ZEROR
    CSGNI = SGN
    if KODE == 1
        @goto line10
    end
    YY = -ZNI
    CPN = DCOS(YY)
    SPN = DSIN(YY)
    (CSGNR,CSGNI) = ZMLT(CSGNR,CSGNI,CPN,SPN,CSGNR,CSGNI)
    @label line10
    INU = INT(SNGL(FNU))
    ARG = (FNU - DBLE(FLOAT(INU))) * SGN
    CPN = DCOS(ARG)
    SPN = DSIN(ARG)
    CSPNR = CPN
    CSPNI = SPN
    if MOD(INU,2) == 0
        @goto line20
    end
    CSPNR = -CSPNR
    CSPNI = -CSPNI
    @label line20
    IUF = 0
    C1R = S1R
    C1I = S1I
    C2R = YR[1]
    C2I = YI[1]
    ASCLE = (1000.0D1MACH1) / TOL
    if KODE == 1
        @goto line30
    end
    (C1R,C1I,C2R,C2I,NW,IUF) = ZS1S2(ZNR,ZNI,C1R,C1I,C2R,C2I,NW,ASCLE,ALIM,IUF)
    NZ = NZ + NW
    SC1R = C1R
    SC1I = C1I
    @label line30
    (STR,STI) = ZMLT(CSPNR,CSPNI,C1R,C1I,STR,STI)
    (PTR,PTI) = ZMLT(CSGNR,CSGNI,C2R,C2I,PTR,PTI)
    YR[1] = STR + PTR
    YI[1] = STI + PTI
    if N == 1
        return NZ
    end
    CSPNR = -CSPNR
    CSPNI = -CSPNI
    S2R = CYR[2]
    S2I = CYI[2]
    C1R = S2R
    C1I = S2I
    C2R = YR[2]
    C2I = YI[2]
    if KODE == 1
        @goto line40
    end
    (C1R,C1I,C2R,C2I,NW,IUF) = ZS1S2(ZNR,ZNI,C1R,C1I,C2R,C2I,NW,ASCLE,ALIM,IUF)
    NZ = NZ + NW
    SC2R = C1R
    SC2I = C1I
    @label line40
    (STR,STI) = ZMLT(CSPNR,CSPNI,C1R,C1I,STR,STI)
    (PTR,PTI) = ZMLT(CSGNR,CSGNI,C2R,C2I,PTR,PTI)
    YR[2] = STR + PTR
    YI[2] = STI + PTI
    if N == 2
        return NZ
    end
    CSPNR = -CSPNR
    CSPNI = -CSPNI
    AZN = ZABS(COMPLEX(ZNR,ZNI))
    RAZN = 1.0 / AZN
    STR = ZNR * RAZN
    STI = -ZNI * RAZN
    RZR = (STR + STR) * RAZN
    RZI = (STI + STI) * RAZN
    FN = FNU + 1.0
    CKR = FN * RZR
    CKI = FN * RZI
    CSCL = 1.0 / TOL
    CSCR = TOL
    CSSR[1] = CSCL
    CSSR[2] = CONER
    CSSR[3] = CSCR
    CSRR[1] = CSCR
    CSRR[2] = CONER
    CSRR[3] = CSCL
    BRY[1] = ASCLE
    BRY[2] = 1.0 / ASCLE
    BRY[3] = D1MACH2
    AS2 = ZABS(COMPLEX(S2R,S2I))
    KFLAG = 2
    if AS2 > BRY[1]
        @goto line50
    end
    KFLAG = 1
    @goto line60
    @label line50
    if AS2 < BRY[2]
        @goto line60
    end
    KFLAG = 3
    @label line60
    BSCLE = BRY[KFLAG]
    S1R = S1R * CSSR[KFLAG]
    S1I = S1I * CSSR[KFLAG]
    S2R = S2R * CSSR[KFLAG]
    S2I = S2I * CSSR[KFLAG]
    CSR = CSRR[KFLAG]
    for I = 3:N
        STR = S2R
        STI = S2I
        S2R = (CKR * STR - CKI * STI) + S1R
        S2I = CKR * STI + CKI * STR + S1I
        S1R = STR
        S1I = STI
        C1R = S2R * CSR
        C1I = S2I * CSR
        STR = C1R
        STI = C1I
        C2R = YR[I]
        C2I = YI[I]
        if KODE == 1
            @goto line70
        end
        if IUF < 0
            @goto line70
        end
        (C1R,C1I,C2R,C2I,NW,IUF) = ZS1S2(ZNR,ZNI,C1R,C1I,C2R,C2I,NW,ASCLE,ALIM,IUF)
        NZ = NZ + NW
        SC1R = SC2R
        SC1I = SC2I
        SC2R = C1R
        SC2I = C1I
        if IUF != 3
            @goto line70
        end
        IUF = -4
        S1R = SC1R * CSSR[KFLAG]
        S1I = SC1I * CSSR[KFLAG]
        S2R = SC2R * CSSR[KFLAG]
        S2I = SC2I * CSSR[KFLAG]
        STR = SC2R
        STI = SC2I
        @label line70
        PTR = CSPNR * C1R - CSPNI * C1I
        PTI = CSPNR * C1I + CSPNI * C1R
        YR[I] = (PTR + CSGNR * C2R) - CSGNI * C2I
        YI[I] = PTI + CSGNR * C2I + CSGNI * C2R
        CKR = CKR + RZR
        CKI = CKI + RZI
        CSPNR = -CSPNR
        CSPNI = -CSPNI
        if KFLAG >= 3
            @goto line80
        end
        PTR = DABS(C1R)
        PTI = DABS(C1I)
        C1M = DMAX1(PTR,PTI)
        if C1M <= BSCLE
            @goto line80
        end
        KFLAG = KFLAG + 1
        BSCLE = BRY[KFLAG]
        S1R = S1R * CSR
        S1I = S1I * CSR
        S2R = STR
        S2I = STI
        S1R = S1R * CSSR[KFLAG]
        S1I = S1I * CSSR[KFLAG]
        S2R = S2R * CSSR[KFLAG]
        S2I = S2I * CSSR[KFLAG]
        CSR = CSRR[KFLAG]
        @label line80
    end
    return NZ
    @label line90
    NZ = -1
    if NW == -2
        NZ = -2
    end
    return NZ
end
