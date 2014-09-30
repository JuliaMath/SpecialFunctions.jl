const _ZUNK1_ZETA2R = Array(Float64,2)
const _ZUNK1_ZETA2I = Array(Float64,2)
const _ZUNK1_ZETA1R = Array(Float64,2)
const _ZUNK1_ZETA1I = Array(Float64,2)
const _ZUNK1_SUMR = Array(Float64,2)
const _ZUNK1_SUMI = Array(Float64,2)
const _ZUNK1_PHIR = Array(Float64,2)
const _ZUNK1_PHII = Array(Float64,2)
const _ZUNK1_INIT = Array(Int32,2)
const _ZUNK1_CYR = Array(Float64,2)
const _ZUNK1_CYI = Array(Float64,2)
const _ZUNK1_CWRKR = Array(Float64,16,3)
const _ZUNK1_CWRKI = Array(Float64,16,3)
const _ZUNK1_CSSR = Array(Float64,3)
const _ZUNK1_CSRR = Array(Float64,3)
const _ZUNK1_BRY = Array(Float64,3)
function ZUNK1(ZR::Float64,ZI::Float64,FNU::Float64,KODE::Integer,MR::Integer,N::Integer,YR::AbstractArray{Float64},YI::AbstractArray{Float64},NZ::Integer,TOL::Float64,ELIM::Float64,ALIM::Float64)
    ANG::Float64 = 0
    APHI::Float64 = 0
    ASC::Float64 = 0
    ASCLE::Float64 = 0
    const BRY = _ZUNK1_BRY
    C1I::Float64 = 0
    C1R::Float64 = 0
    C2I::Float64 = 0
    C2M::Float64 = 0
    C2R::Float64 = 0
    CKI::Float64 = 0
    CKR::Float64 = 0
    CONER::Float64 = 0
    CRSC::Float64 = 0
    CSCL::Float64 = 0
    CSGNI::Float64 = 0
    CSPNI::Float64 = 0
    CSPNR::Float64 = 0
    CSR::Float64 = 0
    const CSRR = _ZUNK1_CSRR
    const CSSR = _ZUNK1_CSSR
    const CWRKI = _ZUNK1_CWRKI
    const CWRKR = _ZUNK1_CWRKR
    const CYI = _ZUNK1_CYI
    const CYR = _ZUNK1_CYR
    FMR::Float64 = 0
    FN::Float64 = 0
    FNF::Float64 = 0
    I::Int32 = 0
    IB::Int32 = 0
    IC::Int32 = 0
    IFLAG::Int32 = 0
    IFN::Int32 = 0
    IL::Int32 = 0
    const INIT = _ZUNK1_INIT
    INITD::Int32 = 0
    INU::Int32 = 0
    IPARD::Int32 = 0
    IUF::Int32 = 0
    J::Int32 = 0
    K::Int32 = 0
    KDFLG::Int32 = 0
    KFLAG::Int32 = 0
    KK::Int32 = 0
    NW::Int32 = 0
    PHIDI::Float64 = 0
    PHIDR::Float64 = 0
    const PHII = _ZUNK1_PHII
    const PHIR = _ZUNK1_PHIR
    PI::Float64 = 0
    RAST::Float64 = 0
    RAZR::Float64 = 0
    RS1::Float64 = 0
    RZI::Float64 = 0
    RZR::Float64 = 0
    S1I::Float64 = 0
    S1R::Float64 = 0
    S2I::Float64 = 0
    S2R::Float64 = 0
    SGN::Float64 = 0
    STI::Float64 = 0
    STR::Float64 = 0
    SUMDI::Float64 = 0
    SUMDR::Float64 = 0
    const SUMI = _ZUNK1_SUMI
    const SUMR = _ZUNK1_SUMR
    ZEROI::Float64 = 0
    ZEROR::Float64 = 0
    ZET1DI::Float64 = 0
    ZET1DR::Float64 = 0
    ZET2DI::Float64 = 0
    ZET2DR::Float64 = 0
    const ZETA1I = _ZUNK1_ZETA1I
    const ZETA1R = _ZUNK1_ZETA1R
    const ZETA2I = _ZUNK1_ZETA2I
    const ZETA2R = _ZUNK1_ZETA2R
    ZRI::Float64 = 0
    ZRR::Float64 = 0
    begin 
        ZEROR = 0.0
        ZEROI = 0.0
        CONER = 1.0
    end
    begin 
        PI = 3.141592653589793
    end
    KDFLG = 1
    NZ = 0
    CSCL = 1.0 / TOL
    CRSC = TOL
    CSSR[1] = CSCL
    CSSR[2] = CONER
    CSSR[3] = CRSC
    CSRR[1] = CRSC
    CSRR[2] = CONER
    CSRR[3] = CSCL
    BRY[1] = (1000.0D1MACH1) / TOL
    BRY[2] = 1.0 / BRY[1]
    BRY[3] = D1MACH2
    ZRR = ZR
    ZRI = ZI
    if ZR >= 0.0
        @goto line10
    end
    ZRR = -ZR
    ZRI = -ZI
    @label line10
    J = 2
    for I = 1:N
        J = 3 - J
        FN = FNU + DBLE(FLOAT(I - 1))
        INIT[J] = 0
        (INIT[J],PHIR[J],PHII[J],ZETA1R[J],ZETA1I[J],ZETA2R[J],ZETA2I[J],SUMR[J],SUMI[J]) = ZUNIK(ZRR,ZRI,FN,2,0,TOL,INIT[J],PHIR[J],PHII[J],ZETA1R[J],ZETA1I[J],ZETA2R[J],ZETA2I[J],SUMR[J],SUMI[J],slice(CWRKR,1:16,int(J)),slice(CWRKI,1:16,int(J)))
        if KODE == 1
            @goto line20
        end
        STR = ZRR + ZETA2R[J]
        STI = ZRI + ZETA2I[J]
        RAST = FN / ZABS(COMPLEX(STR,STI))
        STR = STR * RAST * RAST
        STI = -STI * RAST * RAST
        S1R = ZETA1R[J] - STR
        S1I = ZETA1I[J] - STI
        @goto line30
        @label line20
        S1R = ZETA1R[J] - ZETA2R[J]
        S1I = ZETA1I[J] - ZETA2I[J]
        @label line30
        RS1 = S1R
        if DABS(RS1) > ELIM
            @goto line60
        end
        if KDFLG == 1
            KFLAG = 2
        end
        if DABS(RS1) < ALIM
            @goto line40
        end
        APHI = ZABS(COMPLEX(PHIR[J],PHII[J]))
        RS1 = RS1 + DLOG(APHI)
        if DABS(RS1) > ELIM
            @goto line60
        end
        if KDFLG == 1
            KFLAG = 1
        end
        if RS1 < 0.0
            @goto line40
        end
        if KDFLG == 1
            KFLAG = 3
        end
        @label line40
        S2R = PHIR[J] * SUMR[J] - PHII[J] * SUMI[J]
        S2I = PHIR[J] * SUMI[J] + PHII[J] * SUMR[J]
        STR = DEXP(S1R) * CSSR[KFLAG]
        S1R = STR * DCOS(S1I)
        S1I = STR * DSIN(S1I)
        STR = S2R * S1R - S2I * S1I
        S2I = S1R * S2I + S2R * S1I
        S2R = STR
        if KFLAG != 1
            @goto line50
        end
        (NW,) = ZUCHK(S2R,S2I,NW,BRY[1],TOL)
        if NW != 0
            @goto line60
        end
        @label line50
        CYR[KDFLG] = S2R
        CYI[KDFLG] = S2I
        YR[I] = S2R * CSRR[KFLAG]
        YI[I] = S2I * CSRR[KFLAG]
        if KDFLG == 2
            @goto line75
        end
        KDFLG = 2
        @goto line70
        @label line60
        if RS1 > 0.0
            @goto line300
        end
        if ZR < 0.0
            @goto line300
        end
        KDFLG = 1
        YR[I] = ZEROR
        YI[I] = ZEROI
        NZ = NZ + 1
        if I == 1
            @goto line70
        end
        if YR[I - 1] == ZEROR && YI[I - 1] == ZEROI
            @goto line70
        end
        YR[I - 1] = ZEROR
        YI[I - 1] = ZEROI
        NZ = NZ + 1
        @label line70
    end
    I = N
    @label line75
    RAZR = 1.0 / ZABS(COMPLEX(ZRR,ZRI))
    STR = ZRR * RAZR
    STI = -ZRI * RAZR
    RZR = (STR + STR) * RAZR
    RZI = (STI + STI) * RAZR
    CKR = FN * RZR
    CKI = FN * RZI
    IB = I + 1
    if N < IB
        @goto line160
    end
    FN = FNU + DBLE(FLOAT(N - 1))
    IPARD = 1
    if MR != 0
        IPARD = 0
    end
    INITD = 0
    (INITD,PHIDR,PHIDI,ZET1DR,ZET1DI,ZET2DR,ZET2DI,SUMDR,SUMDI) = ZUNIK(ZRR,ZRI,FN,2,IPARD,TOL,INITD,PHIDR,PHIDI,ZET1DR,ZET1DI,ZET2DR,ZET2DI,SUMDR,SUMDI,slice(CWRKR,1:16,int(3)),slice(CWRKI,1:16,int(3)))
    if KODE == 1
        @goto line80
    end
    STR = ZRR + ZET2DR
    STI = ZRI + ZET2DI
    RAST = FN / ZABS(COMPLEX(STR,STI))
    STR = STR * RAST * RAST
    STI = -STI * RAST * RAST
    S1R = ZET1DR - STR
    S1I = ZET1DI - STI
    @goto line90
    @label line80
    S1R = ZET1DR - ZET2DR
    S1I = ZET1DI - ZET2DI
    @label line90
    RS1 = S1R
    if DABS(RS1) > ELIM
        @goto line95
    end
    if DABS(RS1) < ALIM
        @goto line100
    end
    APHI = ZABS(COMPLEX(PHIDR,PHIDI))
    RS1 = RS1 + DLOG(APHI)
    if DABS(RS1) < ELIM
        @goto line100
    end
    @label line95
    if DABS(RS1) > 0.0
        @goto line300
    end
    if ZR < 0.0
        @goto line300
    end
    NZ = N
    for I = 1:N
        YR[I] = ZEROR
        YI[I] = ZEROI
        @label line96
    end
    return NZ
    @label line100
    S1R = CYR[1]
    S1I = CYI[1]
    S2R = CYR[2]
    S2I = CYI[2]
    C1R = CSRR[KFLAG]
    ASCLE = BRY[KFLAG]
    for I = IB:N
        C2R = S2R
        C2I = S2I
        S2R = (CKR * C2R - CKI * C2I) + S1R
        S2I = CKR * C2I + CKI * C2R + S1I
        S1R = C2R
        S1I = C2I
        CKR = CKR + RZR
        CKI = CKI + RZI
        C2R = S2R * C1R
        C2I = S2I * C1R
        YR[I] = C2R
        YI[I] = C2I
        if KFLAG >= 3
            @goto line120
        end
        STR = DABS(C2R)
        STI = DABS(C2I)
        C2M = DMAX1(STR,STI)
        if C2M <= ASCLE
            @goto line120
        end
        KFLAG = KFLAG + 1
        ASCLE = BRY[KFLAG]
        S1R = S1R * C1R
        S1I = S1I * C1R
        S2R = C2R
        S2I = C2I
        S1R = S1R * CSSR[KFLAG]
        S1I = S1I * CSSR[KFLAG]
        S2R = S2R * CSSR[KFLAG]
        S2I = S2I * CSSR[KFLAG]
        C1R = CSRR[KFLAG]
        @label line120
    end
    @label line160
    if MR == 0
        return NZ
    end
    NZ = 0
    FMR = DBLE(FLOAT(MR))
    SGN = -(DSIGN(PI,FMR))
    CSGNI = SGN
    INU = INT(SNGL(FNU))
    FNF = FNU - DBLE(FLOAT(INU))
    IFN = (INU + N) - 1
    ANG = FNF * SGN
    CSPNR = DCOS(ANG)
    CSPNI = DSIN(ANG)
    if MOD(IFN,2) == 0
        @goto line170
    end
    CSPNR = -CSPNR
    CSPNI = -CSPNI
    @label line170
    ASC = BRY[1]
    IUF = 0
    KK = N
    KDFLG = 1
    IB = IB - 1
    IC = IB - 1
    for K = 1:N
        FN = FNU + DBLE(FLOAT(KK - 1))
        M = 3
        if N > 2
            @goto line175
        end
        @label line172
        INITD = INIT[J]
        PHIDR = PHIR[J]
        PHIDI = PHII[J]
        ZET1DR = ZETA1R[J]
        ZET1DI = ZETA1I[J]
        ZET2DR = ZETA2R[J]
        ZET2DI = ZETA2I[J]
        SUMDR = SUMR[J]
        SUMDI = SUMI[J]
        M = J
        J = 3 - J
        @goto line180
        @label line175
        if KK == N && IB < N
            @goto line180
        end
        if KK == IB || KK == IC
            @goto line172
        end
        INITD = 0
        @label line180
        (INITD,PHIDR,PHIDI,ZET1DR,ZET1DI,ZET2DR,ZET2DI,SUMDR,SUMDI) = ZUNIK(ZRR,ZRI,FN,1,0,TOL,INITD,PHIDR,PHIDI,ZET1DR,ZET1DI,ZET2DR,ZET2DI,SUMDR,SUMDI,slice(CWRKR,1:16,int(M)),slice(CWRKI,1:16,int(M)))
        if KODE == 1
            @goto line200
        end
        STR = ZRR + ZET2DR
        STI = ZRI + ZET2DI
        RAST = FN / ZABS(COMPLEX(STR,STI))
        STR = STR * RAST * RAST
        STI = -STI * RAST * RAST
        S1R = -ZET1DR + STR
        S1I = -ZET1DI + STI
        @goto line210
        @label line200
        S1R = -ZET1DR + ZET2DR
        S1I = -ZET1DI + ZET2DI
        @label line210
        RS1 = S1R
        if DABS(RS1) > ELIM
            @goto line260
        end
        if KDFLG == 1
            IFLAG = 2
        end
        if DABS(RS1) < ALIM
            @goto line220
        end
        APHI = ZABS(COMPLEX(PHIDR,PHIDI))
        RS1 = RS1 + DLOG(APHI)
        if DABS(RS1) > ELIM
            @goto line260
        end
        if KDFLG == 1
            IFLAG = 1
        end
        if RS1 < 0.0
            @goto line220
        end
        if KDFLG == 1
            IFLAG = 3
        end
        @label line220
        STR = PHIDR * SUMDR - PHIDI * SUMDI
        STI = PHIDR * SUMDI + PHIDI * SUMDR
        S2R = -CSGNI * STI
        S2I = CSGNI * STR
        STR = DEXP(S1R) * CSSR[IFLAG]
        S1R = STR * DCOS(S1I)
        S1I = STR * DSIN(S1I)
        STR = S2R * S1R - S2I * S1I
        S2I = S2R * S1I + S2I * S1R
        S2R = STR
        if IFLAG != 1
            @goto line230
        end
        (NW,) = ZUCHK(S2R,S2I,NW,BRY[1],TOL)
        if NW == 0
            @goto line230
        end
        S2R = ZEROR
        S2I = ZEROI
        @label line230
        CYR[KDFLG] = S2R
        CYI[KDFLG] = S2I
        C2R = S2R
        C2I = S2I
        S2R = S2R * CSRR[IFLAG]
        S2I = S2I * CSRR[IFLAG]
        S1R = YR[KK]
        S1I = YI[KK]
        if KODE == 1
            @goto line250
        end
        (S1R,S1I,S2R,S2I,NW,IUF) = ZS1S2(ZRR,ZRI,S1R,S1I,S2R,S2I,NW,ASC,ALIM,IUF)
        NZ = NZ + NW
        @label line250
        YR[KK] = (S1R * CSPNR - S1I * CSPNI) + S2R
        YI[KK] = CSPNR * S1I + CSPNI * S1R + S2I
        KK = KK - 1
        CSPNR = -CSPNR
        CSPNI = -CSPNI
        if C2R != 0.0 || C2I != 0.0
            @goto line255
        end
        KDFLG = 1
        @goto line270
        @label line255
        if KDFLG == 2
            @goto line275
        end
        KDFLG = 2
        @goto line270
        @label line260
        if RS1 > 0.0
            @goto line300
        end
        S2R = ZEROR
        S2I = ZEROI
        @goto line230
        @label line270
    end
    K = N
    @label line275
    IL = N - K
    if IL == 0
        return NZ
    end
    S1R = CYR[1]
    S1I = CYI[1]
    S2R = CYR[2]
    S2I = CYI[2]
    CSR = CSRR[IFLAG]
    ASCLE = BRY[IFLAG]
    FN = DBLE(FLOAT(INU + IL))
    for I = 1:IL
        C2R = S2R
        C2I = S2I
        S2R = S1R + (FN + FNF) * (RZR * C2R - RZI * C2I)
        S2I = S1I + (FN + FNF) * (RZR * C2I + RZI * C2R)
        S1R = C2R
        S1I = C2I
        FN = FN - 1.0
        C2R = S2R * CSR
        C2I = S2I * CSR
        CKR = C2R
        CKI = C2I
        C1R = YR[KK]
        C1I = YI[KK]
        if KODE == 1
            @goto line280
        end
        (C1R,C1I,C2R,C2I,NW,IUF) = ZS1S2(ZRR,ZRI,C1R,C1I,C2R,C2I,NW,ASC,ALIM,IUF)
        NZ = NZ + NW
        @label line280
        YR[KK] = (C1R * CSPNR - C1I * CSPNI) + C2R
        YI[KK] = C1R * CSPNI + C1I * CSPNR + C2I
        KK = KK - 1
        CSPNR = -CSPNR
        CSPNI = -CSPNI
        if IFLAG >= 3
            @goto line290
        end
        C2R = DABS(CKR)
        C2I = DABS(CKI)
        C2M = DMAX1(C2R,C2I)
        if C2M <= ASCLE
            @goto line290
        end
        IFLAG = IFLAG + 1
        ASCLE = BRY[IFLAG]
        S1R = S1R * CSR
        S1I = S1I * CSR
        S2R = CKR
        S2I = CKI
        S1R = S1R * CSSR[IFLAG]
        S1I = S1I * CSSR[IFLAG]
        S2R = S2R * CSSR[IFLAG]
        S2I = S2I * CSSR[IFLAG]
        CSR = CSRR[IFLAG]
        @label line290
    end
    return NZ
    @label line300
    NZ = -1
    return NZ
end
