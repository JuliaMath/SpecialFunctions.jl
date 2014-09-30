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
function ZUNK1(ZR::Float64,ZI::Float64,FNU::Float64,KODE::Int32,MR::Int32,N::Int32,YR::AbstractArray{Float64},YI::AbstractArray{Float64},NZ::Int32,TOL::Float64,ELIM::Float64,ALIM::Float64)
    ANG::Float64 = zero(Float64)
    APHI::Float64 = zero(Float64)
    ASC::Float64 = zero(Float64)
    ASCLE::Float64 = zero(Float64)
    const BRY = _ZUNK1_BRY
    C1I::Float64 = zero(Float64)
    C1R::Float64 = zero(Float64)
    C2I::Float64 = zero(Float64)
    C2M::Float64 = zero(Float64)
    C2R::Float64 = zero(Float64)
    CKI::Float64 = zero(Float64)
    CKR::Float64 = zero(Float64)
    CONER::Float64 = zero(Float64)
    CRSC::Float64 = zero(Float64)
    CSCL::Float64 = zero(Float64)
    CSGNI::Float64 = zero(Float64)
    CSPNI::Float64 = zero(Float64)
    CSPNR::Float64 = zero(Float64)
    CSR::Float64 = zero(Float64)
    const CSRR = _ZUNK1_CSRR
    const CSSR = _ZUNK1_CSSR
    const CWRKI = _ZUNK1_CWRKI
    const CWRKR = _ZUNK1_CWRKR
    const CYI = _ZUNK1_CYI
    const CYR = _ZUNK1_CYR
    FMR::Float64 = zero(Float64)
    FN::Float64 = zero(Float64)
    FNF::Float64 = zero(Float64)
    I::Int32 = zero(Int32)
    IB::Int32 = zero(Int32)
    IC::Int32 = zero(Int32)
    IFLAG::Int32 = zero(Int32)
    IFN::Int32 = zero(Int32)
    IL::Int32 = zero(Int32)
    const INIT = _ZUNK1_INIT
    INITD::Int32 = zero(Int32)
    INU::Int32 = zero(Int32)
    IPARD::Int32 = zero(Int32)
    IUF::Int32 = zero(Int32)
    J::Int32 = zero(Int32)
    K::Int32 = zero(Int32)
    KDFLG::Int32 = zero(Int32)
    KFLAG::Int32 = zero(Int32)
    KK::Int32 = zero(Int32)
    NW::Int32 = zero(Int32)
    PHIDI::Float64 = zero(Float64)
    PHIDR::Float64 = zero(Float64)
    const PHII = _ZUNK1_PHII
    const PHIR = _ZUNK1_PHIR
    PI::Float64 = zero(Float64)
    RAST::Float64 = zero(Float64)
    RAZR::Float64 = zero(Float64)
    RS1::Float64 = zero(Float64)
    RZI::Float64 = zero(Float64)
    RZR::Float64 = zero(Float64)
    S1I::Float64 = zero(Float64)
    S1R::Float64 = zero(Float64)
    S2I::Float64 = zero(Float64)
    S2R::Float64 = zero(Float64)
    SGN::Float64 = zero(Float64)
    STI::Float64 = zero(Float64)
    STR::Float64 = zero(Float64)
    SUMDI::Float64 = zero(Float64)
    SUMDR::Float64 = zero(Float64)
    const SUMI = _ZUNK1_SUMI
    const SUMR = _ZUNK1_SUMR
    ZEROI::Float64 = zero(Float64)
    ZEROR::Float64 = zero(Float64)
    ZET1DI::Float64 = zero(Float64)
    ZET1DR::Float64 = zero(Float64)
    ZET2DI::Float64 = zero(Float64)
    ZET2DR::Float64 = zero(Float64)
    const ZETA1I = _ZUNK1_ZETA1I
    const ZETA1R = _ZUNK1_ZETA1R
    const ZETA2I = _ZUNK1_ZETA2I
    const ZETA2R = _ZUNK1_ZETA2R
    ZRI::Float64 = zero(Float64)
    ZRR::Float64 = zero(Float64)
    begin 
        ZEROR = 0.0
        ZEROI = 0.0
        CONER = 1.0
    end
    begin 
        PI = 3.141592653589793
    end
    KDFLG = int32(1)
    NZ = int32(0)
    CSCL = 1.0 / TOL
    CRSC = TOL
    CSSR[int32(1)] = CSCL
    CSSR[int32(2)] = CONER
    CSSR[int32(3)] = CRSC
    CSRR[int32(1)] = CRSC
    CSRR[int32(2)] = CONER
    CSRR[int32(3)] = CSCL
    BRY[int32(1)] = (1000.0D1MACH1) / TOL
    BRY[int32(2)] = 1.0 / BRY[int32(1)]
    BRY[int32(3)] = D1MACH2
    ZRR = ZR
    ZRI = ZI
    if ZR >= 0.0
        @goto line10
    end
    ZRR = -ZR
    ZRI = -ZI
    @label line10
    J = int32(2)
    for I = int32(1):N
        J = int32(3) - J
        FN = FNU + DBLE(FLOAT(I - int32(1)))
        INIT[J] = int32(0)
        (INIT[J],PHIR[J],PHII[J],ZETA1R[J],ZETA1I[J],ZETA2R[J],ZETA2I[J],SUMR[J],SUMI[J]) = ZUNIK(ZRR,ZRI,FN,int32(2),int32(0),TOL,INIT[J],PHIR[J],PHII[J],ZETA1R[J],ZETA1I[J],ZETA2R[J],ZETA2I[J],SUMR[J],SUMI[J],slice(CWRKR,int32(1):16,int(J)),slice(CWRKI,int32(1):16,int(J)))
        if KODE == int32(1)
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
        if KDFLG == int32(1)
            KFLAG = int32(2)
        end
        if DABS(RS1) < ALIM
            @goto line40
        end
        APHI = ZABS(COMPLEX(PHIR[J],PHII[J]))
        RS1 = RS1 + DLOG(APHI)
        if DABS(RS1) > ELIM
            @goto line60
        end
        if KDFLG == int32(1)
            KFLAG = int32(1)
        end
        if RS1 < 0.0
            @goto line40
        end
        if KDFLG == int32(1)
            KFLAG = int32(3)
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
        if KFLAG != int32(1)
            @goto line50
        end
        (NW,) = ZUCHK(S2R,S2I,NW,BRY[int32(1)],TOL)
        if NW != int32(0)
            @goto line60
        end
        @label line50
        CYR[KDFLG] = S2R
        CYI[KDFLG] = S2I
        YR[I] = S2R * CSRR[KFLAG]
        YI[I] = S2I * CSRR[KFLAG]
        if KDFLG == int32(2)
            @goto line75
        end
        KDFLG = int32(2)
        @goto line70
        @label line60
        if RS1 > 0.0
            @goto line300
        end
        if ZR < 0.0
            @goto line300
        end
        KDFLG = int32(1)
        YR[I] = ZEROR
        YI[I] = ZEROI
        NZ = NZ + int32(1)
        if I == int32(1)
            @goto line70
        end
        if YR[I - int32(1)] == ZEROR && YI[I - int32(1)] == ZEROI
            @goto line70
        end
        YR[I - int32(1)] = ZEROR
        YI[I - int32(1)] = ZEROI
        NZ = NZ + int32(1)
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
    IB = I + int32(1)
    if N < IB
        @goto line160
    end
    FN = FNU + DBLE(FLOAT(N - int32(1)))
    IPARD = int32(1)
    if MR != int32(0)
        IPARD = int32(0)
    end
    INITD = int32(0)
    (INITD,PHIDR,PHIDI,ZET1DR,ZET1DI,ZET2DR,ZET2DI,SUMDR,SUMDI) = ZUNIK(ZRR,ZRI,FN,int32(2),IPARD,TOL,INITD,PHIDR,PHIDI,ZET1DR,ZET1DI,ZET2DR,ZET2DI,SUMDR,SUMDI,slice(CWRKR,int32(1):16,int(int32(3))),slice(CWRKI,int32(1):16,int(int32(3))))
    if KODE == int32(1)
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
    for I = int32(1):N
        YR[I] = ZEROR
        YI[I] = ZEROI
        @label line96
    end
    return NZ
    @label line100
    S1R = CYR[int32(1)]
    S1I = CYI[int32(1)]
    S2R = CYR[int32(2)]
    S2I = CYI[int32(2)]
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
        if KFLAG >= int32(3)
            @goto line120
        end
        STR = DABS(C2R)
        STI = DABS(C2I)
        C2M = DMAX1(STR,STI)
        if C2M <= ASCLE
            @goto line120
        end
        KFLAG = KFLAG + int32(1)
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
    if MR == int32(0)
        return NZ
    end
    NZ = int32(0)
    FMR = DBLE(FLOAT(MR))
    SGN = -(DSIGN(PI,FMR))
    CSGNI = SGN
    INU = INT(SNGL(FNU))
    FNF = FNU - DBLE(FLOAT(INU))
    IFN = (INU + N) - int32(1)
    ANG = FNF * SGN
    CSPNR = DCOS(ANG)
    CSPNI = DSIN(ANG)
    if MOD(IFN,int32(2)) == int32(0)
        @goto line170
    end
    CSPNR = -CSPNR
    CSPNI = -CSPNI
    @label line170
    ASC = BRY[int32(1)]
    IUF = int32(0)
    KK = N
    KDFLG = int32(1)
    IB = IB - int32(1)
    IC = IB - int32(1)
    for K = int32(1):N
        FN = FNU + DBLE(FLOAT(KK - int32(1)))
        M = int32(3)
        if N > int32(2)
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
        J = int32(3) - J
        @goto line180
        @label line175
        if KK == N && IB < N
            @goto line180
        end
        if KK == IB || KK == IC
            @goto line172
        end
        INITD = int32(0)
        @label line180
        (INITD,PHIDR,PHIDI,ZET1DR,ZET1DI,ZET2DR,ZET2DI,SUMDR,SUMDI) = ZUNIK(ZRR,ZRI,FN,int32(1),int32(0),TOL,INITD,PHIDR,PHIDI,ZET1DR,ZET1DI,ZET2DR,ZET2DI,SUMDR,SUMDI,slice(CWRKR,int32(1):16,int(M)),slice(CWRKI,int32(1):16,int(M)))
        if KODE == int32(1)
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
        if KDFLG == int32(1)
            IFLAG = int32(2)
        end
        if DABS(RS1) < ALIM
            @goto line220
        end
        APHI = ZABS(COMPLEX(PHIDR,PHIDI))
        RS1 = RS1 + DLOG(APHI)
        if DABS(RS1) > ELIM
            @goto line260
        end
        if KDFLG == int32(1)
            IFLAG = int32(1)
        end
        if RS1 < 0.0
            @goto line220
        end
        if KDFLG == int32(1)
            IFLAG = int32(3)
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
        if IFLAG != int32(1)
            @goto line230
        end
        (NW,) = ZUCHK(S2R,S2I,NW,BRY[int32(1)],TOL)
        if NW == int32(0)
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
        if KODE == int32(1)
            @goto line250
        end
        (S1R,S1I,S2R,S2I,NW,IUF) = ZS1S2(ZRR,ZRI,S1R,S1I,S2R,S2I,NW,ASC,ALIM,IUF)
        NZ = NZ + NW
        @label line250
        YR[KK] = (S1R * CSPNR - S1I * CSPNI) + S2R
        YI[KK] = CSPNR * S1I + CSPNI * S1R + S2I
        KK = KK - int32(1)
        CSPNR = -CSPNR
        CSPNI = -CSPNI
        if C2R != 0.0 || C2I != 0.0
            @goto line255
        end
        KDFLG = int32(1)
        @goto line270
        @label line255
        if KDFLG == int32(2)
            @goto line275
        end
        KDFLG = int32(2)
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
    if IL == int32(0)
        return NZ
    end
    S1R = CYR[int32(1)]
    S1I = CYI[int32(1)]
    S2R = CYR[int32(2)]
    S2I = CYI[int32(2)]
    CSR = CSRR[IFLAG]
    ASCLE = BRY[IFLAG]
    FN = DBLE(FLOAT(INU + IL))
    for I = int32(1):IL
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
        if KODE == int32(1)
            @goto line280
        end
        (C1R,C1I,C2R,C2I,NW,IUF) = ZS1S2(ZRR,ZRI,C1R,C1I,C2R,C2I,NW,ASC,ALIM,IUF)
        NZ = NZ + NW
        @label line280
        YR[KK] = (C1R * CSPNR - C1I * CSPNI) + C2R
        YI[KK] = C1R * CSPNI + C1I * CSPNR + C2I
        KK = KK - int32(1)
        CSPNR = -CSPNR
        CSPNI = -CSPNI
        if IFLAG >= int32(3)
            @goto line290
        end
        C2R = DABS(CKR)
        C2I = DABS(CKI)
        C2M = DMAX1(C2R,C2I)
        if C2M <= ASCLE
            @goto line290
        end
        IFLAG = IFLAG + int32(1)
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
    NZ = int32(-1)
    return NZ
end
