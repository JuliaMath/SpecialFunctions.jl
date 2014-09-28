function ZUNK2(ZR::Float64,ZI::Float64,FNU::Float64,KODE::Integer,MR::Integer,N::Integer,YR::AbstractArray{Float64},YI::AbstractArray{Float64},NZ::Integer,TOL::Float64,ELIM::Float64,ALIM::Float64)
    AARG::Float64 = 0
    AIC::Float64 = 0
    AII::Float64 = 0
    AIR::Float64 = 0
    ANG::Float64 = 0
    APHI::Float64 = 0
    ARGDI::Float64 = 0
    ARGDR::Float64 = 0
    ARGI = Array(Float64,2)
    ARGR = Array(Float64,2)
    ASC::Float64 = 0
    ASCLE::Float64 = 0
    ASUMDI::Float64 = 0
    ASUMDR::Float64 = 0
    ASUMI = Array(Float64,2)
    ASUMR = Array(Float64,2)
    BRY = Array(Float64,3)
    BSUMDI::Float64 = 0
    BSUMDR::Float64 = 0
    BSUMI = Array(Float64,2)
    BSUMR = Array(Float64,2)
    C1I::Float64 = 0
    C1R::Float64 = 0
    C2I::Float64 = 0
    C2M::Float64 = 0
    C2R::Float64 = 0
    CAR::Float64 = 0
    CIPI = Array(Float64,4)
    CIPR = Array(Float64,4)
    CKI::Float64 = 0
    CKR::Float64 = 0
    CONER::Float64 = 0
    CR1I::Float64 = 0
    CR1R::Float64 = 0
    CR2I::Float64 = 0
    CR2R::Float64 = 0
    CRSC::Float64 = 0
    CSCL::Float64 = 0
    CSGNI::Float64 = 0
    CSI::Float64 = 0
    CSPNI::Float64 = 0
    CSPNR::Float64 = 0
    CSR::Float64 = 0
    CSRR = Array(Float64,3)
    CSSR = Array(Float64,3)
    CYI = Array(Float64,2)
    CYR = Array(Float64,2)
    DAII::Float64 = 0
    DAIR::Float64 = 0
    FMR::Float64 = 0
    FN::Float64 = 0
    FNF::Float64 = 0
    HPI::Float64 = 0
    I::Int32 = 0
    IB::Int32 = 0
    IC::Int32 = 0
    IDUM::Int32 = 0
    IFLAG::Int32 = 0
    IFN::Int32 = 0
    IL::Int32 = 0
    IN::Int32 = 0
    INU::Int32 = 0
    IPARD::Int32 = 0
    IUF::Int32 = 0
    J::Int32 = 0
    K::Int32 = 0
    KDFLG::Int32 = 0
    KFLAG::Int32 = 0
    KK::Int32 = 0
    NAI::Int32 = 0
    NDAI::Int32 = 0
    NW::Int32 = 0
    PHIDI::Float64 = 0
    PHIDR::Float64 = 0
    PHII = Array(Float64,2)
    PHIR = Array(Float64,2)
    PI::Float64 = 0
    PTI::Float64 = 0
    PTR::Float64 = 0
    RAST::Float64 = 0
    RAZR::Float64 = 0
    RS1::Float64 = 0
    RZI::Float64 = 0
    RZR::Float64 = 0
    S1I::Float64 = 0
    S1R::Float64 = 0
    S2I::Float64 = 0
    S2R::Float64 = 0
    SAR::Float64 = 0
    SGN::Float64 = 0
    STI::Float64 = 0
    STR::Float64 = 0
    YY::Float64 = 0
    ZBI::Float64 = 0
    ZBR::Float64 = 0
    ZEROI::Float64 = 0
    ZEROR::Float64 = 0
    ZET1DI::Float64 = 0
    ZET1DR::Float64 = 0
    ZET2DI::Float64 = 0
    ZET2DR::Float64 = 0
    ZETA1I = Array(Float64,2)
    ZETA1R = Array(Float64,2)
    ZETA2I = Array(Float64,2)
    ZETA2R = Array(Float64,2)
    ZNI::Float64 = 0
    ZNR::Float64 = 0
    ZRI::Float64 = 0
    ZRR::Float64 = 0
    begin 
        ZEROR = 0.0
        ZEROI = 0.0
        CONER = 1.0
        CR1R = 1.0
        CR1I = 1.7320508075688772
        CR2R = -0.5
        CR2I = -0.8660254037844386
    end
    begin 
        HPI = 1.5707963267948966
        PI = 3.141592653589793
        AIC = 1.2655121234846454
    end
    begin 
        CIPR[1] = 1.0
        CIPI[1] = 0.0
        CIPR[2] = 0.0
        CIPI[2] = -1.0
        CIPR[3] = -1.0
        CIPI[3] = 0.0
        CIPR[4] = 0.0
        CIPI[4] = 1.0
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
    YY = ZRI
    ZNR = ZRI
    ZNI = -ZRR
    ZBR = ZRR
    ZBI = ZRI
    INU = INT(SNGL(FNU))
    FNF = FNU - DBLE(FLOAT(INU))
    ANG = -HPI * FNF
    CAR = DCOS(ANG)
    SAR = DSIN(ANG)
    C2R = HPI * SAR
    C2I = -HPI * CAR
    KK = MOD(INU,4) + 1
    STR = C2R * CIPR[KK] - C2I * CIPI[KK]
    STI = C2R * CIPI[KK] + C2I * CIPR[KK]
    CSR = CR1R * STR - CR1I * STI
    CSI = CR1R * STI + CR1I * STR
    if YY > 0.0
        @goto line20
    end
    ZNR = -ZNR
    ZBI = -ZBI
    @label line20
    J = 2
    for I = 1:N
        J = 3 - J
        FN = FNU + DBLE(FLOAT(I - 1))
        (PHIR[J],PHII[J],ARGR[J],ARGI[J],ZETA1R[J],ZETA1I[J],ZETA2R[J],ZETA2I[J],ASUMR[J],ASUMI[J],BSUMR[J],BSUMI[J]) = ZUNHJ(ZNR,ZNI,FN,0,TOL,PHIR[J],PHII[J],ARGR[J],ARGI[J],ZETA1R[J],ZETA1I[J],ZETA2R[J],ZETA2I[J],ASUMR[J],ASUMI[J],BSUMR[J],BSUMI[J])
        if KODE == 1
            @goto line30
        end
        STR = ZBR + ZETA2R[J]
        STI = ZBI + ZETA2I[J]
        RAST = FN / ZABS(COMPLEX(STR,STI))
        STR = STR * RAST * RAST
        STI = -STI * RAST * RAST
        S1R = ZETA1R[J] - STR
        S1I = ZETA1I[J] - STI
        @goto line40
        @label line30
        S1R = ZETA1R[J] - ZETA2R[J]
        S1I = ZETA1I[J] - ZETA2I[J]
        @label line40
        RS1 = S1R
        if DABS(RS1) > ELIM
            @goto line70
        end
        if KDFLG == 1
            KFLAG = 2
        end
        if DABS(RS1) < ALIM
            @goto line50
        end
        APHI = ZABS(COMPLEX(PHIR[J],PHII[J]))
        AARG = ZABS(COMPLEX(ARGR[J],ARGI[J]))
        RS1 = ((RS1 + DLOG(APHI)) - 0.25 * DLOG(AARG)) - AIC
        if DABS(RS1) > ELIM
            @goto line70
        end
        if KDFLG == 1
            KFLAG = 1
        end
        if RS1 < 0.0
            @goto line50
        end
        if KDFLG == 1
            KFLAG = 3
        end
        @label line50
        C2R = ARGR[J] * CR2R - ARGI[J] * CR2I
        C2I = ARGR[J] * CR2I + ARGI[J] * CR2R
        (AIR,AII,NAI,IDUM) = ZAIRY(C2R,C2I,0,2,AIR,AII,NAI,IDUM)
        (DAIR,DAII,NDAI,IDUM) = ZAIRY(C2R,C2I,1,2,DAIR,DAII,NDAI,IDUM)
        STR = DAIR * BSUMR[J] - DAII * BSUMI[J]
        STI = DAIR * BSUMI[J] + DAII * BSUMR[J]
        PTR = STR * CR2R - STI * CR2I
        PTI = STR * CR2I + STI * CR2R
        STR = PTR + (AIR * ASUMR[J] - AII * ASUMI[J])
        STI = PTI + (AIR * ASUMI[J] + AII * ASUMR[J])
        PTR = STR * PHIR[J] - STI * PHII[J]
        PTI = STR * PHII[J] + STI * PHIR[J]
        S2R = PTR * CSR - PTI * CSI
        S2I = PTR * CSI + PTI * CSR
        STR = DEXP(S1R) * CSSR[KFLAG]
        S1R = STR * DCOS(S1I)
        S1I = STR * DSIN(S1I)
        STR = S2R * S1R - S2I * S1I
        S2I = S1R * S2I + S2R * S1I
        S2R = STR
        if KFLAG != 1
            @goto line60
        end
        (NW,) = ZUCHK(S2R,S2I,NW,BRY[1],TOL)
        if NW != 0
            @goto line70
        end
        @label line60
        if YY <= 0.0
            S2I = -S2I
        end
        CYR[KDFLG] = S2R
        CYI[KDFLG] = S2I
        YR[I] = S2R * CSRR[KFLAG]
        YI[I] = S2I * CSRR[KFLAG]
        STR = CSI
        CSI = -CSR
        CSR = STR
        if KDFLG == 2
            @goto line85
        end
        KDFLG = 2
        @goto line80
        @label line70
        if RS1 > 0.0
            @goto line320
        end
        if ZR < 0.0
            @goto line320
        end
        KDFLG = 1
        YR[I] = ZEROR
        YI[I] = ZEROI
        NZ = NZ + 1
        STR = CSI
        CSI = -CSR
        CSR = STR
        if I == 1
            @goto line80
        end
        if YR[I - 1] == ZEROR && YI[I - 1] == ZEROI
            @goto line80
        end
        YR[I - 1] = ZEROR
        YI[I - 1] = ZEROI
        NZ = NZ + 1
        @label line80
    end
    I = N
    @label line85
    RAZR = 1.0 / ZABS(COMPLEX(ZRR,ZRI))
    STR = ZRR * RAZR
    STI = -ZRI * RAZR
    RZR = (STR + STR) * RAZR
    RZI = (STI + STI) * RAZR
    CKR = FN * RZR
    CKI = FN * RZI
    IB = I + 1
    if N < IB
        @goto line180
    end
    FN = FNU + DBLE(FLOAT(N - 1))
    IPARD = 1
    if MR != 0
        IPARD = 0
    end
    (PHIDR,PHIDI,ARGDR,ARGDI,ZET1DR,ZET1DI,ZET2DR,ZET2DI,ASUMDR,ASUMDI,BSUMDR,BSUMDI) = ZUNHJ(ZNR,ZNI,FN,IPARD,TOL,PHIDR,PHIDI,ARGDR,ARGDI,ZET1DR,ZET1DI,ZET2DR,ZET2DI,ASUMDR,ASUMDI,BSUMDR,BSUMDI)
    if KODE == 1
        @goto line90
    end
    STR = ZBR + ZET2DR
    STI = ZBI + ZET2DI
    RAST = FN / ZABS(COMPLEX(STR,STI))
    STR = STR * RAST * RAST
    STI = -STI * RAST * RAST
    S1R = ZET1DR - STR
    S1I = ZET1DI - STI
    @goto line100
    @label line90
    S1R = ZET1DR - ZET2DR
    S1I = ZET1DI - ZET2DI
    @label line100
    RS1 = S1R
    if DABS(RS1) > ELIM
        @goto line105
    end
    if DABS(RS1) < ALIM
        @goto line120
    end
    APHI = ZABS(COMPLEX(PHIDR,PHIDI))
    RS1 = RS1 + DLOG(APHI)
    if DABS(RS1) < ELIM
        @goto line120
    end
    @label line105
    if RS1 > 0.0
        @goto line320
    end
    if ZR < 0.0
        @goto line320
    end
    NZ = N
    for I = 1:N
        YR[I] = ZEROR
        YI[I] = ZEROI
        @label line106
    end
    return NZ
    @label line120
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
            @goto line130
        end
        STR = DABS(C2R)
        STI = DABS(C2I)
        C2M = DMAX1(STR,STI)
        if C2M <= ASCLE
            @goto line130
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
        @label line130
    end
    @label line180
    if MR == 0
        return NZ
    end
    NZ = 0
    FMR = DBLE(FLOAT(MR))
    SGN = -(DSIGN(PI,FMR))
    CSGNI = SGN
    if YY <= 0.0
        CSGNI = -CSGNI
    end
    IFN = (INU + N) - 1
    ANG = FNF * SGN
    CSPNR = DCOS(ANG)
    CSPNI = DSIN(ANG)
    if MOD(IFN,2) == 0
        @goto line190
    end
    CSPNR = -CSPNR
    CSPNI = -CSPNI
    @label line190
    CSR = SAR * CSGNI
    CSI = CAR * CSGNI
    IN = MOD(IFN,4) + 1
    C2R = CIPR[IN]
    C2I = CIPI[IN]
    STR = CSR * C2R + CSI * C2I
    CSI = -CSR * C2I + CSI * C2R
    CSR = STR
    ASC = BRY[1]
    IUF = 0
    KK = N
    KDFLG = 1
    IB = IB - 1
    IC = IB - 1
    for K = 1:N
        FN = FNU + DBLE(FLOAT(KK - 1))
        if N > 2
            @goto line175
        end
        @label line172
        PHIDR = PHIR[J]
        PHIDI = PHII[J]
        ARGDR = ARGR[J]
        ARGDI = ARGI[J]
        ZET1DR = ZETA1R[J]
        ZET1DI = ZETA1I[J]
        ZET2DR = ZETA2R[J]
        ZET2DI = ZETA2I[J]
        ASUMDR = ASUMR[J]
        ASUMDI = ASUMI[J]
        BSUMDR = BSUMR[J]
        BSUMDI = BSUMI[J]
        J = 3 - J
        @goto line210
        @label line175
        if KK == N && IB < N
            @goto line210
        end
        if KK == IB || KK == IC
            @goto line172
        end
        (PHIDR,PHIDI,ARGDR,ARGDI,ZET1DR,ZET1DI,ZET2DR,ZET2DI,ASUMDR,ASUMDI,BSUMDR,BSUMDI) = ZUNHJ(ZNR,ZNI,FN,0,TOL,PHIDR,PHIDI,ARGDR,ARGDI,ZET1DR,ZET1DI,ZET2DR,ZET2DI,ASUMDR,ASUMDI,BSUMDR,BSUMDI)
        @label line210
        if KODE == 1
            @goto line220
        end
        STR = ZBR + ZET2DR
        STI = ZBI + ZET2DI
        RAST = FN / ZABS(COMPLEX(STR,STI))
        STR = STR * RAST * RAST
        STI = -STI * RAST * RAST
        S1R = -ZET1DR + STR
        S1I = -ZET1DI + STI
        @goto line230
        @label line220
        S1R = -ZET1DR + ZET2DR
        S1I = -ZET1DI + ZET2DI
        @label line230
        RS1 = S1R
        if DABS(RS1) > ELIM
            @goto line280
        end
        if KDFLG == 1
            IFLAG = 2
        end
        if DABS(RS1) < ALIM
            @goto line240
        end
        APHI = ZABS(COMPLEX(PHIDR,PHIDI))
        AARG = ZABS(COMPLEX(ARGDR,ARGDI))
        RS1 = ((RS1 + DLOG(APHI)) - 0.25 * DLOG(AARG)) - AIC
        if DABS(RS1) > ELIM
            @goto line280
        end
        if KDFLG == 1
            IFLAG = 1
        end
        if RS1 < 0.0
            @goto line240
        end
        if KDFLG == 1
            IFLAG = 3
        end
        @label line240
        (AIR,AII,NAI,IDUM) = ZAIRY(ARGDR,ARGDI,0,2,AIR,AII,NAI,IDUM)
        (DAIR,DAII,NDAI,IDUM) = ZAIRY(ARGDR,ARGDI,1,2,DAIR,DAII,NDAI,IDUM)
        STR = DAIR * BSUMDR - DAII * BSUMDI
        STI = DAIR * BSUMDI + DAII * BSUMDR
        STR = STR + (AIR * ASUMDR - AII * ASUMDI)
        STI = STI + (AIR * ASUMDI + AII * ASUMDR)
        PTR = STR * PHIDR - STI * PHIDI
        PTI = STR * PHIDI + STI * PHIDR
        S2R = PTR * CSR - PTI * CSI
        S2I = PTR * CSI + PTI * CSR
        STR = DEXP(S1R) * CSSR[IFLAG]
        S1R = STR * DCOS(S1I)
        S1I = STR * DSIN(S1I)
        STR = S2R * S1R - S2I * S1I
        S2I = S2R * S1I + S2I * S1R
        S2R = STR
        if IFLAG != 1
            @goto line250
        end
        (NW,) = ZUCHK(S2R,S2I,NW,BRY[1],TOL)
        if NW == 0
            @goto line250
        end
        S2R = ZEROR
        S2I = ZEROI
        @label line250
        if YY <= 0.0
            S2I = -S2I
        end
        CYR[KDFLG] = S2R
        CYI[KDFLG] = S2I
        C2R = S2R
        C2I = S2I
        S2R = S2R * CSRR[IFLAG]
        S2I = S2I * CSRR[IFLAG]
        S1R = YR[KK]
        S1I = YI[KK]
        if KODE == 1
            @goto line270
        end
        (S1R,S1I,S2R,S2I,NW,IUF) = ZS1S2(ZRR,ZRI,S1R,S1I,S2R,S2I,NW,ASC,ALIM,IUF)
        NZ = NZ + NW
        @label line270
        YR[KK] = (S1R * CSPNR - S1I * CSPNI) + S2R
        YI[KK] = S1R * CSPNI + S1I * CSPNR + S2I
        KK = KK - 1
        CSPNR = -CSPNR
        CSPNI = -CSPNI
        STR = CSI
        CSI = -CSR
        CSR = STR
        if C2R != 0.0 || C2I != 0.0
            @goto line255
        end
        KDFLG = 1
        @goto line290
        @label line255
        if KDFLG == 2
            @goto line295
        end
        KDFLG = 2
        @goto line290
        @label line280
        if RS1 > 0.0
            @goto line320
        end
        S2R = ZEROR
        S2I = ZEROI
        @goto line250
        @label line290
    end
    K = N
    @label line295
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
            @goto line300
        end
        (C1R,C1I,C2R,C2I,NW,IUF) = ZS1S2(ZRR,ZRI,C1R,C1I,C2R,C2I,NW,ASC,ALIM,IUF)
        NZ = NZ + NW
        @label line300
        YR[KK] = (C1R * CSPNR - C1I * CSPNI) + C2R
        YI[KK] = C1R * CSPNI + C1I * CSPNR + C2I
        KK = KK - 1
        CSPNR = -CSPNR
        CSPNI = -CSPNI
        if IFLAG >= 3
            @goto line310
        end
        C2R = DABS(CKR)
        C2I = DABS(CKI)
        C2M = DMAX1(C2R,C2I)
        if C2M <= ASCLE
            @goto line310
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
        @label line310
    end
    return NZ
    @label line320
    NZ = -1
    return NZ
end
