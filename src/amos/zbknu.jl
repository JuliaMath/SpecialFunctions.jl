function ZBKNU(ZR::Float64,ZI::Float64,FNU::Float64,KODE::Integer,N::Integer,YR::AbstractArray{Float64},YI::AbstractArray{Float64},NZ::Integer,TOL::Float64,ELIM::Float64,ALIM::Float64)
    A1::Float64 = 0
    A2::Float64 = 0
    AA::Float64 = 0
    AK::Float64 = 0
    ALAS::Float64 = 0
    AS::Float64 = 0
    ASCLE::Float64 = 0
    BB::Float64 = 0
    BK::Float64 = 0
    BRY = Array(Float64,3)
    CAZ::Float64 = 0
    CBI::Float64 = 0
    CBR::Float64 = 0
    CC = Array(Float64,8)
    CCHI::Float64 = 0
    CCHR::Float64 = 0
    CELMR::Float64 = 0
    CKI::Float64 = 0
    CKR::Float64 = 0
    COEFI::Float64 = 0
    COEFR::Float64 = 0
    CONEI::Float64 = 0
    CONER::Float64 = 0
    CRSCR::Float64 = 0
    CSCLR::Float64 = 0
    CSHI::Float64 = 0
    CSHR::Float64 = 0
    CSI::Float64 = 0
    CSR::Float64 = 0
    CSRR = Array(Float64,3)
    CSSR = Array(Float64,3)
    CTWOR::Float64 = 0
    CYI = Array(Float64,2)
    CYR = Array(Float64,2)
    CZEROI::Float64 = 0
    CZEROR::Float64 = 0
    CZI::Float64 = 0
    CZR::Float64 = 0
    DNU::Float64 = 0
    DNU2::Float64 = 0
    DPI::Float64 = 0
    ELM::Float64 = 0
    ETEST::Float64 = 0
    FC::Float64 = 0
    FHS::Float64 = 0
    FI::Float64 = 0
    FK::Float64 = 0
    FKS::Float64 = 0
    FMUI::Float64 = 0
    FMUR::Float64 = 0
    FPI::Float64 = 0
    FR::Float64 = 0
    G1::Float64 = 0
    G2::Float64 = 0
    HELIM::Float64 = 0
    HPI::Float64 = 0
    I::Int32 = 0
    IC::Int32 = 0
    IDUM::Int32 = 0
    IFLAG::Int32 = 0
    INU::Int32 = 0
    INUB::Int32 = 0
    J::Int32 = 0
    K::Int32 = 0
    KFLAG::Int32 = 0
    KK::Int32 = 0
    KMAX::Int32 = 0
    KODED::Int32 = 0
    NW::Int32 = 0
    P1I::Float64 = 0
    P1R::Float64 = 0
    P2I::Float64 = 0
    P2M::Float64 = 0
    P2R::Float64 = 0
    PI::Float64 = 0
    PR::Float64 = 0
    PTI::Float64 = 0
    PTR::Float64 = 0
    QI::Float64 = 0
    QR::Float64 = 0
    R1::Float64 = 0
    RAK::Float64 = 0
    RCAZ::Float64 = 0
    RTHPI::Float64 = 0
    RZI::Float64 = 0
    RZR::Float64 = 0
    S::Float64 = 0
    S1I::Float64 = 0
    S1R::Float64 = 0
    S2I::Float64 = 0
    S2R::Float64 = 0
    SMUI::Float64 = 0
    SMUR::Float64 = 0
    SPI::Float64 = 0
    STI::Float64 = 0
    STR::Float64 = 0
    T1::Float64 = 0
    T2::Float64 = 0
    TM::Float64 = 0
    TTH::Float64 = 0
    ZDI::Float64 = 0
    ZDR::Float64 = 0
    begin 
        KMAX = 30
    end
    begin 
        CZEROR = 0.0
        CZEROI = 0.0
        CONER = 1.0
        CONEI = 0.0
        CTWOR = 2.0
        R1 = 2.0
    end
    begin 
        DPI = 3.141592653589793
        RTHPI = 1.2533141373155003
        SPI = 1.909859317102744
        HPI = 1.5707963267948966
        FPI = 1.8976999933151775
        TTH = 0.6666666666666666
    end
    begin 
        CC[1] = 0.5772156649015329
        CC[2] = -0.04200263503409524
        CC[3] = -0.04219773455554433
        CC[4] = 0.0072189432466631
        CC[5] = -0.00021524167411495098
        CC[6] = -2.013485478078824e-5
        CC[7] = 1.133027231981696e-6
        CC[8] = 6.116095104481416e-9
    end
    CAZ = ZABS(COMPLEX(ZR,ZI))
    CSCLR = 1.0 / TOL
    CRSCR = TOL
    CSSR[1] = CSCLR
    CSSR[2] = 1.0
    CSSR[3] = CRSCR
    CSRR[1] = CRSCR
    CSRR[2] = 1.0
    CSRR[3] = CSCLR
    BRY[1] = (1000.0D1MACH1) / TOL
    BRY[2] = 1.0 / BRY[1]
    BRY[3] = D1MACH2
    NZ = 0
    IFLAG = 0
    KODED = KODE
    RCAZ = 1.0 / CAZ
    STR = ZR * RCAZ
    STI = -ZI * RCAZ
    RZR = (STR + STR) * RCAZ
    RZI = (STI + STI) * RCAZ
    INU = INT(SNGL(FNU + 0.5))
    DNU = FNU - DBLE(FLOAT(INU))
    if DABS(DNU) == 0.5
        @goto line110
    end
    DNU2 = 0.0
    if DABS(DNU) > TOL
        DNU2 = DNU * DNU
    end
    if CAZ > R1
        @goto line110
    end
    FC = 1.0
    (SMUR,SMUI,IDUM) = ZLOG(RZR,RZI,SMUR,SMUI,IDUM)
    FMUR = SMUR * DNU
    FMUI = SMUI * DNU
    (CSHR,CSHI,CCHR,CCHI) = ZSHCH(FMUR,FMUI,CSHR,CSHI,CCHR,CCHI)
    if DNU == 0.0
        @goto line10
    end
    FC = DNU * DPI
    FC = FC / DSIN(FC)
    SMUR = CSHR / DNU
    SMUI = CSHI / DNU
    @label line10
    A2 = 1.0 + DNU
    T2 = DEXP(-(DGAMLN(A2,IDUM)))
    T1 = 1.0 / (T2 * FC)
    if DABS(DNU) > 0.1
        @goto line40
    end
    AK = 1.0
    S = CC[1]
    for K = 2:8
        AK = AK * DNU2
        TM = CC[K] * AK
        S = S + TM
        if DABS(TM) < TOL
            @goto line30
        end
        @label line20
    end
    @label line30
    G1 = -S
    @goto line50
    @label line40
    G1 = (T1 - T2) / (DNU + DNU)
    @label line50
    G2 = (T1 + T2) * 0.5
    FR = FC * (CCHR * G1 + SMUR * G2)
    FI = FC * (CCHI * G1 + SMUI * G2)
    (STR,STI) = ZEXP(FMUR,FMUI,STR,STI)
    PR = (0.5STR) / T2
    PI = (0.5STI) / T2
    (PTR,PTI) = ZDIV(0.5,0.0,STR,STI,PTR,PTI)
    QR = PTR / T1
    QI = PTI / T1
    S1R = FR
    S1I = FI
    S2R = PR
    S2I = PI
    AK = 1.0
    A1 = 1.0
    CKR = CONER
    CKI = CONEI
    BK = 1.0 - DNU2
    if INU > 0 || N > 1
        @goto line80
    end
    if CAZ < TOL
        @goto line70
    end
    (CZR,CZI) = ZMLT(ZR,ZI,ZR,ZI,CZR,CZI)
    CZR = 0.25CZR
    CZI = 0.25CZI
    T1 = 0.25 * CAZ * CAZ
    @label line60
    FR = (FR * AK + PR + QR) / BK
    FI = (FI * AK + PI + QI) / BK
    STR = 1.0 / (AK - DNU)
    PR = PR * STR
    PI = PI * STR
    STR = 1.0 / (AK + DNU)
    QR = QR * STR
    QI = QI * STR
    STR = CKR * CZR - CKI * CZI
    RAK = 1.0 / AK
    CKI = (CKR * CZI + CKI * CZR) * RAK
    CKR = STR * RAK
    S1R = (CKR * FR - CKI * FI) + S1R
    S1I = CKR * FI + CKI * FR + S1I
    A1 = A1 * T1 * RAK
    BK = BK + AK + AK + 1.0
    AK = AK + 1.0
    if A1 > TOL
        @goto line60
    end
    @label line70
    YR[1] = S1R
    YI[1] = S1I
    if KODED == 1
        return NZ
    end
    (STR,STI) = ZEXP(ZR,ZI,STR,STI)
    (YR[1],YI[1]) = ZMLT(S1R,S1I,STR,STI,YR[1],YI[1])
    return NZ
    @label line80
    if CAZ < TOL
        @goto line100
    end
    (CZR,CZI) = ZMLT(ZR,ZI,ZR,ZI,CZR,CZI)
    CZR = 0.25CZR
    CZI = 0.25CZI
    T1 = 0.25 * CAZ * CAZ
    @label line90
    FR = (FR * AK + PR + QR) / BK
    FI = (FI * AK + PI + QI) / BK
    STR = 1.0 / (AK - DNU)
    PR = PR * STR
    PI = PI * STR
    STR = 1.0 / (AK + DNU)
    QR = QR * STR
    QI = QI * STR
    STR = CKR * CZR - CKI * CZI
    RAK = 1.0 / AK
    CKI = (CKR * CZI + CKI * CZR) * RAK
    CKR = STR * RAK
    S1R = (CKR * FR - CKI * FI) + S1R
    S1I = CKR * FI + CKI * FR + S1I
    STR = PR - FR * AK
    STI = PI - FI * AK
    S2R = (CKR * STR - CKI * STI) + S2R
    S2I = CKR * STI + CKI * STR + S2I
    A1 = A1 * T1 * RAK
    BK = BK + AK + AK + 1.0
    AK = AK + 1.0
    if A1 > TOL
        @goto line90
    end
    @label line100
    KFLAG = 2
    A1 = FNU + 1.0
    AK = A1 * DABS(SMUR)
    if AK > ALIM
        KFLAG = 3
    end
    STR = CSSR[KFLAG]
    P2R = S2R * STR
    P2I = S2I * STR
    (S2R,S2I) = ZMLT(P2R,P2I,RZR,RZI,S2R,S2I)
    S1R = S1R * STR
    S1I = S1I * STR
    if KODED == 1
        @goto line210
    end
    (FR,FI) = ZEXP(ZR,ZI,FR,FI)
    (S1R,S1I) = ZMLT(S1R,S1I,FR,FI,S1R,S1I)
    (S2R,S2I) = ZMLT(S2R,S2I,FR,FI,S2R,S2I)
    @goto line210
    @label line110
    (STR,STI) = ZSQRT(ZR,ZI,STR,STI)
    (COEFR,COEFI) = ZDIV(RTHPI,CZEROI,STR,STI,COEFR,COEFI)
    KFLAG = 2
    if KODED == 2
        @goto line120
    end
    if ZR > ALIM
        @goto line290
    end
    STR = DEXP(-ZR) * CSSR[KFLAG]
    STI = -STR * DSIN(ZI)
    STR = STR * DCOS(ZI)
    (COEFR,COEFI) = ZMLT(COEFR,COEFI,STR,STI,COEFR,COEFI)
    @label line120
    if DABS(DNU) == 0.5
        @goto line300
    end
    AK = DCOS(DPI * DNU)
    AK = DABS(AK)
    if AK == CZEROR
        @goto line300
    end
    FHS = DABS(0.25 - DNU2)
    if FHS == CZEROR
        @goto line300
    end
    T1 = DBLE(FLOAT(I1MACH14 - 1))
    T1 = T1 * D1MACH5 * 3.321928094
    T1 = DMAX1(T1,12.0)
    T1 = DMIN1(T1,60.0)
    T2 = TTH * T1 - 6.0
    if ZR != 0.0
        @goto line130
    end
    T1 = HPI
    @goto line140
    @label line130
    T1 = DATAN(ZI / ZR)
    T1 = DABS(T1)
    @label line140
    if T2 > CAZ
        @goto line170
    end
    ETEST = AK / (DPI * CAZ * TOL)
    FK = CONER
    if ETEST < CONER
        @goto line180
    end
    FKS = CTWOR
    CKR = CAZ + CAZ + CTWOR
    P1R = CZEROR
    P2R = CONER
    for I = 1:KMAX
        AK = FHS / FKS
        CBR = CKR / (FK + CONER)
        PTR = P2R
        P2R = CBR * P2R - P1R * AK
        P1R = PTR
        CKR = CKR + CTWOR
        FKS = FKS + FK + FK + CTWOR
        FHS = FHS + FK + FK
        FK = FK + CONER
        STR = DABS(P2R) * FK
        if ETEST < STR
            @goto line160
        end
        @label line150
    end
    @goto line310
    @label line160
    FK = FK + SPI * T1 * DSQRT(T2 / CAZ)
    FHS = DABS(0.25 - DNU2)
    @goto line180
    @label line170
    A2 = DSQRT(CAZ)
    AK = (FPI * AK) / (TOL * DSQRT(A2))
    AA = (3.0T1) / (1.0 + CAZ)
    BB = (14.7T1) / (28.0 + CAZ)
    AK = (DLOG(AK) + (CAZ * DCOS(AA)) / (1.0 + 0.008CAZ)) / DCOS(BB)
    FK = (0.12125 * AK * AK) / CAZ + 1.5
    @label line180
    K = INT(SNGL(FK))
    FK = DBLE(FLOAT(K))
    FKS = FK * FK
    P1R = CZEROR
    P1I = CZEROI
    P2R = TOL
    P2I = CZEROI
    CSR = P2R
    CSI = P2I
    for I = 1:K
        A1 = FKS - FK
        AK = (FKS + FK) / (A1 + FHS)
        RAK = 2.0 / (FK + CONER)
        CBR = (FK + ZR) * RAK
        CBI = ZI * RAK
        PTR = P2R
        PTI = P2I
        P2R = ((PTR * CBR - PTI * CBI) - P1R) * AK
        P2I = ((PTI * CBR + PTR * CBI) - P1I) * AK
        P1R = PTR
        P1I = PTI
        CSR = CSR + P2R
        CSI = CSI + P2I
        FKS = (A1 - FK) + CONER
        FK = FK - CONER
        @label line190
    end
    TM = ZABS(COMPLEX(CSR,CSI))
    PTR = 1.0 / TM
    S1R = P2R * PTR
    S1I = P2I * PTR
    CSR = CSR * PTR
    CSI = -CSI * PTR
    (STR,STI) = ZMLT(COEFR,COEFI,S1R,S1I,STR,STI)
    (S1R,S1I) = ZMLT(STR,STI,CSR,CSI,S1R,S1I)
    if INU > 0 || N > 1
        @goto line200
    end
    ZDR = ZR
    ZDI = ZI
    if IFLAG == 1
        @goto line270
    end
    @goto line240
    @label line200
    TM = ZABS(COMPLEX(P2R,P2I))
    PTR = 1.0 / TM
    P1R = P1R * PTR
    P1I = P1I * PTR
    P2R = P2R * PTR
    P2I = -P2I * PTR
    (PTR,PTI) = ZMLT(P1R,P1I,P2R,P2I,PTR,PTI)
    STR = (DNU + 0.5) - PTR
    STI = -PTI
    (STR,STI) = ZDIV(STR,STI,ZR,ZI,STR,STI)
    STR = STR + 1.0
    (S2R,S2I) = ZMLT(STR,STI,S1R,S1I,S2R,S2I)
    @label line210
    STR = DNU + 1.0
    CKR = STR * RZR
    CKI = STR * RZI
    if N == 1
        INU = INU - 1
    end
    if INU > 0
        @goto line220
    end
    if N > 1
        @goto line215
    end
    S1R = S2R
    S1I = S2I
    @label line215
    ZDR = ZR
    ZDI = ZI
    if IFLAG == 1
        @goto line270
    end
    @goto line240
    @label line220
    INUB = 1
    if IFLAG == 1
        @goto line261
    end
    @label line225
    P1R = CSRR[KFLAG]
    ASCLE = BRY[KFLAG]
    for I = INUB:INU
        STR = S2R
        STI = S2I
        S2R = (CKR * STR - CKI * STI) + S1R
        S2I = CKR * STI + CKI * STR + S1I
        S1R = STR
        S1I = STI
        CKR = CKR + RZR
        CKI = CKI + RZI
        if KFLAG >= 3
            @goto line230
        end
        P2R = S2R * P1R
        P2I = S2I * P1R
        STR = DABS(P2R)
        STI = DABS(P2I)
        P2M = DMAX1(STR,STI)
        if P2M <= ASCLE
            @goto line230
        end
        KFLAG = KFLAG + 1
        ASCLE = BRY[KFLAG]
        S1R = S1R * P1R
        S1I = S1I * P1R
        S2R = P2R
        S2I = P2I
        STR = CSSR[KFLAG]
        S1R = S1R * STR
        S1I = S1I * STR
        S2R = S2R * STR
        S2I = S2I * STR
        P1R = CSRR[KFLAG]
        @label line230
    end
    if N != 1
        @goto line240
    end
    S1R = S2R
    S1I = S2I
    @label line240
    STR = CSRR[KFLAG]
    YR[1] = S1R * STR
    YI[1] = S1I * STR
    if N == 1
        return NZ
    end
    YR[2] = S2R * STR
    YI[2] = S2I * STR
    if N == 2
        return NZ
    end
    KK = 2
    @label line250
    KK = KK + 1
    if KK > N
        return NZ
    end
    P1R = CSRR[KFLAG]
    ASCLE = BRY[KFLAG]
    for I = KK:N
        P2R = S2R
        P2I = S2I
        S2R = (CKR * P2R - CKI * P2I) + S1R
        S2I = CKI * P2R + CKR * P2I + S1I
        S1R = P2R
        S1I = P2I
        CKR = CKR + RZR
        CKI = CKI + RZI
        P2R = S2R * P1R
        P2I = S2I * P1R
        YR[I] = P2R
        YI[I] = P2I
        if KFLAG >= 3
            @goto line260
        end
        STR = DABS(P2R)
        STI = DABS(P2I)
        P2M = DMAX1(STR,STI)
        if P2M <= ASCLE
            @goto line260
        end
        KFLAG = KFLAG + 1
        ASCLE = BRY[KFLAG]
        S1R = S1R * P1R
        S1I = S1I * P1R
        S2R = P2R
        S2I = P2I
        STR = CSSR[KFLAG]
        S1R = S1R * STR
        S1I = S1I * STR
        S2R = S2R * STR
        S2I = S2I * STR
        P1R = CSRR[KFLAG]
        @label line260
    end
    return NZ
    @label line261
    HELIM = 0.5ELIM
    ELM = DEXP(-ELIM)
    CELMR = ELM
    ASCLE = BRY[1]
    ZDR = ZR
    ZDI = ZI
    IC = -1
    J = 2
    for I = 1:INU
        STR = S2R
        STI = S2I
        S2R = (STR * CKR - STI * CKI) + S1R
        S2I = STI * CKR + STR * CKI + S1I
        S1R = STR
        S1I = STI
        CKR = CKR + RZR
        CKI = CKI + RZI
        AS = ZABS(COMPLEX(S2R,S2I))
        ALAS = DLOG(AS)
        P2R = -ZDR + ALAS
        if P2R < -ELIM
            @goto line263
        end
        (STR,STI,IDUM) = ZLOG(S2R,S2I,STR,STI,IDUM)
        P2R = -ZDR + STR
        P2I = -ZDI + STI
        P2M = DEXP(P2R) / TOL
        P1R = P2M * DCOS(P2I)
        P1I = P2M * DSIN(P2I)
        (NW,) = ZUCHK(P1R,P1I,NW,ASCLE,TOL)
        if NW != 0
            @goto line263
        end
        J = 3 - J
        CYR[J] = P1R
        CYI[J] = P1I
        if IC == I - 1
            @goto line264
        end
        IC = I
        @goto line262
        @label line263
        if ALAS < HELIM
            @goto line262
        end
        ZDR = ZDR - ELIM
        S1R = S1R * CELMR
        S1I = S1I * CELMR
        S2R = S2R * CELMR
        S2I = S2I * CELMR
        @label line262
    end
    if N != 1
        @goto line270
    end
    S1R = S2R
    S1I = S2I
    @goto line270
    @label line264
    KFLAG = 1
    INUB = I + 1
    S2R = CYR[J]
    S2I = CYI[J]
    J = 3 - J
    S1R = CYR[J]
    S1I = CYI[J]
    if INUB <= INU
        @goto line225
    end
    if N != 1
        @goto line240
    end
    S1R = S2R
    S1I = S2I
    @goto line240
    @label line270
    YR[1] = S1R
    YI[1] = S1I
    if N == 1
        @goto line280
    end
    YR[2] = S2R
    YI[2] = S2I
    @label line280
    ASCLE = BRY[1]
    (NZ,) = ZKSCL(ZDR,ZDI,FNU,N,YR,YI,NZ,RZR,RZI,ASCLE,TOL,ELIM)
    INU = N - NZ
    if INU <= 0
        return NZ
    end
    KK = NZ + 1
    S1R = YR[KK]
    S1I = YI[KK]
    YR[KK] = S1R * CSRR[1]
    YI[KK] = S1I * CSRR[1]
    if INU == 1
        return NZ
    end
    KK = NZ + 2
    S2R = YR[KK]
    S2I = YI[KK]
    YR[KK] = S2R * CSRR[1]
    YI[KK] = S2I * CSRR[1]
    if INU == 2
        return NZ
    end
    T2 = FNU + DBLE(FLOAT(KK - 1))
    CKR = T2 * RZR
    CKI = T2 * RZI
    KFLAG = 1
    @goto line250
    @label line290
    KODED = 2
    IFLAG = 1
    KFLAG = 2
    @goto line120
    @label line300
    S1R = COEFR
    S1I = COEFI
    S2R = COEFR
    S2I = COEFI
    @goto line210
    @label line310
    NZ = -2
    return NZ
end
