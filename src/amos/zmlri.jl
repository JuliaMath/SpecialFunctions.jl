function ZMLRI(ZR::Float64,ZI::Float64,FNU::Float64,KODE::Integer,N::Integer,YR::AbstractArray{Float64},YI::AbstractArray{Float64},NZ::Integer,TOL::Float64)
    ACK::Float64 = 0
    AK::Float64 = 0
    AP::Float64 = 0
    AT::Float64 = 0
    AZ::Float64 = 0
    BK::Float64 = 0
    CKI::Float64 = 0
    CKR::Float64 = 0
    CNORMI::Float64 = 0
    CNORMR::Float64 = 0
    CONEI::Float64 = 0
    CONER::Float64 = 0
    FKAP::Float64 = 0
    FKK::Float64 = 0
    FLAM::Float64 = 0
    FNF::Float64 = 0
    I::Int32 = 0
    IAZ::Int32 = 0
    IDUM::Int32 = 0
    IFNU::Int32 = 0
    INU::Int32 = 0
    ITIME::Int32 = 0
    K::Int32 = 0
    KK::Int32 = 0
    KM::Int32 = 0
    M::Int32 = 0
    P1I::Float64 = 0
    P1R::Float64 = 0
    P2I::Float64 = 0
    P2R::Float64 = 0
    PTI::Float64 = 0
    PTR::Float64 = 0
    RAZ::Float64 = 0
    RHO::Float64 = 0
    RHO2::Float64 = 0
    RZI::Float64 = 0
    RZR::Float64 = 0
    SCLE::Float64 = 0
    STI::Float64 = 0
    STR::Float64 = 0
    SUMI::Float64 = 0
    SUMR::Float64 = 0
    TFNF::Float64 = 0
    TST::Float64 = 0
    ZEROI::Float64 = 0
    ZEROR::Float64 = 0
    begin 
        ZEROR = 0.0
        ZEROI = 0.0
        CONER = 1.0
        CONEI = 0.0
    end
    SCLE = D1MACH1 / TOL
    NZ = 0
    AZ = ZABS(COMPLEX(ZR,ZI))
    IAZ = INT(SNGL(AZ))
    IFNU = INT(SNGL(FNU))
    INU = (IFNU + N) - 1
    AT = DBLE(FLOAT(IAZ)) + 1.0
    RAZ = 1.0 / AZ
    STR = ZR * RAZ
    STI = -ZI * RAZ
    CKR = STR * AT * RAZ
    CKI = STI * AT * RAZ
    RZR = (STR + STR) * RAZ
    RZI = (STI + STI) * RAZ
    P1R = ZEROR
    P1I = ZEROI
    P2R = CONER
    P2I = CONEI
    ACK = (AT + 1.0) * RAZ
    RHO = ACK + DSQRT(ACK * ACK - 1.0)
    RHO2 = RHO * RHO
    TST = (RHO2 + RHO2) / ((RHO2 - 1.0) * (RHO - 1.0))
    TST = TST / TOL
    AK = AT
    for I = 1:80
        PTR = P2R
        PTI = P2I
        P2R = P1R - (CKR * PTR - CKI * PTI)
        P2I = P1I - (CKI * PTR + CKR * PTI)
        P1R = PTR
        P1I = PTI
        CKR = CKR + RZR
        CKI = CKI + RZI
        AP = ZABS(COMPLEX(P2R,P2I))
        if AP > TST * AK * AK
            @goto line20
        end
        AK = AK + 1.0
        @label line10
    end
    @goto line110
    @label line20
    I = I + 1
    K = 0
    if INU < IAZ
        @goto line40
    end
    P1R = ZEROR
    P1I = ZEROI
    P2R = CONER
    P2I = CONEI
    AT = DBLE(FLOAT(INU)) + 1.0
    STR = ZR * RAZ
    STI = -ZI * RAZ
    CKR = STR * AT * RAZ
    CKI = STI * AT * RAZ
    ACK = AT * RAZ
    TST = DSQRT(ACK / TOL)
    ITIME = 1
    for K = 1:80
        PTR = P2R
        PTI = P2I
        P2R = P1R - (CKR * PTR - CKI * PTI)
        P2I = P1I - (CKR * PTI + CKI * PTR)
        P1R = PTR
        P1I = PTI
        CKR = CKR + RZR
        CKI = CKI + RZI
        AP = ZABS(COMPLEX(P2R,P2I))
        if AP < TST
            @goto line30
        end
        if ITIME == 2
            @goto line40
        end
        ACK = ZABS(COMPLEX(CKR,CKI))
        FLAM = ACK + DSQRT(ACK * ACK - 1.0)
        FKAP = AP / ZABS(COMPLEX(P1R,P1I))
        RHO = DMIN1(FLAM,FKAP)
        TST = TST * DSQRT(RHO / (RHO * RHO - 1.0))
        ITIME = 2
        @label line30
    end
    @goto line110
    @label line40
    K = K + 1
    KK = MAX0(I + IAZ,K + INU)
    FKK = DBLE(FLOAT(KK))
    P1R = ZEROR
    P1I = ZEROI
    P2R = SCLE
    P2I = ZEROI
    FNF = FNU - DBLE(FLOAT(IFNU))
    TFNF = FNF + FNF
    BK = (DGAMLN(FKK + TFNF + 1.0,IDUM) - DGAMLN(FKK + 1.0,IDUM)) - DGAMLN(TFNF + 1.0,IDUM)
    BK = DEXP(BK)
    SUMR = ZEROR
    SUMI = ZEROI
    KM = KK - INU
    for I = 1:KM
        PTR = P2R
        PTI = P2I
        P2R = P1R + (FKK + FNF) * (RZR * PTR - RZI * PTI)
        P2I = P1I + (FKK + FNF) * (RZI * PTR + RZR * PTI)
        P1R = PTR
        P1I = PTI
        AK = 1.0 - TFNF / (FKK + TFNF)
        ACK = BK * AK
        SUMR = SUMR + (ACK + BK) * P1R
        SUMI = SUMI + (ACK + BK) * P1I
        BK = ACK
        FKK = FKK - 1.0
        @label line50
    end
    YR[N] = P2R
    YI[N] = P2I
    if N == 1
        @goto line70
    end
    for I = 2:N
        PTR = P2R
        PTI = P2I
        P2R = P1R + (FKK + FNF) * (RZR * PTR - RZI * PTI)
        P2I = P1I + (FKK + FNF) * (RZI * PTR + RZR * PTI)
        P1R = PTR
        P1I = PTI
        AK = 1.0 - TFNF / (FKK + TFNF)
        ACK = BK * AK
        SUMR = SUMR + (ACK + BK) * P1R
        SUMI = SUMI + (ACK + BK) * P1I
        BK = ACK
        FKK = FKK - 1.0
        M = (N - I) + 1
        YR[M] = P2R
        YI[M] = P2I
        @label line60
    end
    @label line70
    if IFNU <= 0
        @goto line90
    end
    for I = 1:IFNU
        PTR = P2R
        PTI = P2I
        P2R = P1R + (FKK + FNF) * (RZR * PTR - RZI * PTI)
        P2I = P1I + (FKK + FNF) * (RZR * PTI + RZI * PTR)
        P1R = PTR
        P1I = PTI
        AK = 1.0 - TFNF / (FKK + TFNF)
        ACK = BK * AK
        SUMR = SUMR + (ACK + BK) * P1R
        SUMI = SUMI + (ACK + BK) * P1I
        BK = ACK
        FKK = FKK - 1.0
        @label line80
    end
    @label line90
    PTR = ZR
    PTI = ZI
    if KODE == 2
        PTR = ZEROR
    end
    (STR,STI,IDUM) = ZLOG(RZR,RZI,STR,STI,IDUM)
    P1R = -FNF * STR + PTR
    P1I = -FNF * STI + PTI
    AP = DGAMLN(1.0 + FNF,IDUM)
    PTR = P1R - AP
    PTI = P1I
    P2R = P2R + SUMR
    P2I = P2I + SUMI
    AP = ZABS(COMPLEX(P2R,P2I))
    P1R = 1.0 / AP
    (STR,STI) = ZEXP(PTR,PTI,STR,STI)
    CKR = STR * P1R
    CKI = STI * P1R
    PTR = P2R * P1R
    PTI = -P2I * P1R
    (CNORMR,CNORMI) = ZMLT(CKR,CKI,PTR,PTI,CNORMR,CNORMI)
    for I = 1:N
        STR = YR[I] * CNORMR - YI[I] * CNORMI
        YI[I] = YR[I] * CNORMI + YI[I] * CNORMR
        YR[I] = STR
        @label line100
    end
    return NZ
    @label line110
    NZ = -2
    return NZ
end
