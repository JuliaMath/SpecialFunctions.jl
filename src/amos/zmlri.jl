function ZMLRI(ZR::Float64,ZI::Float64,FNU::Float64,KODE::Int32,N::Int32,YR::AbstractArray{Float64},YI::AbstractArray{Float64},NZ::Int32,TOL::Float64)
    ACK::Float64 = zero(Float64)
    AK::Float64 = zero(Float64)
    AP::Float64 = zero(Float64)
    AT::Float64 = zero(Float64)
    AZ::Float64 = zero(Float64)
    BK::Float64 = zero(Float64)
    CKI::Float64 = zero(Float64)
    CKR::Float64 = zero(Float64)
    CNORMI::Float64 = zero(Float64)
    CNORMR::Float64 = zero(Float64)
    CONEI::Float64 = zero(Float64)
    CONER::Float64 = zero(Float64)
    FKAP::Float64 = zero(Float64)
    FKK::Float64 = zero(Float64)
    FLAM::Float64 = zero(Float64)
    FNF::Float64 = zero(Float64)
    I::Int32 = zero(Int32)
    IAZ::Int32 = zero(Int32)
    IDUM::Int32 = zero(Int32)
    IFNU::Int32 = zero(Int32)
    INU::Int32 = zero(Int32)
    ITIME::Int32 = zero(Int32)
    K::Int32 = zero(Int32)
    KK::Int32 = zero(Int32)
    KM::Int32 = zero(Int32)
    M::Int32 = zero(Int32)
    P1I::Float64 = zero(Float64)
    P1R::Float64 = zero(Float64)
    P2I::Float64 = zero(Float64)
    P2R::Float64 = zero(Float64)
    PTI::Float64 = zero(Float64)
    PTR::Float64 = zero(Float64)
    RAZ::Float64 = zero(Float64)
    RHO::Float64 = zero(Float64)
    RHO2::Float64 = zero(Float64)
    RZI::Float64 = zero(Float64)
    RZR::Float64 = zero(Float64)
    SCLE::Float64 = zero(Float64)
    STI::Float64 = zero(Float64)
    STR::Float64 = zero(Float64)
    SUMI::Float64 = zero(Float64)
    SUMR::Float64 = zero(Float64)
    TFNF::Float64 = zero(Float64)
    TST::Float64 = zero(Float64)
    ZEROI::Float64 = zero(Float64)
    ZEROR::Float64 = zero(Float64)
    begin 
        ZEROR = 0.0
        ZEROI = 0.0
        CONER = 1.0
        CONEI = 0.0
    end
    SCLE = D1MACH1 / TOL
    NZ = Int32(0)
    AZ = abs(COMPLEX(ZR,ZI))
    IAZ = INT(SNGL(AZ))
    IFNU = INT(SNGL(FNU))
    INU = (IFNU + N) - Int32(1)
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
    for I = Int32(1):Int32(80)
        PTR = P2R
        PTI = P2I
        P2R = P1R - (CKR * PTR - CKI * PTI)
        P2I = P1I - (CKI * PTR + CKR * PTI)
        P1R = PTR
        P1I = PTI
        CKR = CKR + RZR
        CKI = CKI + RZI
        AP = abs(COMPLEX(P2R,P2I))
        if AP > TST * AK * AK
            @goto line20
        end
        AK = AK + 1.0
        @label line10
    end
    @goto line110
    @label line20
    I = I + Int32(1)
    K = Int32(0)
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
    ITIME = Int32(1)
    for K = Int32(1):Int32(80)
        PTR = P2R
        PTI = P2I
        P2R = P1R - (CKR * PTR - CKI * PTI)
        P2I = P1I - (CKR * PTI + CKI * PTR)
        P1R = PTR
        P1I = PTI
        CKR = CKR + RZR
        CKI = CKI + RZI
        AP = abs(COMPLEX(P2R,P2I))
        if AP < TST
            @goto line30
        end
        if ITIME == Int32(2)
            @goto line40
        end
        ACK = abs(COMPLEX(CKR,CKI))
        FLAM = ACK + DSQRT(ACK * ACK - 1.0)
        FKAP = AP / abs(COMPLEX(P1R,P1I))
        RHO = DMIN1(FLAM,FKAP)
        TST = TST * DSQRT(RHO / (RHO * RHO - 1.0))
        ITIME = Int32(2)
        @label line30
    end
    @goto line110
    @label line40
    K = K + Int32(1)
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
    for I = Int32(1):KM
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
    if N == Int32(1)
        @goto line70
    end
    for I = Int32(2):N
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
        M = (N - I) + Int32(1)
        YR[M] = P2R
        YI[M] = P2I
        @label line60
    end
    @label line70
    if IFNU <= Int32(0)
        @goto line90
    end
    for I = Int32(1):IFNU
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
    if KODE == Int32(2)
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
    AP = abs(COMPLEX(P2R,P2I))
    P1R = 1.0 / AP
    (STR,STI) = ZEXP(PTR,PTI,STR,STI)
    CKR = STR * P1R
    CKI = STI * P1R
    PTR = P2R * P1R
    PTI = -P2I * P1R
    (CNORMR,CNORMI) = ZMLT(CKR,CKI,PTR,PTI,CNORMR,CNORMI)
    for I = Int32(1):N
        STR = YR[I] * CNORMR - YI[I] * CNORMI
        YI[I] = YR[I] * CNORMI + YI[I] * CNORMR
        YR[I] = STR
        @label line100
    end
    return NZ
    @label line110
    NZ = Int32(-2)
    return NZ
end
