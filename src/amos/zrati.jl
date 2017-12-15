function ZRATI(ZR::Float64,ZI::Float64,FNU::Float64,N::Int32,CYR::AbstractArray{Float64},CYI::AbstractArray{Float64},TOL::Float64)
    AK::Float64 = zero(Float64)
    AMAGZ::Float64 = zero(Float64)
    AP1::Float64 = zero(Float64)
    AP2::Float64 = zero(Float64)
    ARG::Float64 = zero(Float64)
    AZ::Float64 = zero(Float64)
    CDFNUI::Float64 = zero(Float64)
    CDFNUR::Float64 = zero(Float64)
    CONEI::Float64 = zero(Float64)
    CONER::Float64 = zero(Float64)
    CZEROI::Float64 = zero(Float64)
    CZEROR::Float64 = zero(Float64)
    DFNU::Float64 = zero(Float64)
    FDNU::Float64 = zero(Float64)
    FLAM::Float64 = zero(Float64)
    FNUP::Float64 = zero(Float64)
    I::Int32 = zero(Int32)
    ID::Int32 = zero(Int32)
    IDNU::Int32 = zero(Int32)
    INU::Int32 = zero(Int32)
    ITIME::Int32 = zero(Int32)
    K::Int32 = zero(Int32)
    KK::Int32 = zero(Int32)
    MAGZ::Int32 = zero(Int32)
    P1I::Float64 = zero(Float64)
    P1R::Float64 = zero(Float64)
    P2I::Float64 = zero(Float64)
    P2R::Float64 = zero(Float64)
    PTI::Float64 = zero(Float64)
    PTR::Float64 = zero(Float64)
    RAK::Float64 = zero(Float64)
    RAP1::Float64 = zero(Float64)
    RHO::Float64 = zero(Float64)
    RT2::Float64 = zero(Float64)
    RZI::Float64 = zero(Float64)
    RZR::Float64 = zero(Float64)
    T1I::Float64 = zero(Float64)
    T1R::Float64 = zero(Float64)
    TEST::Float64 = zero(Float64)
    TEST1::Float64 = zero(Float64)
    TTI::Float64 = zero(Float64)
    TTR::Float64 = zero(Float64)
    begin 
        CZEROR = 0.0
        CZEROI = 0.0
        CONER = 1.0
        CONEI = 0.0
        RT2 = 1.4142135623730951
    end
    AZ = abs(complex(ZR,ZI))
    INU = trunc(Int32,SNGL(FNU))
    IDNU = (INU + N) - Int32(1)
    MAGZ = trunc(Int32,SNGL(AZ))
    AMAGZ = DBLE(FLOAT(MAGZ + Int32(1)))
    FDNU = DBLE(FLOAT(IDNU))
    FNUP = max(AMAGZ,FDNU)
    ID = (IDNU - MAGZ) - Int32(1)
    ITIME = Int32(1)
    K = Int32(1)
    PTR = 1.0 / AZ
    RZR = PTR * (ZR + ZR) * PTR
    RZI = -PTR * (ZI + ZI) * PTR
    T1R = RZR * FNUP
    T1I = RZI * FNUP
    P2R = -T1R
    P2I = -T1I
    P1R = CONER
    P1I = CONEI
    T1R = T1R + RZR
    T1I = T1I + RZI
    if ID > Int32(0)
        ID = Int32(0)
    end
    AP2 = abs(complex(P2R,P2I))
    AP1 = abs(complex(P1R,P1I))
    ARG = (AP2 + AP2) / (AP1 * TOL)
    TEST1 = sqrt(ARG)
    TEST = TEST1
    RAP1 = 1.0 / AP1
    P1R = P1R * RAP1
    P1I = P1I * RAP1
    P2R = P2R * RAP1
    P2I = P2I * RAP1
    AP2 = AP2 * RAP1
    @label line10
    K = K + Int32(1)
    AP1 = AP2
    PTR = P2R
    PTI = P2I
    P2R = P1R - (T1R * PTR - T1I * PTI)
    P2I = P1I - (T1R * PTI + T1I * PTR)
    P1R = PTR
    P1I = PTI
    T1R = T1R + RZR
    T1I = T1I + RZI
    AP2 = abs(complex(P2R,P2I))
    if AP1 <= TEST
        @goto line10
    end
    if ITIME == Int32(2)
        @goto line20
    end
    AK = abs(complex(T1R,T1I) * 0.5)
    FLAM = AK + sqrt(AK * AK - 1.0)
    RHO = min(AP2 / AP1,FLAM)
    TEST = TEST1 * sqrt(RHO / (RHO * RHO - 1.0))
    ITIME = Int32(2)
    @goto line10
    @label line20
    KK = (K + Int32(1)) - ID
    AK = DBLE(FLOAT(KK))
    T1R = AK
    T1I = CZEROI
    DFNU = FNU + DBLE(FLOAT(N - Int32(1)))
    P1R = 1.0 / AP2
    P1I = CZEROI
    P2R = CZEROR
    P2I = CZEROI
    for I = Int32(1):KK
        PTR = P1R
        PTI = P1I
        RAP1 = DFNU + T1R
        TTR = RZR * RAP1
        TTI = RZI * RAP1
        P1R = (PTR * TTR - PTI * TTI) + P2R
        P1I = (PTR * TTI + PTI * TTR) + P2I
        P2R = PTR
        P2I = PTI
        T1R = T1R - CONER
        @label line30
    end
    if P1R != CZEROR || P1I != CZEROI
        @goto line40
    end
    P1R = TOL
    P1I = TOL
    @label line40
    (CYR[N],CYI[N]) = reim(complex(P2R,P2I) / complex(P1R,P1I))
    if N == Int32(1)
        return
    end
    K = N - Int32(1)
    AK = DBLE(FLOAT(K))
    T1R = AK
    T1I = CZEROI
    CDFNUR = FNU * RZR
    CDFNUI = FNU * RZI
    for I = Int32(2):N
        PTR = CDFNUR + (T1R * RZR - T1I * RZI) + CYR[K + Int32(1)]
        PTI = CDFNUI + (T1R * RZI + T1I * RZR) + CYI[K + Int32(1)]
        AK = abs(complex(PTR,PTI))
        if AK != CZEROR
            @goto line50
        end
        PTR = TOL
        PTI = TOL
        AK = TOL * RT2
        @label line50
        RAK = CONER / AK
        CYR[K] = RAK * PTR * RAK
        CYI[K] = -RAK * PTI * RAK
        T1R = T1R - CONER
        K = K - Int32(1)
        @label line60
    end
    return
end
