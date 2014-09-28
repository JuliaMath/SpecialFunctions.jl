function ZRATI(ZR::Float64,ZI::Float64,FNU::Float64,N::Integer,CYR::AbstractArray{Float64},CYI::AbstractArray{Float64},TOL::Float64)
    AK::Float64 = 0
    AMAGZ::Float64 = 0
    AP1::Float64 = 0
    AP2::Float64 = 0
    ARG::Float64 = 0
    AZ::Float64 = 0
    CDFNUI::Float64 = 0
    CDFNUR::Float64 = 0
    CONEI::Float64 = 0
    CONER::Float64 = 0
    CZEROI::Float64 = 0
    CZEROR::Float64 = 0
    DFNU::Float64 = 0
    FDNU::Float64 = 0
    FLAM::Float64 = 0
    FNUP::Float64 = 0
    I::Int32 = 0
    ID::Int32 = 0
    IDNU::Int32 = 0
    INU::Int32 = 0
    ITIME::Int32 = 0
    K::Int32 = 0
    KK::Int32 = 0
    MAGZ::Int32 = 0
    P1I::Float64 = 0
    P1R::Float64 = 0
    P2I::Float64 = 0
    P2R::Float64 = 0
    PTI::Float64 = 0
    PTR::Float64 = 0
    RAK::Float64 = 0
    RAP1::Float64 = 0
    RHO::Float64 = 0
    RT2::Float64 = 0
    RZI::Float64 = 0
    RZR::Float64 = 0
    T1I::Float64 = 0
    T1R::Float64 = 0
    TEST::Float64 = 0
    TEST1::Float64 = 0
    TTI::Float64 = 0
    TTR::Float64 = 0
    begin 
        CZEROR = 0.0
        CZEROI = 0.0
        CONER = 1.0
        CONEI = 0.0
        RT2 = 1.4142135623730951
    end
    AZ = ZABS(COMPLEX(ZR,ZI))
    INU = INT(SNGL(FNU))
    IDNU = (INU + N) - 1
    MAGZ = INT(SNGL(AZ))
    AMAGZ = DBLE(FLOAT(MAGZ + 1))
    FDNU = DBLE(FLOAT(IDNU))
    FNUP = DMAX1(AMAGZ,FDNU)
    ID = (IDNU - MAGZ) - 1
    ITIME = 1
    K = 1
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
    if ID > 0
        ID = 0
    end
    AP2 = ZABS(COMPLEX(P2R,P2I))
    AP1 = ZABS(COMPLEX(P1R,P1I))
    ARG = (AP2 + AP2) / (AP1 * TOL)
    TEST1 = DSQRT(ARG)
    TEST = TEST1
    RAP1 = 1.0 / AP1
    P1R = P1R * RAP1
    P1I = P1I * RAP1
    P2R = P2R * RAP1
    P2I = P2I * RAP1
    AP2 = AP2 * RAP1
    @label line10
    K = K + 1
    AP1 = AP2
    PTR = P2R
    PTI = P2I
    P2R = P1R - (T1R * PTR - T1I * PTI)
    P2I = P1I - (T1R * PTI + T1I * PTR)
    P1R = PTR
    P1I = PTI
    T1R = T1R + RZR
    T1I = T1I + RZI
    AP2 = ZABS(COMPLEX(P2R,P2I))
    if AP1 <= TEST
        @goto line10
    end
    if ITIME == 2
        @goto line20
    end
    AK = ZABS(COMPLEX(T1R,T1I) * 0.5)
    FLAM = AK + DSQRT(AK * AK - 1.0)
    RHO = DMIN1(AP2 / AP1,FLAM)
    TEST = TEST1 * DSQRT(RHO / (RHO * RHO - 1.0))
    ITIME = 2
    @goto line10
    @label line20
    KK = (K + 1) - ID
    AK = DBLE(FLOAT(KK))
    T1R = AK
    T1I = CZEROI
    DFNU = FNU + DBLE(FLOAT(N - 1))
    P1R = 1.0 / AP2
    P1I = CZEROI
    P2R = CZEROR
    P2I = CZEROI
    for I = 1:KK
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
    (CYR[N],CYI[N]) = ZDIV(P2R,P2I,P1R,P1I,CYR[N],CYI[N])
    if N == 1
        return
    end
    K = N - 1
    AK = DBLE(FLOAT(K))
    T1R = AK
    T1I = CZEROI
    CDFNUR = FNU * RZR
    CDFNUI = FNU * RZI
    for I = 2:N
        PTR = CDFNUR + (T1R * RZR - T1I * RZI) + CYR[K + 1]
        PTI = CDFNUI + (T1R * RZI + T1I * RZR) + CYI[K + 1]
        AK = ZABS(COMPLEX(PTR,PTI))
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
        K = K - 1
        @label line60
    end
    return
end
