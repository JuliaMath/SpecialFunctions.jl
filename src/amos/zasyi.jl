function ZASYI(ZR::Float64,ZI::Float64,FNU::Float64,KODE::Integer,N::Integer,YR::AbstractArray{Float64},YI::AbstractArray{Float64},NZ::Integer,RL::Float64,TOL::Float64,ELIM::Float64,ALIM::Float64)
    AA::Float64 = 0
    AEZ::Float64 = 0
    AK::Float64 = 0
    AK1I::Float64 = 0
    AK1R::Float64 = 0
    ARG::Float64 = 0
    ARM::Float64 = 0
    ATOL::Float64 = 0
    AZ::Float64 = 0
    BB::Float64 = 0
    BK::Float64 = 0
    CKI::Float64 = 0
    CKR::Float64 = 0
    CONEI::Float64 = 0
    CONER::Float64 = 0
    CS1I::Float64 = 0
    CS1R::Float64 = 0
    CS2I::Float64 = 0
    CS2R::Float64 = 0
    CZI::Float64 = 0
    CZR::Float64 = 0
    DFNU::Float64 = 0
    DKI::Float64 = 0
    DKR::Float64 = 0
    DNU2::Float64 = 0
    EZI::Float64 = 0
    EZR::Float64 = 0
    FDN::Float64 = 0
    I::Int32 = 0
    IB::Int32 = 0
    IL::Int32 = 0
    INU::Int32 = 0
    J::Int32 = 0
    JL::Int32 = 0
    K::Int32 = 0
    KODED::Int32 = 0
    M::Int32 = 0
    NN::Int32 = 0
    P1I::Float64 = 0
    P1R::Float64 = 0
    PI::Float64 = 0
    RAZ::Float64 = 0
    RTPI::Float64 = 0
    RTR1::Float64 = 0
    RZI::Float64 = 0
    RZR::Float64 = 0
    S::Float64 = 0
    S2I::Float64 = 0
    S2R::Float64 = 0
    SGN::Float64 = 0
    SQK::Float64 = 0
    STI::Float64 = 0
    STR::Float64 = 0
    TZI::Float64 = 0
    TZR::Float64 = 0
    ZEROI::Float64 = 0
    ZEROR::Float64 = 0
    begin 
        PI = 3.141592653589793
        RTPI = 0.15915494309189535
    end
    begin 
        ZEROR = 0.0
        ZEROI = 0.0
        CONER = 1.0
        CONEI = 0.0
    end
    NZ = 0
    AZ = ZABS(COMPLEX(ZR,ZI))
    ARM = 1000.0D1MACH1
    RTR1 = DSQRT(ARM)
    IL = MIN0(2,N)
    DFNU = FNU + DBLE(FLOAT(N - IL))
    RAZ = 1.0 / AZ
    STR = ZR * RAZ
    STI = -ZI * RAZ
    AK1R = RTPI * STR * RAZ
    AK1I = RTPI * STI * RAZ
    (AK1R,AK1I) = ZSQRT(AK1R,AK1I,AK1R,AK1I)
    CZR = ZR
    CZI = ZI
    if KODE != 2
        @goto line10
    end
    CZR = ZEROR
    CZI = ZI
    @label line10
    if DABS(CZR) > ELIM
        @goto line100
    end
    DNU2 = DFNU + DFNU
    KODED = 1
    if DABS(CZR) > ALIM && N > 2
        @goto line20
    end
    KODED = 0
    (STR,STI) = ZEXP(CZR,CZI,STR,STI)
    (AK1R,AK1I) = ZMLT(AK1R,AK1I,STR,STI,AK1R,AK1I)
    @label line20
    FDN = 0.0
    if DNU2 > RTR1
        FDN = DNU2 * DNU2
    end
    EZR = ZR * 8.0
    EZI = ZI * 8.0
    AEZ = 8.0AZ
    S = TOL / AEZ
    JL = INT(SNGL(RL + RL)) + 2
    P1R = ZEROR
    P1I = ZEROI
    if ZI == 0.0
        @goto line30
    end
    INU = INT(SNGL(FNU))
    ARG = (FNU - DBLE(FLOAT(INU))) * PI
    INU = (INU + N) - IL
    AK = -(DSIN(ARG))
    BK = DCOS(ARG)
    if ZI < 0.0
        BK = -BK
    end
    P1R = AK
    P1I = BK
    if MOD(INU,2) == 0
        @goto line30
    end
    P1R = -P1R
    P1I = -P1I
    @label line30
    for K = 1:IL
        SQK = FDN - 1.0
        ATOL = S * DABS(SQK)
        SGN = 1.0
        CS1R = CONER
        CS1I = CONEI
        CS2R = CONER
        CS2I = CONEI
        CKR = CONER
        CKI = CONEI
        AK = 0.0
        AA = 1.0
        BB = AEZ
        DKR = EZR
        DKI = EZI
        for J = 1:JL
            (STR,STI) = ZDIV(CKR,CKI,DKR,DKI,STR,STI)
            CKR = STR * SQK
            CKI = STI * SQK
            CS2R = CS2R + CKR
            CS2I = CS2I + CKI
            SGN = -SGN
            CS1R = CS1R + CKR * SGN
            CS1I = CS1I + CKI * SGN
            DKR = DKR + EZR
            DKI = DKI + EZI
            AA = (AA * DABS(SQK)) / BB
            BB = BB + AEZ
            AK = AK + 8.0
            SQK = SQK - AK
            if AA <= ATOL
                @goto line50
            end
            @label line40
        end
        @goto line110
        @label line50
        S2R = CS1R
        S2I = CS1I
        if ZR + ZR >= ELIM
            @goto line60
        end
        TZR = ZR + ZR
        TZI = ZI + ZI
        (STR,STI) = ZEXP(-TZR,-TZI,STR,STI)
        (STR,STI) = ZMLT(STR,STI,P1R,P1I,STR,STI)
        (STR,STI) = ZMLT(STR,STI,CS2R,CS2I,STR,STI)
        S2R = S2R + STR
        S2I = S2I + STI
        @label line60
        FDN = FDN + 8.0DFNU + 4.0
        P1R = -P1R
        P1I = -P1I
        M = (N - IL) + K
        YR[M] = S2R * AK1R - S2I * AK1I
        YI[M] = S2R * AK1I + S2I * AK1R
        @label line70
    end
    if N <= 2
        return NZ
    end
    NN = N
    K = NN - 2
    AK = DBLE(FLOAT(K))
    STR = ZR * RAZ
    STI = -ZI * RAZ
    RZR = (STR + STR) * RAZ
    RZI = (STI + STI) * RAZ
    IB = 3
    for I = IB:NN
        YR[K] = (AK + FNU) * (RZR * YR[K + 1] - RZI * YI[K + 1]) + YR[K + 2]
        YI[K] = (AK + FNU) * (RZR * YI[K + 1] + RZI * YR[K + 1]) + YI[K + 2]
        AK = AK - 1.0
        K = K - 1
        @label line80
    end
    if KODED == 0
        return NZ
    end
    (CKR,CKI) = ZEXP(CZR,CZI,CKR,CKI)
    for I = 1:NN
        STR = YR[I] * CKR - YI[I] * CKI
        YI[I] = YR[I] * CKI + YI[I] * CKR
        YR[I] = STR
        @label line90
    end
    return NZ
    @label line100
    NZ = -1
    return NZ
    @label line110
    NZ = -2
    return NZ
end
