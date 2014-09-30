function ZASYI(ZR::Float64,ZI::Float64,FNU::Float64,KODE::Int32,N::Int32,YR::AbstractArray{Float64},YI::AbstractArray{Float64},NZ::Int32,RL::Float64,TOL::Float64,ELIM::Float64,ALIM::Float64)
    AA::Float64 = zero(Float64)
    AEZ::Float64 = zero(Float64)
    AK::Float64 = zero(Float64)
    AK1I::Float64 = zero(Float64)
    AK1R::Float64 = zero(Float64)
    ARG::Float64 = zero(Float64)
    ARM::Float64 = zero(Float64)
    ATOL::Float64 = zero(Float64)
    AZ::Float64 = zero(Float64)
    BB::Float64 = zero(Float64)
    BK::Float64 = zero(Float64)
    CKI::Float64 = zero(Float64)
    CKR::Float64 = zero(Float64)
    CONEI::Float64 = zero(Float64)
    CONER::Float64 = zero(Float64)
    CS1I::Float64 = zero(Float64)
    CS1R::Float64 = zero(Float64)
    CS2I::Float64 = zero(Float64)
    CS2R::Float64 = zero(Float64)
    CZI::Float64 = zero(Float64)
    CZR::Float64 = zero(Float64)
    DFNU::Float64 = zero(Float64)
    DKI::Float64 = zero(Float64)
    DKR::Float64 = zero(Float64)
    DNU2::Float64 = zero(Float64)
    EZI::Float64 = zero(Float64)
    EZR::Float64 = zero(Float64)
    FDN::Float64 = zero(Float64)
    I::Int32 = zero(Int32)
    IB::Int32 = zero(Int32)
    IL::Int32 = zero(Int32)
    INU::Int32 = zero(Int32)
    J::Int32 = zero(Int32)
    JL::Int32 = zero(Int32)
    K::Int32 = zero(Int32)
    KODED::Int32 = zero(Int32)
    M::Int32 = zero(Int32)
    NN::Int32 = zero(Int32)
    P1I::Float64 = zero(Float64)
    P1R::Float64 = zero(Float64)
    PI::Float64 = zero(Float64)
    RAZ::Float64 = zero(Float64)
    RTPI::Float64 = zero(Float64)
    RTR1::Float64 = zero(Float64)
    RZI::Float64 = zero(Float64)
    RZR::Float64 = zero(Float64)
    S::Float64 = zero(Float64)
    S2I::Float64 = zero(Float64)
    S2R::Float64 = zero(Float64)
    SGN::Float64 = zero(Float64)
    SQK::Float64 = zero(Float64)
    STI::Float64 = zero(Float64)
    STR::Float64 = zero(Float64)
    TZI::Float64 = zero(Float64)
    TZR::Float64 = zero(Float64)
    ZEROI::Float64 = zero(Float64)
    ZEROR::Float64 = zero(Float64)
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
    NZ = int32(0)
    AZ = ZABS(COMPLEX(ZR,ZI))
    ARM = 1000.0D1MACH1
    RTR1 = DSQRT(ARM)
    IL = MIN0(int32(2),N)
    DFNU = FNU + DBLE(FLOAT(N - IL))
    RAZ = 1.0 / AZ
    STR = ZR * RAZ
    STI = -ZI * RAZ
    AK1R = RTPI * STR * RAZ
    AK1I = RTPI * STI * RAZ
    (AK1R,AK1I) = ZSQRT(AK1R,AK1I,AK1R,AK1I)
    CZR = ZR
    CZI = ZI
    if KODE != int32(2)
        @goto line10
    end
    CZR = ZEROR
    CZI = ZI
    @label line10
    if DABS(CZR) > ELIM
        @goto line100
    end
    DNU2 = DFNU + DFNU
    KODED = int32(1)
    if DABS(CZR) > ALIM && N > int32(2)
        @goto line20
    end
    KODED = int32(0)
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
    JL = INT(SNGL(RL + RL)) + int32(2)
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
    if MOD(INU,int32(2)) == int32(0)
        @goto line30
    end
    P1R = -P1R
    P1I = -P1I
    @label line30
    for K = int32(1):IL
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
        for J = int32(1):JL
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
    if N <= int32(2)
        return NZ
    end
    NN = N
    K = NN - int32(2)
    AK = DBLE(FLOAT(K))
    STR = ZR * RAZ
    STI = -ZI * RAZ
    RZR = (STR + STR) * RAZ
    RZI = (STI + STI) * RAZ
    IB = int32(3)
    for I = IB:NN
        YR[K] = (AK + FNU) * (RZR * YR[K + int32(1)] - RZI * YI[K + int32(1)]) + YR[K + int32(2)]
        YI[K] = (AK + FNU) * (RZR * YI[K + int32(1)] + RZI * YR[K + int32(1)]) + YI[K + int32(2)]
        AK = AK - 1.0
        K = K - int32(1)
        @label line80
    end
    if KODED == int32(0)
        return NZ
    end
    (CKR,CKI) = ZEXP(CZR,CZI,CKR,CKI)
    for I = int32(1):NN
        STR = YR[I] * CKR - YI[I] * CKI
        YI[I] = YR[I] * CKI + YI[I] * CKR
        YR[I] = STR
        @label line90
    end
    return NZ
    @label line100
    NZ = int32(-1)
    return NZ
    @label line110
    NZ = int32(-2)
    return NZ
end
