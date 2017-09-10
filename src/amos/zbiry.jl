 _ZBIRY_CYR = Array{Float64}(2)
 _ZBIRY_CYI = Array{Float64}(2)
function ZBIRY(ZR::Float64,ZI::Float64,ID::Int32,KODE::Int32,BIR::Float64,BII::Float64,IERR::Int32)
    AA::Float64 = zero(Float64)
    AD::Float64 = zero(Float64)
    AK::Float64 = zero(Float64)
    ALIM::Float64 = zero(Float64)
    ATRM::Float64 = zero(Float64)
    AZ::Float64 = zero(Float64)
    AZ3::Float64 = zero(Float64)
    BB::Float64 = zero(Float64)
    BK::Float64 = zero(Float64)
    C1::Float64 = zero(Float64)
    C2::Float64 = zero(Float64)
    CC::Float64 = zero(Float64)
    CK::Float64 = zero(Float64)
    COEF::Float64 = zero(Float64)
    CONEI::Float64 = zero(Float64)
    CONER::Float64 = zero(Float64)
    CSQI::Float64 = zero(Float64)
    CSQR::Float64 = zero(Float64)
     CYI = _ZBIRY_CYI
     CYR = _ZBIRY_CYR
    D1::Float64 = zero(Float64)
    D2::Float64 = zero(Float64)
    DIG::Float64 = zero(Float64)
    DK::Float64 = zero(Float64)
    EAA::Float64 = zero(Float64)
    ELIM::Float64 = zero(Float64)
    FID::Float64 = zero(Float64)
    FMR::Float64 = zero(Float64)
    FNU::Float64 = zero(Float64)
    FNUL::Float64 = zero(Float64)
    K::Int32 = zero(Int32)
    K1::Int32 = zero(Int32)
    K2::Int32 = zero(Int32)
    NZ::Int32 = zero(Int32)
    PI::Float64 = zero(Float64)
    R1M5::Float64 = zero(Float64)
    RL::Float64 = zero(Float64)
    S1I::Float64 = zero(Float64)
    S1R::Float64 = zero(Float64)
    S2I::Float64 = zero(Float64)
    S2R::Float64 = zero(Float64)
    SFAC::Float64 = zero(Float64)
    STI::Float64 = zero(Float64)
    STR::Float64 = zero(Float64)
    TOL::Float64 = zero(Float64)
    TRM1I::Float64 = zero(Float64)
    TRM1R::Float64 = zero(Float64)
    TRM2I::Float64 = zero(Float64)
    TRM2R::Float64 = zero(Float64)
    TTH::Float64 = zero(Float64)
    Z3I::Float64 = zero(Float64)
    Z3R::Float64 = zero(Float64)
    ZTAI::Float64 = zero(Float64)
    ZTAR::Float64 = zero(Float64)
    begin 
        TTH = 0.6666666666666666
        C1 = 0.6149266274460007
        C2 = 0.4482883573538264
        COEF = 0.5773502691896257
        PI = 3.141592653589793
    end
    begin 
        CONER = 1.0
        CONEI = 0.0
    end
    IERR = int32(0)
    NZ = int32(0)
    if ID < int32(0) || ID > int32(1)
        IERR = int32(1)
    end
    if KODE < int32(1) || KODE > int32(2)
        IERR = int32(1)
    end
    if IERR != int32(0)
        return (BIR,BII,IERR)
    end
    AZ = ZABS(COMPLEX(ZR,ZI))
    TOL = DMAX1(D1MACH4,1.0e-18)
    FID = DBLE(FLOAT(ID))
    if AZ > 1.0
        @goto line70
    end
    S1R = CONER
    S1I = CONEI
    S2R = CONER
    S2I = CONEI
    if AZ < TOL
        @goto line130
    end
    AA = AZ * AZ
    if AA < TOL / AZ
        @goto line40
    end
    TRM1R = CONER
    TRM1I = CONEI
    TRM2R = CONER
    TRM2I = CONEI
    ATRM = 1.0
    STR = ZR * ZR - ZI * ZI
    STI = ZR * ZI + ZI * ZR
    Z3R = STR * ZR - STI * ZI
    Z3I = STR * ZI + STI * ZR
    AZ3 = AZ * AA
    AK = 2.0 + FID
    BK = (3.0 - FID) - FID
    CK = 4.0 - FID
    DK = 3.0 + FID + FID
    D1 = AK * DK
    D2 = BK * CK
    AD = DMIN1(D1,D2)
    AK = 24.0 + 9.0FID
    BK = 30.0 - 9.0FID
    for K = int32(1):int32(25)
        STR = (TRM1R * Z3R - TRM1I * Z3I) / D1
        TRM1I = (TRM1R * Z3I + TRM1I * Z3R) / D1
        TRM1R = STR
        S1R = S1R + TRM1R
        S1I = S1I + TRM1I
        STR = (TRM2R * Z3R - TRM2I * Z3I) / D2
        TRM2I = (TRM2R * Z3I + TRM2I * Z3R) / D2
        TRM2R = STR
        S2R = S2R + TRM2R
        S2I = S2I + TRM2I
        ATRM = (ATRM * AZ3) / AD
        D1 = D1 + AK
        D2 = D2 + BK
        AD = DMIN1(D1,D2)
        if ATRM < TOL * AD
            @goto line40
        end
        AK = AK + 18.0
        BK = BK + 18.0
        @label line30
    end
    @label line40
    if ID == int32(1)
        @goto line50
    end
    BIR = C1 * S1R + C2 * (ZR * S2R - ZI * S2I)
    BII = C1 * S1I + C2 * (ZR * S2I + ZI * S2R)
    if KODE == int32(1)
        return (BIR,BII,IERR)
    end
    (STR,STI) = ZSQRT(ZR,ZI,STR,STI)
    ZTAR = TTH * (ZR * STR - ZI * STI)
    ZTAI = TTH * (ZR * STI + ZI * STR)
    AA = ZTAR
    AA = -(DABS(AA))
    EAA = DEXP(AA)
    BIR = BIR * EAA
    BII = BII * EAA
    return (BIR,BII,IERR)
    @label line50
    BIR = S2R * C2
    BII = S2I * C2
    if AZ <= TOL
        @goto line60
    end
    CC = C1 / (1.0 + FID)
    STR = S1R * ZR - S1I * ZI
    STI = S1R * ZI + S1I * ZR
    BIR = BIR + CC * (STR * ZR - STI * ZI)
    BII = BII + CC * (STR * ZI + STI * ZR)
    @label line60
    if KODE == int32(1)
        return (BIR,BII,IERR)
    end
    (STR,STI) = ZSQRT(ZR,ZI,STR,STI)
    ZTAR = TTH * (ZR * STR - ZI * STI)
    ZTAI = TTH * (ZR * STI + ZI * STR)
    AA = ZTAR
    AA = -(DABS(AA))
    EAA = DEXP(AA)
    BIR = BIR * EAA
    BII = BII * EAA
    return (BIR,BII,IERR)
    @label line70
    FNU = (1.0 + FID) / 3.0
    K1 = I1MACH15
    K2 = I1MACH16
    R1M5 = D1MACH5
    K = MIN0(IABS(K1),IABS(K2))
    ELIM = 2.303 * (DBLE(FLOAT(K)) * R1M5 - 3.0)
    K1 = I1MACH14 - int32(1)
    AA = R1M5 * DBLE(FLOAT(K1))
    DIG = DMIN1(AA,18.0)
    AA = AA * 2.303
    ALIM = ELIM + DMAX1(-AA,-41.45)
    RL = 1.2DIG + 3.0
    FNUL = 10.0 + 6.0 * (DIG - 3.0)
    AA = 0.5 / TOL
    BB = DBLE(FLOAT(I1MACH9)) * 0.5
    AA = DMIN1(AA,BB)
    AA = AA^TTH
    if AZ > AA
        @goto line260
    end
    AA = DSQRT(AA)
    if AZ > AA
        IERR = int32(3)
    end
    (CSQR,CSQI) = ZSQRT(ZR,ZI,CSQR,CSQI)
    ZTAR = TTH * (ZR * CSQR - ZI * CSQI)
    ZTAI = TTH * (ZR * CSQI + ZI * CSQR)
    SFAC = 1.0
    AK = ZTAI
    if ZR >= 0.0
        @goto line80
    end
    BK = ZTAR
    CK = -(DABS(BK))
    ZTAR = CK
    ZTAI = AK
    @label line80
    if ZI != 0.0 || ZR > 0.0
        @goto line90
    end
    ZTAR = 0.0
    ZTAI = AK
    @label line90
    AA = ZTAR
    if KODE == int32(2)
        @goto line100
    end
    BB = DABS(AA)
    if BB < ALIM
        @goto line100
    end
    BB = BB + 0.25 * DLOG(AZ)
    SFAC = TOL
    if BB > ELIM
        @goto line190
    end
    @label line100
    FMR = 0.0
    if AA >= 0.0 && ZR > 0.0
        @goto line110
    end
    FMR = PI
    if ZI < 0.0
        FMR = -PI
    end
    ZTAR = -ZTAR
    ZTAI = -ZTAI
    @label line110
    (NZ,) = ZBINU(ZTAR,ZTAI,FNU,KODE,int32(1),CYR,CYI,NZ,RL,FNUL,TOL,ELIM,ALIM)
    if NZ < int32(0)
        @goto line200
    end
    AA = FMR * FNU
    Z3R = SFAC
    STR = DCOS(AA)
    STI = DSIN(AA)
    S1R = (STR * CYR[int32(1)] - STI * CYI[int32(1)]) * Z3R
    S1I = (STR * CYI[int32(1)] + STI * CYR[int32(1)]) * Z3R
    FNU = (2.0 - FID) / 3.0
    (NZ,) = ZBINU(ZTAR,ZTAI,FNU,KODE,int32(2),CYR,CYI,NZ,RL,FNUL,TOL,ELIM,ALIM)
    CYR[int32(1)] = CYR[int32(1)] * Z3R
    CYI[int32(1)] = CYI[int32(1)] * Z3R
    CYR[int32(2)] = CYR[int32(2)] * Z3R
    CYI[int32(2)] = CYI[int32(2)] * Z3R
    (STR,STI) = ZDIV(CYR[int32(1)],CYI[int32(1)],ZTAR,ZTAI,STR,STI)
    S2R = (FNU + FNU) * STR + CYR[int32(2)]
    S2I = (FNU + FNU) * STI + CYI[int32(2)]
    AA = FMR * (FNU - 1.0)
    STR = DCOS(AA)
    STI = DSIN(AA)
    S1R = COEF * ((S1R + S2R * STR) - S2I * STI)
    S1I = COEF * (S1I + S2R * STI + S2I * STR)
    if ID == int32(1)
        @goto line120
    end
    STR = CSQR * S1R - CSQI * S1I
    S1I = CSQR * S1I + CSQI * S1R
    S1R = STR
    BIR = S1R / SFAC
    BII = S1I / SFAC
    return (BIR,BII,IERR)
    @label line120
    STR = ZR * S1R - ZI * S1I
    S1I = ZR * S1I + ZI * S1R
    S1R = STR
    BIR = S1R / SFAC
    BII = S1I / SFAC
    return (BIR,BII,IERR)
    @label line130
    AA = C1 * (1.0 - FID) + FID * C2
    BIR = AA
    BII = 0.0
    return (BIR,BII,IERR)
    @label line190
    IERR = int32(2)
    NZ = int32(0)
    return (BIR,BII,IERR)
    @label line200
    if NZ == int32(-1)
        @goto line190
    end
    NZ = int32(0)
    IERR = int32(5)
    return (BIR,BII,IERR)
    @label line260
    IERR = int32(4)
    NZ = int32(0)
    return (BIR,BII,IERR)
end
