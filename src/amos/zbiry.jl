function ZBIRY(ZR::Float64,ZI::Float64,ID::Integer,KODE::Integer,BIR::Float64,BII::Float64,IERR::Integer)
    AA::Float64 = 0
    AD::Float64 = 0
    AK::Float64 = 0
    ALIM::Float64 = 0
    ATRM::Float64 = 0
    AZ::Float64 = 0
    AZ3::Float64 = 0
    BB::Float64 = 0
    BK::Float64 = 0
    C1::Float64 = 0
    C2::Float64 = 0
    CC::Float64 = 0
    CK::Float64 = 0
    COEF::Float64 = 0
    CONEI::Float64 = 0
    CONER::Float64 = 0
    CSQI::Float64 = 0
    CSQR::Float64 = 0
    CYI = Array(Float64,2)
    CYR = Array(Float64,2)
    D1::Float64 = 0
    D2::Float64 = 0
    DIG::Float64 = 0
    DK::Float64 = 0
    EAA::Float64 = 0
    ELIM::Float64 = 0
    FID::Float64 = 0
    FMR::Float64 = 0
    FNU::Float64 = 0
    FNUL::Float64 = 0
    K::Int32 = 0
    K1::Int32 = 0
    K2::Int32 = 0
    NZ::Int32 = 0
    PI::Float64 = 0
    R1M5::Float64 = 0
    RL::Float64 = 0
    S1I::Float64 = 0
    S1R::Float64 = 0
    S2I::Float64 = 0
    S2R::Float64 = 0
    SFAC::Float64 = 0
    STI::Float64 = 0
    STR::Float64 = 0
    TOL::Float64 = 0
    TRM1I::Float64 = 0
    TRM1R::Float64 = 0
    TRM2I::Float64 = 0
    TRM2R::Float64 = 0
    TTH::Float64 = 0
    Z3I::Float64 = 0
    Z3R::Float64 = 0
    ZTAI::Float64 = 0
    ZTAR::Float64 = 0
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
    IERR = 0
    NZ = 0
    if ID < 0 || ID > 1
        IERR = 1
    end
    if KODE < 1 || KODE > 2
        IERR = 1
    end
    if IERR != 0
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
    for K = 1:25
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
    if ID == 1
        @goto line50
    end
    BIR = C1 * S1R + C2 * (ZR * S2R - ZI * S2I)
    BII = C1 * S1I + C2 * (ZR * S2I + ZI * S2R)
    if KODE == 1
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
    if KODE == 1
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
    K1 = I1MACH14 - 1
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
        IERR = 3
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
    if KODE == 2
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
    (NZ,) = ZBINU(ZTAR,ZTAI,FNU,KODE,1,CYR,CYI,NZ,RL,FNUL,TOL,ELIM,ALIM)
    if NZ < 0
        @goto line200
    end
    AA = FMR * FNU
    Z3R = SFAC
    STR = DCOS(AA)
    STI = DSIN(AA)
    S1R = (STR * CYR[1] - STI * CYI[1]) * Z3R
    S1I = (STR * CYI[1] + STI * CYR[1]) * Z3R
    FNU = (2.0 - FID) / 3.0
    (NZ,) = ZBINU(ZTAR,ZTAI,FNU,KODE,2,CYR,CYI,NZ,RL,FNUL,TOL,ELIM,ALIM)
    CYR[1] = CYR[1] * Z3R
    CYI[1] = CYI[1] * Z3R
    CYR[2] = CYR[2] * Z3R
    CYI[2] = CYI[2] * Z3R
    (STR,STI) = ZDIV(CYR[1],CYI[1],ZTAR,ZTAI,STR,STI)
    S2R = (FNU + FNU) * STR + CYR[2]
    S2I = (FNU + FNU) * STI + CYI[2]
    AA = FMR * (FNU - 1.0)
    STR = DCOS(AA)
    STI = DSIN(AA)
    S1R = COEF * ((S1R + S2R * STR) - S2I * STI)
    S1I = COEF * (S1I + S2R * STI + S2I * STR)
    if ID == 1
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
    IERR = 2
    NZ = 0
    return (BIR,BII,IERR)
    @label line200
    if NZ == -1
        @goto line190
    end
    NZ = 0
    IERR = 5
    return (BIR,BII,IERR)
    @label line260
    IERR = 4
    NZ = 0
    return (BIR,BII,IERR)
end
