const _ZAIRY_CYR = Array(Float64,1)
const _ZAIRY_CYI = Array(Float64,1)
function ZAIRY(ZR::Float64,ZI::Float64,ID::Integer,KODE::Integer,AIR::Float64,AII::Float64,NZ::Integer,IERR::Integer)
    AA::Float64 = 0
    AD::Float64 = 0
    AK::Float64 = 0
    ALAZ::Float64 = 0
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
    const CYI = _ZAIRY_CYI
    const CYR = _ZAIRY_CYR
    D1::Float64 = 0
    D2::Float64 = 0
    DIG::Float64 = 0
    DK::Float64 = 0
    ELIM::Float64 = 0
    FID::Float64 = 0
    FNU::Float64 = 0
    IFLAG::Int32 = 0
    K::Int32 = 0
    K1::Int32 = 0
    K2::Int32 = 0
    MR::Int32 = 0
    NN::Int32 = 0
    PTR::Float64 = 0
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
    ZEROI::Float64 = 0
    ZEROR::Float64 = 0
    ZTAI::Float64 = 0
    ZTAR::Float64 = 0
    begin 
        TTH = 0.6666666666666666
        C1 = 0.3550280538878172
        C2 = 0.2588194037928068
        COEF = 0.18377629847393068
    end
    begin 
        ZEROR = 0.0
        ZEROI = 0.0
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
        return (AIR,AII,NZ,IERR)
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
        @goto line170
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
    AIR = S1R * C1 - C2 * (ZR * S2R - ZI * S2I)
    AII = S1I * C1 - C2 * (ZR * S2I + ZI * S2R)
    if KODE == 1
        return (AIR,AII,NZ,IERR)
    end
    (STR,STI) = ZSQRT(ZR,ZI,STR,STI)
    ZTAR = TTH * (ZR * STR - ZI * STI)
    ZTAI = TTH * (ZR * STI + ZI * STR)
    (STR,STI) = ZEXP(ZTAR,ZTAI,STR,STI)
    PTR = AIR * STR - AII * STI
    AII = AIR * STI + AII * STR
    AIR = PTR
    return (AIR,AII,NZ,IERR)
    @label line50
    AIR = -S2R * C2
    AII = -S2I * C2
    if AZ <= TOL
        @goto line60
    end
    STR = ZR * S1R - ZI * S1I
    STI = ZR * S1I + ZI * S1R
    CC = C1 / (1.0 + FID)
    AIR = AIR + CC * (STR * ZR - STI * ZI)
    AII = AII + CC * (STR * ZI + STI * ZR)
    @label line60
    if KODE == 1
        return (AIR,AII,NZ,IERR)
    end
    (STR,STI) = ZSQRT(ZR,ZI,STR,STI)
    ZTAR = TTH * (ZR * STR - ZI * STI)
    ZTAI = TTH * (ZR * STI + ZI * STR)
    (STR,STI) = ZEXP(ZTAR,ZTAI,STR,STI)
    PTR = STR * AIR - STI * AII
    AII = STR * AII + STI * AIR
    AIR = PTR
    return (AIR,AII,NZ,IERR)
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
    ALAZ = DLOG(AZ)
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
    IFLAG = 0
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
    if ZI != 0.0
        @goto line90
    end
    if ZR > 0.0
        @goto line90
    end
    ZTAR = 0.0
    ZTAI = AK
    @label line90
    AA = ZTAR
    if AA >= 0.0 && ZR > 0.0
        @goto line110
    end
    if KODE == 2
        @goto line100
    end
    if AA > -ALIM
        @goto line100
    end
    AA = -AA + 0.25ALAZ
    IFLAG = 1
    SFAC = TOL
    if AA > ELIM
        @goto line270
    end
    @label line100
    MR = 1
    if ZI < 0.0
        MR = -1
    end
    (NN,) = ZACAI(ZTAR,ZTAI,FNU,KODE,MR,1,CYR,CYI,NN,RL,TOL,ELIM,ALIM)
    if NN < 0
        @goto line280
    end
    NZ = NZ + NN
    @goto line130
    @label line110
    if KODE == 2
        @goto line120
    end
    if AA < ALIM
        @goto line120
    end
    AA = -AA - 0.25ALAZ
    IFLAG = 2
    SFAC = 1.0 / TOL
    if AA < -ELIM
        @goto line210
    end
    @label line120
    (NZ,) = ZBKNU(ZTAR,ZTAI,FNU,KODE,1,CYR,CYI,NZ,TOL,ELIM,ALIM)
    @label line130
    S1R = CYR[1] * COEF
    S1I = CYI[1] * COEF
    if IFLAG != 0
        @goto line150
    end
    if ID == 1
        @goto line140
    end
    AIR = CSQR * S1R - CSQI * S1I
    AII = CSQR * S1I + CSQI * S1R
    return (AIR,AII,NZ,IERR)
    @label line140
    AIR = -((ZR * S1R - ZI * S1I))
    AII = -((ZR * S1I + ZI * S1R))
    return (AIR,AII,NZ,IERR)
    @label line150
    S1R = S1R * SFAC
    S1I = S1I * SFAC
    if ID == 1
        @goto line160
    end
    STR = S1R * CSQR - S1I * CSQI
    S1I = S1R * CSQI + S1I * CSQR
    S1R = STR
    AIR = S1R / SFAC
    AII = S1I / SFAC
    return (AIR,AII,NZ,IERR)
    @label line160
    STR = -((S1R * ZR - S1I * ZI))
    S1I = -((S1R * ZI + S1I * ZR))
    S1R = STR
    AIR = S1R / SFAC
    AII = S1I / SFAC
    return (AIR,AII,NZ,IERR)
    @label line170
    AA = 1000.0D1MACH1
    S1R = ZEROR
    S1I = ZEROI
    if ID == 1
        @goto line190
    end
    if AZ <= AA
        @goto line180
    end
    S1R = C2 * ZR
    S1I = C2 * ZI
    @label line180
    AIR = C1 - S1R
    AII = -S1I
    return (AIR,AII,NZ,IERR)
    @label line190
    AIR = -C2
    AII = 0.0
    AA = DSQRT(AA)
    if AZ <= AA
        @goto line200
    end
    S1R = 0.5 * (ZR * ZR - ZI * ZI)
    S1I = ZR * ZI
    @label line200
    AIR = AIR + C1 * S1R
    AII = AII + C1 * S1I
    return (AIR,AII,NZ,IERR)
    @label line210
    NZ = 1
    AIR = ZEROR
    AII = ZEROI
    return (AIR,AII,NZ,IERR)
    @label line270
    NZ = 0
    IERR = 2
    return (AIR,AII,NZ,IERR)
    @label line280
    if NN == -1
        @goto line270
    end
    NZ = 0
    IERR = 5
    return (AIR,AII,NZ,IERR)
    @label line260
    IERR = 4
    NZ = 0
    return (AIR,AII,NZ,IERR)
end
