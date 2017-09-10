 _ZUOIK_CWRKR = Array{Float64}(16)
 _ZUOIK_CWRKI = Array{Float64}(16)
function ZUOIK(ZR::Float64,ZI::Float64,FNU::Float64,KODE::Int32,IKFLG::Int32,N::Int32,YR::AbstractArray{Float64},YI::AbstractArray{Float64},NUF::Int32,TOL::Float64,ELIM::Float64,ALIM::Float64)
    AARG::Float64 = zero(Float64)
    AIC::Float64 = zero(Float64)
    APHI::Float64 = zero(Float64)
    ARGI::Float64 = zero(Float64)
    ARGR::Float64 = zero(Float64)
    ASCLE::Float64 = zero(Float64)
    ASUMI::Float64 = zero(Float64)
    ASUMR::Float64 = zero(Float64)
    AX::Float64 = zero(Float64)
    AY::Float64 = zero(Float64)
    BSUMI::Float64 = zero(Float64)
    BSUMR::Float64 = zero(Float64)
     CWRKI = _ZUOIK_CWRKI
     CWRKR = _ZUOIK_CWRKR
    CZI::Float64 = zero(Float64)
    CZR::Float64 = zero(Float64)
    FNN::Float64 = zero(Float64)
    GNN::Float64 = zero(Float64)
    GNU::Float64 = zero(Float64)
    I::Int32 = zero(Int32)
    IDUM::Int32 = zero(Int32)
    IFORM::Int32 = zero(Int32)
    INIT::Int32 = zero(Int32)
    NN::Int32 = zero(Int32)
    NW::Int32 = zero(Int32)
    PHII::Float64 = zero(Float64)
    PHIR::Float64 = zero(Float64)
    RCZ::Float64 = zero(Float64)
    STI::Float64 = zero(Float64)
    STR::Float64 = zero(Float64)
    SUMI::Float64 = zero(Float64)
    SUMR::Float64 = zero(Float64)
    ZBI::Float64 = zero(Float64)
    ZBR::Float64 = zero(Float64)
    ZEROI::Float64 = zero(Float64)
    ZEROR::Float64 = zero(Float64)
    ZETA1I::Float64 = zero(Float64)
    ZETA1R::Float64 = zero(Float64)
    ZETA2I::Float64 = zero(Float64)
    ZETA2R::Float64 = zero(Float64)
    ZNI::Float64 = zero(Float64)
    ZNR::Float64 = zero(Float64)
    ZRI::Float64 = zero(Float64)
    ZRR::Float64 = zero(Float64)
    begin 
        ZEROR = 0.0
        ZEROI = 0.0
    end
    begin 
        AIC = 1.2655121234846454
    end
    NUF = int32(0)
    NN = N
    ZRR = ZR
    ZRI = ZI
    if ZR >= 0.0
        @goto line10
    end
    ZRR = -ZR
    ZRI = -ZI
    @label line10
    ZBR = ZRR
    ZBI = ZRI
    AX = DABS(ZR) * 1.7321
    AY = DABS(ZI)
    IFORM = int32(1)
    if AY > AX
        IFORM = int32(2)
    end
    GNU = DMAX1(FNU,1.0)
    if IKFLG == int32(1)
        @goto line20
    end
    FNN = DBLE(FLOAT(NN))
    GNN = (FNU + FNN) - 1.0
    GNU = DMAX1(GNN,FNN)
    @label line20
    if IFORM == int32(2)
        @goto line30
    end
    INIT = int32(0)
    (INIT,PHIR,PHII,ZETA1R,ZETA1I,ZETA2R,ZETA2I,SUMR,SUMI) = ZUNIK(ZRR,ZRI,GNU,IKFLG,int32(1),TOL,INIT,PHIR,PHII,ZETA1R,ZETA1I,ZETA2R,ZETA2I,SUMR,SUMI,CWRKR,CWRKI)
    CZR = -ZETA1R + ZETA2R
    CZI = -ZETA1I + ZETA2I
    @goto line50
    @label line30
    ZNR = ZRI
    ZNI = -ZRR
    if ZI > 0.0
        @goto line40
    end
    ZNR = -ZNR
    @label line40
    (PHIR,PHII,ARGR,ARGI,ZETA1R,ZETA1I,ZETA2R,ZETA2I,ASUMR,ASUMI,BSUMR,BSUMI) = ZUNHJ(ZNR,ZNI,GNU,int32(1),TOL,PHIR,PHII,ARGR,ARGI,ZETA1R,ZETA1I,ZETA2R,ZETA2I,ASUMR,ASUMI,BSUMR,BSUMI)
    CZR = -ZETA1R + ZETA2R
    CZI = -ZETA1I + ZETA2I
    AARG = ZABS(COMPLEX(ARGR,ARGI))
    @label line50
    if KODE == int32(1)
        @goto line60
    end
    CZR = CZR - ZBR
    CZI = CZI - ZBI
    @label line60
    if IKFLG == int32(1)
        @goto line70
    end
    CZR = -CZR
    CZI = -CZI
    @label line70
    APHI = ZABS(COMPLEX(PHIR,PHII))
    RCZ = CZR
    if RCZ > ELIM
        @goto line210
    end
    if RCZ < ALIM
        @goto line80
    end
    RCZ = RCZ + DLOG(APHI)
    if IFORM == int32(2)
        RCZ = (RCZ - 0.25 * DLOG(AARG)) - AIC
    end
    if RCZ > ELIM
        @goto line210
    end
    @goto line130
    @label line80
    if RCZ < -ELIM
        @goto line90
    end
    if RCZ > -ALIM
        @goto line130
    end
    RCZ = RCZ + DLOG(APHI)
    if IFORM == int32(2)
        RCZ = (RCZ - 0.25 * DLOG(AARG)) - AIC
    end
    if RCZ > -ELIM
        @goto line110
    end
    @label line90
    for I = int32(1):NN
        YR[I] = ZEROR
        YI[I] = ZEROI
        @label line100
    end
    NUF = NN
    return NUF
    @label line110
    ASCLE = (1000.0D1MACH1) / TOL
    (STR,STI,IDUM) = ZLOG(PHIR,PHII,STR,STI,IDUM)
    CZR = CZR + STR
    CZI = CZI + STI
    if IFORM == int32(1)
        @goto line120
    end
    (STR,STI,IDUM) = ZLOG(ARGR,ARGI,STR,STI,IDUM)
    CZR = (CZR - 0.25STR) - AIC
    CZI = CZI - 0.25STI
    @label line120
    AX = DEXP(RCZ) / TOL
    AY = CZI
    CZR = AX * DCOS(AY)
    CZI = AX * DSIN(AY)
    (NW,) = ZUCHK(CZR,CZI,NW,ASCLE,TOL)
    if NW != int32(0)
        @goto line90
    end
    @label line130
    if IKFLG == int32(2)
        return NUF
    end
    if N == int32(1)
        return NUF
    end
    @label line140
    GNU = FNU + DBLE(FLOAT(NN - int32(1)))
    if IFORM == int32(2)
        @goto line150
    end
    INIT = int32(0)
    (INIT,PHIR,PHII,ZETA1R,ZETA1I,ZETA2R,ZETA2I,SUMR,SUMI) = ZUNIK(ZRR,ZRI,GNU,IKFLG,int32(1),TOL,INIT,PHIR,PHII,ZETA1R,ZETA1I,ZETA2R,ZETA2I,SUMR,SUMI,CWRKR,CWRKI)
    CZR = -ZETA1R + ZETA2R
    CZI = -ZETA1I + ZETA2I
    @goto line160
    @label line150
    (PHIR,PHII,ARGR,ARGI,ZETA1R,ZETA1I,ZETA2R,ZETA2I,ASUMR,ASUMI,BSUMR,BSUMI) = ZUNHJ(ZNR,ZNI,GNU,int32(1),TOL,PHIR,PHII,ARGR,ARGI,ZETA1R,ZETA1I,ZETA2R,ZETA2I,ASUMR,ASUMI,BSUMR,BSUMI)
    CZR = -ZETA1R + ZETA2R
    CZI = -ZETA1I + ZETA2I
    AARG = ZABS(COMPLEX(ARGR,ARGI))
    @label line160
    if KODE == int32(1)
        @goto line170
    end
    CZR = CZR - ZBR
    CZI = CZI - ZBI
    @label line170
    APHI = ZABS(COMPLEX(PHIR,PHII))
    RCZ = CZR
    if RCZ < -ELIM
        @goto line180
    end
    if RCZ > -ALIM
        return NUF
    end
    RCZ = RCZ + DLOG(APHI)
    if IFORM == int32(2)
        RCZ = (RCZ - 0.25 * DLOG(AARG)) - AIC
    end
    if RCZ > -ELIM
        @goto line190
    end
    @label line180
    YR[NN] = ZEROR
    YI[NN] = ZEROI
    NN = NN - int32(1)
    NUF = NUF + int32(1)
    if NN == int32(0)
        return NUF
    end
    @goto line140
    @label line190
    ASCLE = (1000.0D1MACH1) / TOL
    (STR,STI,IDUM) = ZLOG(PHIR,PHII,STR,STI,IDUM)
    CZR = CZR + STR
    CZI = CZI + STI
    if IFORM == int32(1)
        @goto line200
    end
    (STR,STI,IDUM) = ZLOG(ARGR,ARGI,STR,STI,IDUM)
    CZR = (CZR - 0.25STR) - AIC
    CZI = CZI - 0.25STI
    @label line200
    AX = DEXP(RCZ) / TOL
    AY = CZI
    CZR = AX * DCOS(AY)
    CZI = AX * DSIN(AY)
    (NW,) = ZUCHK(CZR,CZI,NW,ASCLE,TOL)
    if NW != int32(0)
        @goto line180
    end
    return NUF
    @label line210
    NUF = int32(-1)
    return NUF
end
