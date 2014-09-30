const _ZACAI_CYR = Array(Float64,2)
const _ZACAI_CYI = Array(Float64,2)
function ZACAI(ZR::Float64,ZI::Float64,FNU::Float64,KODE::Int32,MR::Int32,N::Int32,YR::AbstractArray{Float64},YI::AbstractArray{Float64},NZ::Int32,RL::Float64,TOL::Float64,ELIM::Float64,ALIM::Float64)
    ARG::Float64 = zero(Float64)
    ASCLE::Float64 = zero(Float64)
    AZ::Float64 = zero(Float64)
    C1I::Float64 = zero(Float64)
    C1R::Float64 = zero(Float64)
    C2I::Float64 = zero(Float64)
    C2R::Float64 = zero(Float64)
    CSGNI::Float64 = zero(Float64)
    CSGNR::Float64 = zero(Float64)
    CSPNI::Float64 = zero(Float64)
    CSPNR::Float64 = zero(Float64)
    const CYI = _ZACAI_CYI
    const CYR = _ZACAI_CYR
    DFNU::Float64 = zero(Float64)
    FMR::Float64 = zero(Float64)
    INU::Int32 = zero(Int32)
    IUF::Int32 = zero(Int32)
    NN::Int32 = zero(Int32)
    NW::Int32 = zero(Int32)
    PI::Float64 = zero(Float64)
    SGN::Float64 = zero(Float64)
    YY::Float64 = zero(Float64)
    ZNI::Float64 = zero(Float64)
    ZNR::Float64 = zero(Float64)
    begin 
        PI = 3.141592653589793
    end
    NZ = int32(0)
    ZNR = -ZR
    ZNI = -ZI
    AZ = ZABS(COMPLEX(ZR,ZI))
    NN = N
    DFNU = FNU + DBLE(FLOAT(N - int32(1)))
    if AZ <= 2.0
        @goto line10
    end
    if AZ * AZ * 0.25 > DFNU + 1.0
        @goto line20
    end
    @label line10
    (NW,) = ZSERI(ZNR,ZNI,FNU,KODE,NN,YR,YI,NW,TOL,ELIM,ALIM)
    @goto line40
    @label line20
    if AZ < RL
        @goto line30
    end
    (NW,) = ZASYI(ZNR,ZNI,FNU,KODE,NN,YR,YI,NW,RL,TOL,ELIM,ALIM)
    if NW < int32(0)
        @goto line80
    end
    @goto line40
    @label line30
    (NW,) = ZMLRI(ZNR,ZNI,FNU,KODE,NN,YR,YI,NW,TOL)
    if NW < int32(0)
        @goto line80
    end
    @label line40
    (NW,) = ZBKNU(ZNR,ZNI,FNU,KODE,int32(1),CYR,CYI,NW,TOL,ELIM,ALIM)
    if NW != int32(0)
        @goto line80
    end
    FMR = DBLE(FLOAT(MR))
    SGN = -(DSIGN(PI,FMR))
    CSGNR = 0.0
    CSGNI = SGN
    if KODE == int32(1)
        @goto line50
    end
    YY = -ZNI
    CSGNR = -CSGNI * DSIN(YY)
    CSGNI = CSGNI * DCOS(YY)
    @label line50
    INU = INT(SNGL(FNU))
    ARG = (FNU - DBLE(FLOAT(INU))) * SGN
    CSPNR = DCOS(ARG)
    CSPNI = DSIN(ARG)
    if MOD(INU,int32(2)) == int32(0)
        @goto line60
    end
    CSPNR = -CSPNR
    CSPNI = -CSPNI
    @label line60
    C1R = CYR[int32(1)]
    C1I = CYI[int32(1)]
    C2R = YR[int32(1)]
    C2I = YI[int32(1)]
    if KODE == int32(1)
        @goto line70
    end
    IUF = int32(0)
    ASCLE = (1000.0D1MACH1) / TOL
    (C1R,C1I,C2R,C2I,NW,IUF) = ZS1S2(ZNR,ZNI,C1R,C1I,C2R,C2I,NW,ASCLE,ALIM,IUF)
    NZ = NZ + NW
    @label line70
    YR[int32(1)] = ((CSPNR * C1R - CSPNI * C1I) + CSGNR * C2R) - CSGNI * C2I
    YI[int32(1)] = CSPNR * C1I + CSPNI * C1R + CSGNR * C2I + CSGNI * C2R
    return NZ
    @label line80
    NZ = int32(-1)
    if NW == int32(-2)
        NZ = int32(-2)
    end
    return NZ
end
