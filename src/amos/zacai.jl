function ZACAI(ZR::Float64,ZI::Float64,FNU::Float64,KODE::Integer,MR::Integer,N::Integer,YR::AbstractArray{Float64},YI::AbstractArray{Float64},NZ::Integer,RL::Float64,TOL::Float64,ELIM::Float64,ALIM::Float64)
    ARG::Float64 = 0
    ASCLE::Float64 = 0
    AZ::Float64 = 0
    C1I::Float64 = 0
    C1R::Float64 = 0
    C2I::Float64 = 0
    C2R::Float64 = 0
    CSGNI::Float64 = 0
    CSGNR::Float64 = 0
    CSPNI::Float64 = 0
    CSPNR::Float64 = 0
    CYI = Array(Float64,2)
    CYR = Array(Float64,2)
    DFNU::Float64 = 0
    FMR::Float64 = 0
    INU::Int32 = 0
    IUF::Int32 = 0
    NN::Int32 = 0
    NW::Int32 = 0
    PI::Float64 = 0
    SGN::Float64 = 0
    YY::Float64 = 0
    ZNI::Float64 = 0
    ZNR::Float64 = 0
    begin 
        PI = 3.141592653589793
    end
    NZ = 0
    ZNR = -ZR
    ZNI = -ZI
    AZ = ZABS(COMPLEX(ZR,ZI))
    NN = N
    DFNU = FNU + DBLE(FLOAT(N - 1))
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
    if NW < 0
        @goto line80
    end
    @goto line40
    @label line30
    (NW,) = ZMLRI(ZNR,ZNI,FNU,KODE,NN,YR,YI,NW,TOL)
    if NW < 0
        @goto line80
    end
    @label line40
    (NW,) = ZBKNU(ZNR,ZNI,FNU,KODE,1,CYR,CYI,NW,TOL,ELIM,ALIM)
    if NW != 0
        @goto line80
    end
    FMR = DBLE(FLOAT(MR))
    SGN = -(DSIGN(PI,FMR))
    CSGNR = 0.0
    CSGNI = SGN
    if KODE == 1
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
    if MOD(INU,2) == 0
        @goto line60
    end
    CSPNR = -CSPNR
    CSPNI = -CSPNI
    @label line60
    C1R = CYR[1]
    C1I = CYI[1]
    C2R = YR[1]
    C2I = YI[1]
    if KODE == 1
        @goto line70
    end
    IUF = 0
    ASCLE = (1000.0D1MACH1) / TOL
    (C1R,C1I,C2R,C2I,NW,IUF) = ZS1S2(ZNR,ZNI,C1R,C1I,C2R,C2I,NW,ASCLE,ALIM,IUF)
    NZ = NZ + NW
    @label line70
    YR[1] = ((CSPNR * C1R - CSPNI * C1I) + CSGNR * C2R) - CSGNI * C2I
    YI[1] = CSPNR * C1I + CSPNI * C1R + CSGNR * C2I + CSGNI * C2R
    return NZ
    @label line80
    NZ = -1
    if NW == -2
        NZ = -2
    end
    return NZ
end
