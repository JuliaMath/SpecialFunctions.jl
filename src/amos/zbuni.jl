const _ZBUNI_CYR = Array(Float64,2)
const _ZBUNI_CYI = Array(Float64,2)
const _ZBUNI_BRY = Array(Float64,3)
function ZBUNI(ZR::Float64,ZI::Float64,FNU::Float64,KODE::Integer,N::Integer,YR::AbstractArray{Float64},YI::AbstractArray{Float64},NZ::Integer,NUI::Integer,NLAST::Integer,FNUL::Float64,TOL::Float64,ELIM::Float64,ALIM::Float64)
    ASCLE::Float64 = 0
    AX::Float64 = 0
    AY::Float64 = 0
    const BRY = _ZBUNI_BRY
    C1I::Float64 = 0
    C1M::Float64 = 0
    C1R::Float64 = 0
    CSCLR::Float64 = 0
    CSCRR::Float64 = 0
    const CYI = _ZBUNI_CYI
    const CYR = _ZBUNI_CYR
    DFNU::Float64 = 0
    FNUI::Float64 = 0
    GNU::Float64 = 0
    I::Int32 = 0
    IFLAG::Int32 = 0
    IFORM::Int32 = 0
    K::Int32 = 0
    NL::Int32 = 0
    NW::Int32 = 0
    RAZ::Float64 = 0
    RZI::Float64 = 0
    RZR::Float64 = 0
    S1I::Float64 = 0
    S1R::Float64 = 0
    S2I::Float64 = 0
    S2R::Float64 = 0
    STI::Float64 = 0
    STR::Float64 = 0
    NZ = 0
    AX = DABS(ZR) * 1.7321
    AY = DABS(ZI)
    IFORM = 1
    if AY > AX
        IFORM = 2
    end
    if NUI == 0
        @goto line60
    end
    FNUI = DBLE(FLOAT(NUI))
    DFNU = FNU + DBLE(FLOAT(N - 1))
    GNU = DFNU + FNUI
    if IFORM == 2
        @goto line10
    end
    (NW,NLAST) = ZUNI1(ZR,ZI,GNU,KODE,2,CYR,CYI,NW,NLAST,FNUL,TOL,ELIM,ALIM)
    @goto line20
    @label line10
    (NW,NLAST) = ZUNI2(ZR,ZI,GNU,KODE,2,CYR,CYI,NW,NLAST,FNUL,TOL,ELIM,ALIM)
    @label line20
    if NW < 0
        @goto line50
    end
    if NW != 0
        @goto line90
    end
    STR = ZABS(COMPLEX(CYR[1],CYI[1]))
    BRY[1] = (1000.0D1MACH1) / TOL
    BRY[2] = 1.0 / BRY[1]
    BRY[3] = BRY[2]
    IFLAG = 2
    ASCLE = BRY[2]
    CSCLR = 1.0
    if STR > BRY[1]
        @goto line21
    end
    IFLAG = 1
    ASCLE = BRY[1]
    CSCLR = 1.0 / TOL
    @goto line25
    @label line21
    if STR < BRY[2]
        @goto line25
    end
    IFLAG = 3
    ASCLE = BRY[3]
    CSCLR = TOL
    @label line25
    CSCRR = 1.0 / CSCLR
    S1R = CYR[2] * CSCLR
    S1I = CYI[2] * CSCLR
    S2R = CYR[1] * CSCLR
    S2I = CYI[1] * CSCLR
    RAZ = 1.0 / ZABS(COMPLEX(ZR,ZI))
    STR = ZR * RAZ
    STI = -ZI * RAZ
    RZR = (STR + STR) * RAZ
    RZI = (STI + STI) * RAZ
    for I = 1:NUI
        STR = S2R
        STI = S2I
        S2R = (DFNU + FNUI) * (RZR * STR - RZI * STI) + S1R
        S2I = (DFNU + FNUI) * (RZR * STI + RZI * STR) + S1I
        S1R = STR
        S1I = STI
        FNUI = FNUI - 1.0
        if IFLAG >= 3
            @goto line30
        end
        STR = S2R * CSCRR
        STI = S2I * CSCRR
        C1R = DABS(STR)
        C1I = DABS(STI)
        C1M = DMAX1(C1R,C1I)
        if C1M <= ASCLE
            @goto line30
        end
        IFLAG = IFLAG + 1
        ASCLE = BRY[IFLAG]
        S1R = S1R * CSCRR
        S1I = S1I * CSCRR
        S2R = STR
        S2I = STI
        CSCLR = CSCLR * TOL
        CSCRR = 1.0 / CSCLR
        S1R = S1R * CSCLR
        S1I = S1I * CSCLR
        S2R = S2R * CSCLR
        S2I = S2I * CSCLR
        @label line30
    end
    YR[N] = S2R * CSCRR
    YI[N] = S2I * CSCRR
    if N == 1
        return (NZ,NLAST)
    end
    NL = N - 1
    FNUI = DBLE(FLOAT(NL))
    K = NL
    for I = 1:NL
        STR = S2R
        STI = S2I
        S2R = (FNU + FNUI) * (RZR * STR - RZI * STI) + S1R
        S2I = (FNU + FNUI) * (RZR * STI + RZI * STR) + S1I
        S1R = STR
        S1I = STI
        STR = S2R * CSCRR
        STI = S2I * CSCRR
        YR[K] = STR
        YI[K] = STI
        FNUI = FNUI - 1.0
        K = K - 1
        if IFLAG >= 3
            @goto line40
        end
        C1R = DABS(STR)
        C1I = DABS(STI)
        C1M = DMAX1(C1R,C1I)
        if C1M <= ASCLE
            @goto line40
        end
        IFLAG = IFLAG + 1
        ASCLE = BRY[IFLAG]
        S1R = S1R * CSCRR
        S1I = S1I * CSCRR
        S2R = STR
        S2I = STI
        CSCLR = CSCLR * TOL
        CSCRR = 1.0 / CSCLR
        S1R = S1R * CSCLR
        S1I = S1I * CSCLR
        S2R = S2R * CSCLR
        S2I = S2I * CSCLR
        @label line40
    end
    return (NZ,NLAST)
    @label line50
    NZ = -1
    if NW == -2
        NZ = -2
    end
    return (NZ,NLAST)
    @label line60
    if IFORM == 2
        @goto line70
    end
    (NW,NLAST) = ZUNI1(ZR,ZI,FNU,KODE,N,YR,YI,NW,NLAST,FNUL,TOL,ELIM,ALIM)
    @goto line80
    @label line70
    (NW,NLAST) = ZUNI2(ZR,ZI,FNU,KODE,N,YR,YI,NW,NLAST,FNUL,TOL,ELIM,ALIM)
    @label line80
    if NW < 0
        @goto line50
    end
    NZ = NW
    return (NZ,NLAST)
    @label line90
    NLAST = N
    return (NZ,NLAST)
end
