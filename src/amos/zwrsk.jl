function ZWRSK(ZRR::Float64,ZRI::Float64,FNU::Float64,KODE::Integer,N::Integer,YR::AbstractArray{Float64},YI::AbstractArray{Float64},NZ::Integer,CWR::AbstractArray{Float64},CWI::AbstractArray{Float64},TOL::Float64,ELIM::Float64,ALIM::Float64)
    ACT::Float64 = 0
    ACW::Float64 = 0
    ASCLE::Float64 = 0
    C1I::Float64 = 0
    C1R::Float64 = 0
    C2I::Float64 = 0
    C2R::Float64 = 0
    CINUI::Float64 = 0
    CINUR::Float64 = 0
    CSCLR::Float64 = 0
    CTI::Float64 = 0
    CTR::Float64 = 0
    I::Int32 = 0
    NW::Int32 = 0
    PTI::Float64 = 0
    PTR::Float64 = 0
    RACT::Float64 = 0
    STI::Float64 = 0
    STR::Float64 = 0
    NZ = 0
    (NW,) = ZBKNU(ZRR,ZRI,FNU,KODE,2,CWR,CWI,NW,TOL,ELIM,ALIM)
    if NW != 0
        @goto line50
    end
    ZRATI(ZRR,ZRI,FNU,N,YR,YI,TOL)
    CINUR = 1.0
    CINUI = 0.0
    if KODE == 1
        @goto line10
    end
    CINUR = DCOS(ZRI)
    CINUI = DSIN(ZRI)
    @label line10
    ACW = ZABS(COMPLEX(CWR[2],CWI[2]))
    ASCLE = (1000.0D1MACH1) / TOL
    CSCLR = 1.0
    if ACW > ASCLE
        @goto line20
    end
    CSCLR = 1.0 / TOL
    @goto line30
    @label line20
    ASCLE = 1.0 / ASCLE
    if ACW < ASCLE
        @goto line30
    end
    CSCLR = TOL
    @label line30
    C1R = CWR[1] * CSCLR
    C1I = CWI[1] * CSCLR
    C2R = CWR[2] * CSCLR
    C2I = CWI[2] * CSCLR
    STR = YR[1]
    STI = YI[1]
    PTR = STR * C1R - STI * C1I
    PTI = STR * C1I + STI * C1R
    PTR = PTR + C2R
    PTI = PTI + C2I
    CTR = ZRR * PTR - ZRI * PTI
    CTI = ZRR * PTI + ZRI * PTR
    ACT = ZABS(COMPLEX(CTR,CTI))
    RACT = 1.0 / ACT
    CTR = CTR * RACT
    CTI = -CTI * RACT
    PTR = CINUR * RACT
    PTI = CINUI * RACT
    CINUR = PTR * CTR - PTI * CTI
    CINUI = PTR * CTI + PTI * CTR
    YR[1] = CINUR * CSCLR
    YI[1] = CINUI * CSCLR
    if N == 1
        return NZ
    end
    for I = 2:N
        PTR = STR * CINUR - STI * CINUI
        CINUI = STR * CINUI + STI * CINUR
        CINUR = PTR
        STR = YR[I]
        STI = YI[I]
        YR[I] = CINUR * CSCLR
        YI[I] = CINUI * CSCLR
        @label line40
    end
    return NZ
    @label line50
    NZ = -1
    if NW == -2
        NZ = -2
    end
    return NZ
end
