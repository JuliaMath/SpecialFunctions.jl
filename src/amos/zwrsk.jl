function ZWRSK(ZRR::Float64,ZRI::Float64,FNU::Float64,KODE::Int32,N::Int32,YR::AbstractArray{Float64},YI::AbstractArray{Float64},NZ::Int32,CWR::AbstractArray{Float64},CWI::AbstractArray{Float64},TOL::Float64,ELIM::Float64,ALIM::Float64)
    ACT::Float64 = zero(Float64)
    ACW::Float64 = zero(Float64)
    ASCLE::Float64 = zero(Float64)
    C1I::Float64 = zero(Float64)
    C1R::Float64 = zero(Float64)
    C2I::Float64 = zero(Float64)
    C2R::Float64 = zero(Float64)
    CINUI::Float64 = zero(Float64)
    CINUR::Float64 = zero(Float64)
    CSCLR::Float64 = zero(Float64)
    CTI::Float64 = zero(Float64)
    CTR::Float64 = zero(Float64)
    I::Int32 = zero(Int32)
    NW::Int32 = zero(Int32)
    PTI::Float64 = zero(Float64)
    PTR::Float64 = zero(Float64)
    RACT::Float64 = zero(Float64)
    STI::Float64 = zero(Float64)
    STR::Float64 = zero(Float64)
    NZ = Int32(0)
    (NW,) = ZBKNU(ZRR,ZRI,FNU,KODE,Int32(2),CWR,CWI,NW,TOL,ELIM,ALIM)
    if NW != Int32(0)
        @goto line50
    end
    ZRATI(ZRR,ZRI,FNU,N,YR,YI,TOL)
    CINUR = 1.0
    CINUI = 0.0
    if KODE == Int32(1)
        @goto line10
    end
    CINUR = cos(ZRI)
    CINUI = sin(ZRI)
    @label line10
    ACW = abs(complex(CWR[Int32(2)],CWI[Int32(2)]))
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
    C1R = CWR[Int32(1)] * CSCLR
    C1I = CWI[Int32(1)] * CSCLR
    C2R = CWR[Int32(2)] * CSCLR
    C2I = CWI[Int32(2)] * CSCLR
    STR = YR[Int32(1)]
    STI = YI[Int32(1)]
    PTR = STR * C1R - STI * C1I
    PTI = STR * C1I + STI * C1R
    PTR = PTR + C2R
    PTI = PTI + C2I
    CTR = ZRR * PTR - ZRI * PTI
    CTI = ZRR * PTI + ZRI * PTR
    ACT = abs(complex(CTR,CTI))
    RACT = 1.0 / ACT
    CTR = CTR * RACT
    CTI = -CTI * RACT
    PTR = CINUR * RACT
    PTI = CINUI * RACT
    CINUR = PTR * CTR - PTI * CTI
    CINUI = PTR * CTI + PTI * CTR
    YR[Int32(1)] = CINUR * CSCLR
    YI[Int32(1)] = CINUI * CSCLR
    if N == Int32(1)
        return NZ
    end
    for I = Int32(2):N
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
    NZ = Int32(-1)
    if NW == Int32(-2)
        NZ = Int32(-2)
    end
    return NZ
end
