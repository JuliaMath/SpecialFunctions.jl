const _ZACON_CYR = Array{Float64}(2)
const _ZACON_CYI = Array{Float64}(2)
const _ZACON_CSSR = Array{Float64}(3)
const _ZACON_CSRR = Array{Float64}(3)
const _ZACON_BRY = Array{Float64}(3)
function ZACON(ZR::Float64,ZI::Float64,FNU::Float64,KODE::Int32,MR::Int32,N::Int32,YR::AbstractArray{Float64},YI::AbstractArray{Float64},NZ::Int32,RL::Float64,FNUL::Float64,TOL::Float64,ELIM::Float64,ALIM::Float64)
    ARG::Float64 = zero(Float64)
    AS2::Float64 = zero(Float64)
    ASCLE::Float64 = zero(Float64)
    AZN::Float64 = zero(Float64)
    const BRY = _ZACON_BRY
    BSCLE::Float64 = zero(Float64)
    C1I::Float64 = zero(Float64)
    C1M::Float64 = zero(Float64)
    C1R::Float64 = zero(Float64)
    C2I::Float64 = zero(Float64)
    C2R::Float64 = zero(Float64)
    CKI::Float64 = zero(Float64)
    CKR::Float64 = zero(Float64)
    CONER::Float64 = zero(Float64)
    CPN::Float64 = zero(Float64)
    CSCL::Float64 = zero(Float64)
    CSCR::Float64 = zero(Float64)
    CSGNI::Float64 = zero(Float64)
    CSGNR::Float64 = zero(Float64)
    CSPNI::Float64 = zero(Float64)
    CSPNR::Float64 = zero(Float64)
    CSR::Float64 = zero(Float64)
    const CSRR = _ZACON_CSRR
    const CSSR = _ZACON_CSSR
    const CYI = _ZACON_CYI
    const CYR = _ZACON_CYR
    FMR::Float64 = zero(Float64)
    FN::Float64 = zero(Float64)
    I::Int32 = zero(Int32)
    INU::Int32 = zero(Int32)
    IUF::Int32 = zero(Int32)
    KFLAG::Int32 = zero(Int32)
    NN::Int32 = zero(Int32)
    NW::Int32 = zero(Int32)
    PI::Float64 = zero(Float64)
    PTI::Float64 = zero(Float64)
    PTR::Float64 = zero(Float64)
    RAZN::Float64 = zero(Float64)
    RZI::Float64 = zero(Float64)
    RZR::Float64 = zero(Float64)
    S1I::Float64 = zero(Float64)
    S1R::Float64 = zero(Float64)
    S2I::Float64 = zero(Float64)
    S2R::Float64 = zero(Float64)
    SC1I::Float64 = zero(Float64)
    SC1R::Float64 = zero(Float64)
    SC2I::Float64 = zero(Float64)
    SC2R::Float64 = zero(Float64)
    SGN::Float64 = zero(Float64)
    SPN::Float64 = zero(Float64)
    STI::Float64 = zero(Float64)
    STR::Float64 = zero(Float64)
    YY::Float64 = zero(Float64)
    ZEROR::Float64 = zero(Float64)
    ZNI::Float64 = zero(Float64)
    ZNR::Float64 = zero(Float64)
    begin 
        PI = 3.141592653589793
    end
    begin 
        ZEROR = 0.0
        CONER = 1.0
    end
    NZ = Int32(0)
    ZNR = -ZR
    ZNI = -ZI
    NN = N
    (NW,) = ZBINU(ZNR,ZNI,FNU,KODE,NN,YR,YI,NW,RL,FNUL,TOL,ELIM,ALIM)
    if NW < Int32(0)
        @goto line90
    end
    NN = min(Int32(2),N)
    (NW,) = ZBKNU(ZNR,ZNI,FNU,KODE,NN,CYR,CYI,NW,TOL,ELIM,ALIM)
    if NW != Int32(0)
        @goto line90
    end
    S1R = CYR[Int32(1)]
    S1I = CYI[Int32(1)]
    FMR = DBLE(FLOAT(MR))
    SGN = -(copysign(PI,FMR))
    CSGNR = ZEROR
    CSGNI = SGN
    if KODE == Int32(1)
        @goto line10
    end
    YY = -ZNI
    CPN = cos(YY)
    SPN = sin(YY)
    (CSGNR,CSGNI) = ZMLT(CSGNR,CSGNI,CPN,SPN,CSGNR,CSGNI)
    @label line10
    INU = trunc(Int32,SNGL(FNU))
    ARG = (FNU - DBLE(FLOAT(INU))) * SGN
    CPN = cos(ARG)
    SPN = sin(ARG)
    CSPNR = CPN
    CSPNI = SPN
    if mod(INU,Int32(2)) == Int32(0)
        @goto line20
    end
    CSPNR = -CSPNR
    CSPNI = -CSPNI
    @label line20
    IUF = Int32(0)
    C1R = S1R
    C1I = S1I
    C2R = YR[Int32(1)]
    C2I = YI[Int32(1)]
    ASCLE = (1000.0D1MACH1) / TOL
    if KODE == Int32(1)
        @goto line30
    end
    (C1R,C1I,C2R,C2I,NW,IUF) = ZS1S2(ZNR,ZNI,C1R,C1I,C2R,C2I,NW,ASCLE,ALIM,IUF)
    NZ = NZ + NW
    SC1R = C1R
    SC1I = C1I
    @label line30
    (STR,STI) = ZMLT(CSPNR,CSPNI,C1R,C1I,STR,STI)
    (PTR,PTI) = ZMLT(CSGNR,CSGNI,C2R,C2I,PTR,PTI)
    YR[Int32(1)] = STR + PTR
    YI[Int32(1)] = STI + PTI
    if N == Int32(1)
        return NZ
    end
    CSPNR = -CSPNR
    CSPNI = -CSPNI
    S2R = CYR[Int32(2)]
    S2I = CYI[Int32(2)]
    C1R = S2R
    C1I = S2I
    C2R = YR[Int32(2)]
    C2I = YI[Int32(2)]
    if KODE == Int32(1)
        @goto line40
    end
    (C1R,C1I,C2R,C2I,NW,IUF) = ZS1S2(ZNR,ZNI,C1R,C1I,C2R,C2I,NW,ASCLE,ALIM,IUF)
    NZ = NZ + NW
    SC2R = C1R
    SC2I = C1I
    @label line40
    (STR,STI) = ZMLT(CSPNR,CSPNI,C1R,C1I,STR,STI)
    (PTR,PTI) = ZMLT(CSGNR,CSGNI,C2R,C2I,PTR,PTI)
    YR[Int32(2)] = STR + PTR
    YI[Int32(2)] = STI + PTI
    if N == Int32(2)
        return NZ
    end
    CSPNR = -CSPNR
    CSPNI = -CSPNI
    AZN = abs(complex(ZNR,ZNI))
    RAZN = 1.0 / AZN
    STR = ZNR * RAZN
    STI = -ZNI * RAZN
    RZR = (STR + STR) * RAZN
    RZI = (STI + STI) * RAZN
    FN = FNU + 1.0
    CKR = FN * RZR
    CKI = FN * RZI
    CSCL = 1.0 / TOL
    CSCR = TOL
    CSSR[Int32(1)] = CSCL
    CSSR[Int32(2)] = CONER
    CSSR[Int32(3)] = CSCR
    CSRR[Int32(1)] = CSCR
    CSRR[Int32(2)] = CONER
    CSRR[Int32(3)] = CSCL
    BRY[Int32(1)] = ASCLE
    BRY[Int32(2)] = 1.0 / ASCLE
    BRY[Int32(3)] = D1MACH2
    AS2 = abs(complex(S2R,S2I))
    KFLAG = Int32(2)
    if AS2 > BRY[Int32(1)]
        @goto line50
    end
    KFLAG = Int32(1)
    @goto line60
    @label line50
    if AS2 < BRY[Int32(2)]
        @goto line60
    end
    KFLAG = Int32(3)
    @label line60
    BSCLE = BRY[KFLAG]
    S1R = S1R * CSSR[KFLAG]
    S1I = S1I * CSSR[KFLAG]
    S2R = S2R * CSSR[KFLAG]
    S2I = S2I * CSSR[KFLAG]
    CSR = CSRR[KFLAG]
    for I = Int32(3):N
        STR = S2R
        STI = S2I
        S2R = (CKR * STR - CKI * STI) + S1R
        S2I = CKR * STI + CKI * STR + S1I
        S1R = STR
        S1I = STI
        C1R = S2R * CSR
        C1I = S2I * CSR
        STR = C1R
        STI = C1I
        C2R = YR[I]
        C2I = YI[I]
        if KODE == Int32(1)
            @goto line70
        end
        if IUF < Int32(0)
            @goto line70
        end
        (C1R,C1I,C2R,C2I,NW,IUF) = ZS1S2(ZNR,ZNI,C1R,C1I,C2R,C2I,NW,ASCLE,ALIM,IUF)
        NZ = NZ + NW
        SC1R = SC2R
        SC1I = SC2I
        SC2R = C1R
        SC2I = C1I
        if IUF != Int32(3)
            @goto line70
        end
        IUF = Int32(-4)
        S1R = SC1R * CSSR[KFLAG]
        S1I = SC1I * CSSR[KFLAG]
        S2R = SC2R * CSSR[KFLAG]
        S2I = SC2I * CSSR[KFLAG]
        STR = SC2R
        STI = SC2I
        @label line70
        PTR = CSPNR * C1R - CSPNI * C1I
        PTI = CSPNR * C1I + CSPNI * C1R
        YR[I] = (PTR + CSGNR * C2R) - CSGNI * C2I
        YI[I] = PTI + CSGNR * C2I + CSGNI * C2R
        CKR = CKR + RZR
        CKI = CKI + RZI
        CSPNR = -CSPNR
        CSPNI = -CSPNI
        if KFLAG >= Int32(3)
            @goto line80
        end
        PTR = abs(C1R)
        PTI = abs(C1I)
        C1M = max(PTR,PTI)
        if C1M <= BSCLE
            @goto line80
        end
        KFLAG = KFLAG + Int32(1)
        BSCLE = BRY[KFLAG]
        S1R = S1R * CSR
        S1I = S1I * CSR
        S2R = STR
        S2I = STI
        S1R = S1R * CSSR[KFLAG]
        S1I = S1I * CSSR[KFLAG]
        S2R = S2R * CSSR[KFLAG]
        S2I = S2I * CSSR[KFLAG]
        CSR = CSRR[KFLAG]
        @label line80
    end
    return NZ
    @label line90
    NZ = Int32(-1)
    if NW == Int32(-2)
        NZ = Int32(-2)
    end
    return NZ
end
