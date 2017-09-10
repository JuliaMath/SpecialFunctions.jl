include("dgamln_constants.jl")

function DGAMLN(Z::Float64,IERR::Int32)
    __DGAMLN__::Float64 = zero(Float64)
    FLN::Float64 = zero(Float64)
    FZ::Float64 = zero(Float64)
    I::Int32 = zero(Int32)
    I1M::Int32 = zero(Int32)
    K::Int32 = zero(Int32)
    MZ::Int32 = zero(Int32)
    NZ::Int32 = zero(Int32)
    RLN::Float64 = zero(Float64)
    S::Float64 = zero(Float64)
    T1::Float64 = zero(Float64)
    TLG::Float64 = zero(Float64)
    TRM::Float64 = zero(Float64)
    TST::Float64 = zero(Float64)
    WDTOL::Float64 = zero(Float64)
    ZDMY::Float64 = zero(Float64)
    ZINC::Float64 = zero(Float64)
    ZM::Float64 = zero(Float64)
    ZMIN::Float64 = zero(Float64)
    ZP::Float64 = zero(Float64)
    ZSQ::Float64 = zero(Float64)
    IERR = Int32(0)
    if Z <= 0.0
        @goto line70
    end
    if Z > 101.0
        @goto line10
    end
    NZ = INT(SNGL(Z))
    FZ = Z - FLOAT(NZ)
    if FZ > 0.0
        @goto line10
    end
    if NZ > Int32(100)
        @goto line10
    end
    __DGAMLN__ = GLN[NZ]
    return __DGAMLN__
    @label line10
    WDTOL = D1MACH4
    WDTOL = DMAX1(WDTOL,5.0e-19)
    I1M = I1MACH14
    RLN = D1MACH5 * FLOAT(I1M)
    FLN = DMIN1(RLN,20.0)
    FLN = DMAX1(FLN,3.0)
    FLN = FLN - 3.0
    ZM = 1.8 + 0.3875FLN
    MZ = INT(SNGL(ZM)) + Int32(1)
    ZMIN = FLOAT(MZ)
    ZDMY = Z
    ZINC = 0.0
    if Z >= ZMIN
        @goto line20
    end
    ZINC = ZMIN - FLOAT(NZ)
    ZDMY = Z + ZINC
    @label line20
    ZP = 1.0 / ZDMY
    T1 = CF[Int32(1)] * ZP
    S = T1
    if ZP < WDTOL
        @goto line40
    end
    ZSQ = ZP * ZP
    TST = T1 * WDTOL
    for K = Int32(2):Int32(22)
        ZP = ZP * ZSQ
        TRM = CF[K] * ZP
        if DABS(TRM) < TST
            @goto line40
        end
        S = S + TRM
        @label line30
    end
    @label line40
    if ZINC != 0.0
        @goto line50
    end
    TLG = DLOG(Z)
    __DGAMLN__ = Z * (TLG - 1.0) + 0.5 * (CON - TLG) + S
    return __DGAMLN__
    @label line50
    ZP = 1.0
    NZ = INT(SNGL(ZINC))
    for I = Int32(1):NZ
        ZP = ZP * (Z + FLOAT(I - Int32(1)))
        @label line60
    end
    TLG = DLOG(ZDMY)
    __DGAMLN__ = (ZDMY * (TLG - 1.0) - DLOG(ZP)) + 0.5 * (CON - TLG) + S
    return __DGAMLN__
    @label line70
    IERR = Int32(1)
    return __DGAMLN__
end
