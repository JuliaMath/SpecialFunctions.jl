const _ZKSCL_CYR = Array{Float64}(2)
const _ZKSCL_CYI = Array{Float64}(2)
function ZKSCL(ZRR::Float64,ZRI::Float64,FNU::Float64,N::Int32,YR::AbstractArray{Float64},YI::AbstractArray{Float64},NZ::Int32,RZR::Float64,RZI::Float64,ASCLE::Float64,TOL::Float64,ELIM::Float64)
    ACS::Float64 = zero(Float64)
    ALAS::Float64 = zero(Float64)
    AS::Float64 = zero(Float64)
    CELMR::Float64 = zero(Float64)
    CKI::Float64 = zero(Float64)
    CKR::Float64 = zero(Float64)
    CSI::Float64 = zero(Float64)
    CSR::Float64 = zero(Float64)
    const CYI = _ZKSCL_CYI
    const CYR = _ZKSCL_CYR
    ELM::Float64 = zero(Float64)
    FN::Float64 = zero(Float64)
    HELIM::Float64 = zero(Float64)
    I::Int32 = zero(Int32)
    IC::Int32 = zero(Int32)
    IDUM::Int32 = zero(Int32)
    KK::Int32 = zero(Int32)
    NN::Int32 = zero(Int32)
    NW::Int32 = zero(Int32)
    S1I::Float64 = zero(Float64)
    S1R::Float64 = zero(Float64)
    S2I::Float64 = zero(Float64)
    S2R::Float64 = zero(Float64)
    STR::Float64 = zero(Float64)
    ZDI::Float64 = zero(Float64)
    ZDR::Float64 = zero(Float64)
    ZEROI::Float64 = zero(Float64)
    ZEROR::Float64 = zero(Float64)
    begin 
        ZEROR = 0.0
        ZEROI = 0.0
    end
    NZ = Int32(0)
    IC = Int32(0)
    NN = MIN0(Int32(2),N)
    for I = Int32(1):NN
        S1R = YR[I]
        S1I = YI[I]
        CYR[I] = S1R
        CYI[I] = S1I
        AS = abs(COMPLEX(S1R,S1I))
        ACS = -ZRR + DLOG(AS)
        NZ = NZ + Int32(1)
        YR[I] = ZEROR
        YI[I] = ZEROI
        if ACS < -ELIM
            @goto line10
        end
        (CSR,CSI,IDUM) = ZLOG(S1R,S1I,CSR,CSI,IDUM)
        CSR = CSR - ZRR
        CSI = CSI - ZRI
        STR = DEXP(CSR) / TOL
        CSR = STR * DCOS(CSI)
        CSI = STR * DSIN(CSI)
        (NW,) = ZUCHK(CSR,CSI,NW,ASCLE,TOL)
        if NW != Int32(0)
            @goto line10
        end
        YR[I] = CSR
        YI[I] = CSI
        IC = I
        NZ = NZ - Int32(1)
        @label line10
    end
    if N == Int32(1)
        return NZ
    end
    if IC > Int32(1)
        @goto line20
    end
    YR[Int32(1)] = ZEROR
    YI[Int32(1)] = ZEROI
    NZ = Int32(2)
    @label line20
    if N == Int32(2)
        return NZ
    end
    if NZ == Int32(0)
        return NZ
    end
    FN = FNU + 1.0
    CKR = FN * RZR
    CKI = FN * RZI
    S1R = CYR[Int32(1)]
    S1I = CYI[Int32(1)]
    S2R = CYR[Int32(2)]
    S2I = CYI[Int32(2)]
    HELIM = 0.5ELIM
    ELM = DEXP(-ELIM)
    CELMR = ELM
    ZDR = ZRR
    ZDI = ZRI
    for I = Int32(3):N
        KK = I
        CSR = S2R
        CSI = S2I
        S2R = (CKR * CSR - CKI * CSI) + S1R
        S2I = CKI * CSR + CKR * CSI + S1I
        S1R = CSR
        S1I = CSI
        CKR = CKR + RZR
        CKI = CKI + RZI
        AS = abs(COMPLEX(S2R,S2I))
        ALAS = DLOG(AS)
        ACS = -ZDR + ALAS
        NZ = NZ + Int32(1)
        YR[I] = ZEROR
        YI[I] = ZEROI
        if ACS < -ELIM
            @goto line25
        end
        (CSR,CSI,IDUM) = ZLOG(S2R,S2I,CSR,CSI,IDUM)
        CSR = CSR - ZDR
        CSI = CSI - ZDI
        STR = DEXP(CSR) / TOL
        CSR = STR * DCOS(CSI)
        CSI = STR * DSIN(CSI)
        (NW,) = ZUCHK(CSR,CSI,NW,ASCLE,TOL)
        if NW != Int32(0)
            @goto line25
        end
        YR[I] = CSR
        YI[I] = CSI
        NZ = NZ - Int32(1)
        if IC == KK - Int32(1)
            @goto line40
        end
        IC = KK
        @goto line30
        @label line25
        if ALAS < HELIM
            @goto line30
        end
        ZDR = ZDR - ELIM
        S1R = S1R * CELMR
        S1I = S1I * CELMR
        S2R = S2R * CELMR
        S2I = S2I * CELMR
        @label line30
    end
    NZ = N
    if IC == N
        NZ = N - Int32(1)
    end
    @goto line45
    @label line40
    NZ = KK - Int32(2)
    @label line45
    for I = Int32(1):NZ
        YR[I] = ZEROR
        YI[I] = ZEROI
        @label line50
    end
    return NZ
end
