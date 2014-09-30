const _ZKSCL_CYR = Array(Float64,2)
const _ZKSCL_CYI = Array(Float64,2)
function ZKSCL(ZRR::Float64,ZRI::Float64,FNU::Float64,N::Integer,YR::AbstractArray{Float64},YI::AbstractArray{Float64},NZ::Integer,RZR::Float64,RZI::Float64,ASCLE::Float64,TOL::Float64,ELIM::Float64)
    ACS::Float64 = 0
    ALAS::Float64 = 0
    AS::Float64 = 0
    CELMR::Float64 = 0
    CKI::Float64 = 0
    CKR::Float64 = 0
    CSI::Float64 = 0
    CSR::Float64 = 0
    const CYI = _ZKSCL_CYI
    const CYR = _ZKSCL_CYR
    ELM::Float64 = 0
    FN::Float64 = 0
    HELIM::Float64 = 0
    I::Int32 = 0
    IC::Int32 = 0
    IDUM::Int32 = 0
    KK::Int32 = 0
    NN::Int32 = 0
    NW::Int32 = 0
    S1I::Float64 = 0
    S1R::Float64 = 0
    S2I::Float64 = 0
    S2R::Float64 = 0
    STR::Float64 = 0
    ZDI::Float64 = 0
    ZDR::Float64 = 0
    ZEROI::Float64 = 0
    ZEROR::Float64 = 0
    begin 
        ZEROR = 0.0
        ZEROI = 0.0
    end
    NZ = 0
    IC = 0
    NN = MIN0(2,N)
    for I = 1:NN
        S1R = YR[I]
        S1I = YI[I]
        CYR[I] = S1R
        CYI[I] = S1I
        AS = ZABS(COMPLEX(S1R,S1I))
        ACS = -ZRR + DLOG(AS)
        NZ = NZ + 1
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
        if NW != 0
            @goto line10
        end
        YR[I] = CSR
        YI[I] = CSI
        IC = I
        NZ = NZ - 1
        @label line10
    end
    if N == 1
        return NZ
    end
    if IC > 1
        @goto line20
    end
    YR[1] = ZEROR
    YI[1] = ZEROI
    NZ = 2
    @label line20
    if N == 2
        return NZ
    end
    if NZ == 0
        return NZ
    end
    FN = FNU + 1.0
    CKR = FN * RZR
    CKI = FN * RZI
    S1R = CYR[1]
    S1I = CYI[1]
    S2R = CYR[2]
    S2I = CYI[2]
    HELIM = 0.5ELIM
    ELM = DEXP(-ELIM)
    CELMR = ELM
    ZDR = ZRR
    ZDI = ZRI
    for I = 3:N
        KK = I
        CSR = S2R
        CSI = S2I
        S2R = (CKR * CSR - CKI * CSI) + S1R
        S2I = CKI * CSR + CKR * CSI + S1I
        S1R = CSR
        S1I = CSI
        CKR = CKR + RZR
        CKI = CKI + RZI
        AS = ZABS(COMPLEX(S2R,S2I))
        ALAS = DLOG(AS)
        ACS = -ZDR + ALAS
        NZ = NZ + 1
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
        if NW != 0
            @goto line25
        end
        YR[I] = CSR
        YI[I] = CSI
        NZ = NZ - 1
        if IC == KK - 1
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
        NZ = N - 1
    end
    @goto line45
    @label line40
    NZ = KK - 2
    @label line45
    for I = 1:NZ
        YR[I] = ZEROR
        YI[I] = ZEROI
        @label line50
    end
    return NZ
end
