function ZS1S2(ZRR::Float64,ZRI::Float64,S1R::Float64,S1I::Float64,S2R::Float64,S2I::Float64,NZ::Int32,ASCLE::Float64,ALIM::Float64,IUF::Int32)
    AA::Float64 = zero(Float64)
    ALN::Float64 = zero(Float64)
    AS1::Float64 = zero(Float64)
    AS2::Float64 = zero(Float64)
    C1I::Float64 = zero(Float64)
    C1R::Float64 = zero(Float64)
    IDUM::Int32 = zero(Int32)
    S1DI::Float64 = zero(Float64)
    S1DR::Float64 = zero(Float64)
    ZEROI::Float64 = zero(Float64)
    ZEROR::Float64 = zero(Float64)
    begin 
        ZEROR = 0.0
        ZEROI = 0.0
    end
    NZ = int32(0)
    AS1 = ZABS(COMPLEX(S1R,S1I))
    AS2 = ZABS(COMPLEX(S2R,S2I))
    if S1R == 0.0 && S1I == 0.0
        @goto line10
    end
    if AS1 == 0.0
        @goto line10
    end
    ALN = (-ZRR - ZRR) + DLOG(AS1)
    S1DR = S1R
    S1DI = S1I
    S1R = ZEROR
    S1I = ZEROI
    AS1 = ZEROR
    if ALN < -ALIM
        @goto line10
    end
    (C1R,C1I,IDUM) = ZLOG(S1DR,S1DI,C1R,C1I,IDUM)
    C1R = (C1R - ZRR) - ZRR
    C1I = (C1I - ZRI) - ZRI
    (S1R,S1I) = ZEXP(C1R,C1I,S1R,S1I)
    AS1 = ZABS(COMPLEX(S1R,S1I))
    IUF = IUF + int32(1)
    @label line10
    AA = DMAX1(AS1,AS2)
    if AA > ASCLE
        return (S1R,S1I,S2R,S2I,NZ,IUF)
    end
    S1R = ZEROR
    S1I = ZEROI
    S2R = ZEROR
    S2I = ZEROI
    NZ = int32(1)
    IUF = int32(0)
    return (S1R,S1I,S2R,S2I,NZ,IUF)
end
