function ZBESK(ZR::Float64,ZI::Float64,FNU::Float64,KODE::Int32,N::Int32,CYR::AbstractArray{Float64},CYI::AbstractArray{Float64},NZ::Int32,IERR::Int32)
    AA::Float64 = zero(Float64)
    ALIM::Float64 = zero(Float64)
    ALN::Float64 = zero(Float64)
    ARG::Float64 = zero(Float64)
    AZ::Float64 = zero(Float64)
    BB::Float64 = zero(Float64)
    DIG::Float64 = zero(Float64)
    ELIM::Float64 = zero(Float64)
    FN::Float64 = zero(Float64)
    FNUL::Float64 = zero(Float64)
    K::Int32 = zero(Int32)
    K1::Int32 = zero(Int32)
    K2::Int32 = zero(Int32)
    MR::Int32 = zero(Int32)
    NN::Int32 = zero(Int32)
    NUF::Int32 = zero(Int32)
    NW::Int32 = zero(Int32)
    R1M5::Float64 = zero(Float64)
    RL::Float64 = zero(Float64)
    TOL::Float64 = zero(Float64)
    UFL::Float64 = zero(Float64)
    IERR = Int32(0)
    NZ = Int32(0)
    if ZI == 0.0 && ZR == 0.0
        IERR = Int32(1)
    end
    if FNU < 0.0
        IERR = Int32(1)
    end
    if KODE < Int32(1) || KODE > Int32(2)
        IERR = Int32(1)
    end
    if N < Int32(1)
        IERR = Int32(1)
    end
    if IERR != Int32(0)
        return (NZ,IERR)
    end
    NN = N
    TOL = DMAX1(D1MACH4,1.0e-18)
    K1 = I1MACH15
    K2 = I1MACH16
    R1M5 = D1MACH5
    K = MIN0(IABS(K1),IABS(K2))
    ELIM = 2.303 * (DBLE(FLOAT(K)) * R1M5 - 3.0)
    K1 = I1MACH14 - Int32(1)
    AA = R1M5 * DBLE(FLOAT(K1))
    DIG = DMIN1(AA,18.0)
    AA = AA * 2.303
    ALIM = ELIM + DMAX1(-AA,-41.45)
    FNUL = 10.0 + 6.0 * (DIG - 3.0)
    RL = 1.2DIG + 3.0
    AZ = abs(COMPLEX(ZR,ZI))
    FN = FNU + DBLE(FLOAT(NN - Int32(1)))
    AA = 0.5 / TOL
    BB = DBLE(FLOAT(I1MACH9)) * 0.5
    AA = DMIN1(AA,BB)
    if AZ > AA
        @goto line260
    end
    if FN > AA
        @goto line260
    end
    AA = DSQRT(AA)
    if AZ > AA
        IERR = Int32(3)
    end
    if FN > AA
        IERR = Int32(3)
    end
    UFL = D1MACH1 * 1000.0
    if AZ < UFL
        @goto line180
    end
    if FNU > FNUL
        @goto line80
    end
    if FN <= 1.0
        @goto line60
    end
    if FN > 2.0
        @goto line50
    end
    if AZ > TOL
        @goto line60
    end
    ARG = 0.5AZ
    ALN = -FN * DLOG(ARG)
    if ALN > ELIM
        @goto line180
    end
    @goto line60
    @label line50
    (NUF,) = ZUOIK(ZR,ZI,FNU,KODE,Int32(2),NN,CYR,CYI,NUF,TOL,ELIM,ALIM)
    if NUF < Int32(0)
        @goto line180
    end
    NZ = NZ + NUF
    NN = NN - NUF
    if NN == Int32(0)
        @goto line100
    end
    @label line60
    if ZR < 0.0
        @goto line70
    end
    (NW,) = ZBKNU(ZR,ZI,FNU,KODE,NN,CYR,CYI,NW,TOL,ELIM,ALIM)
    if NW < Int32(0)
        @goto line200
    end
    NZ = NW
    return (NZ,IERR)
    @label line70
    if NZ != Int32(0)
        @goto line180
    end
    MR = Int32(1)
    if ZI < 0.0
        MR = Int32(-1)
    end
    (NW,) = ZACON(ZR,ZI,FNU,KODE,MR,NN,CYR,CYI,NW,RL,FNUL,TOL,ELIM,ALIM)
    if NW < Int32(0)
        @goto line200
    end
    NZ = NW
    return (NZ,IERR)
    @label line80
    MR = Int32(0)
    if ZR >= 0.0
        @goto line90
    end
    MR = Int32(1)
    if ZI < 0.0
        MR = Int32(-1)
    end
    @label line90
    (NW,) = ZBUNK(ZR,ZI,FNU,KODE,MR,NN,CYR,CYI,NW,TOL,ELIM,ALIM)
    if NW < Int32(0)
        @goto line200
    end
    NZ = NZ + NW
    return (NZ,IERR)
    @label line100
    if ZR < 0.0
        @goto line180
    end
    return (NZ,IERR)
    @label line180
    NZ = Int32(0)
    IERR = Int32(2)
    return (NZ,IERR)
    @label line200
    if NW == Int32(-1)
        @goto line180
    end
    NZ = Int32(0)
    IERR = Int32(5)
    return (NZ,IERR)
    @label line260
    NZ = Int32(0)
    IERR = Int32(4)
    return (NZ,IERR)
end
