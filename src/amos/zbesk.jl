function ZBESK(ZR::Float64,ZI::Float64,FNU::Float64,KODE::Integer,N::Integer,CYR::AbstractArray{Float64},CYI::AbstractArray{Float64},NZ::Integer,IERR::Integer)
    AA::Float64 = 0
    ALIM::Float64 = 0
    ALN::Float64 = 0
    ARG::Float64 = 0
    AZ::Float64 = 0
    BB::Float64 = 0
    DIG::Float64 = 0
    ELIM::Float64 = 0
    FN::Float64 = 0
    FNUL::Float64 = 0
    K::Int32 = 0
    K1::Int32 = 0
    K2::Int32 = 0
    MR::Int32 = 0
    NN::Int32 = 0
    NUF::Int32 = 0
    NW::Int32 = 0
    R1M5::Float64 = 0
    RL::Float64 = 0
    TOL::Float64 = 0
    UFL::Float64 = 0
    IERR = 0
    NZ = 0
    if ZI == 0.0 && ZR == 0.0
        IERR = 1
    end
    if FNU < 0.0
        IERR = 1
    end
    if KODE < 1 || KODE > 2
        IERR = 1
    end
    if N < 1
        IERR = 1
    end
    if IERR != 0
        return (NZ,IERR)
    end
    NN = N
    TOL = DMAX1(D1MACH4,1.0e-18)
    K1 = I1MACH15
    K2 = I1MACH16
    R1M5 = D1MACH5
    K = MIN0(IABS(K1),IABS(K2))
    ELIM = 2.303 * (DBLE(FLOAT(K)) * R1M5 - 3.0)
    K1 = I1MACH14 - 1
    AA = R1M5 * DBLE(FLOAT(K1))
    DIG = DMIN1(AA,18.0)
    AA = AA * 2.303
    ALIM = ELIM + DMAX1(-AA,-41.45)
    FNUL = 10.0 + 6.0 * (DIG - 3.0)
    RL = 1.2DIG + 3.0
    AZ = ZABS(COMPLEX(ZR,ZI))
    FN = FNU + DBLE(FLOAT(NN - 1))
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
        IERR = 3
    end
    if FN > AA
        IERR = 3
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
    (NUF,) = ZUOIK(ZR,ZI,FNU,KODE,2,NN,CYR,CYI,NUF,TOL,ELIM,ALIM)
    if NUF < 0
        @goto line180
    end
    NZ = NZ + NUF
    NN = NN - NUF
    if NN == 0
        @goto line100
    end
    @label line60
    if ZR < 0.0
        @goto line70
    end
    (NW,) = ZBKNU(ZR,ZI,FNU,KODE,NN,CYR,CYI,NW,TOL,ELIM,ALIM)
    if NW < 0
        @goto line200
    end
    NZ = NW
    return (NZ,IERR)
    @label line70
    if NZ != 0
        @goto line180
    end
    MR = 1
    if ZI < 0.0
        MR = -1
    end
    (NW,) = ZACON(ZR,ZI,FNU,KODE,MR,NN,CYR,CYI,NW,RL,FNUL,TOL,ELIM,ALIM)
    if NW < 0
        @goto line200
    end
    NZ = NW
    return (NZ,IERR)
    @label line80
    MR = 0
    if ZR >= 0.0
        @goto line90
    end
    MR = 1
    if ZI < 0.0
        MR = -1
    end
    @label line90
    (NW,) = ZBUNK(ZR,ZI,FNU,KODE,MR,NN,CYR,CYI,NW,TOL,ELIM,ALIM)
    if NW < 0
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
    NZ = 0
    IERR = 2
    return (NZ,IERR)
    @label line200
    if NW == -1
        @goto line180
    end
    NZ = 0
    IERR = 5
    return (NZ,IERR)
    @label line260
    NZ = 0
    IERR = 4
    return (NZ,IERR)
end
