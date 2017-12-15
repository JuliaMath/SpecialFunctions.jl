function ZBESH(ZR::Float64,ZI::Float64,FNU::Float64,KODE::Int32,M::Int32,N::Int32,CYR::AbstractArray{Float64},CYI::AbstractArray{Float64},NZ::Int32,IERR::Int32)
    AA::Float64 = zero(Float64)
    ALIM::Float64 = zero(Float64)
    ALN::Float64 = zero(Float64)
    ARG::Float64 = zero(Float64)
    ASCLE::Float64 = zero(Float64)
    ATOL::Float64 = zero(Float64)
    AZ::Float64 = zero(Float64)
    BB::Float64 = zero(Float64)
    CSGNI::Float64 = zero(Float64)
    CSGNR::Float64 = zero(Float64)
    DIG::Float64 = zero(Float64)
    ELIM::Float64 = zero(Float64)
    FMM::Float64 = zero(Float64)
    FN::Float64 = zero(Float64)
    FNUL::Float64 = zero(Float64)
    HPI::Float64 = zero(Float64)
    I::Int32 = zero(Int32)
    INU::Int32 = zero(Int32)
    INUH::Int32 = zero(Int32)
    IR::Int32 = zero(Int32)
    K::Int32 = zero(Int32)
    K1::Int32 = zero(Int32)
    K2::Int32 = zero(Int32)
    MM::Int32 = zero(Int32)
    MR::Int32 = zero(Int32)
    NN::Int32 = zero(Int32)
    NUF::Int32 = zero(Int32)
    NW::Int32 = zero(Int32)
    R1M5::Float64 = zero(Float64)
    RHPI::Float64 = zero(Float64)
    RL::Float64 = zero(Float64)
    RTOL::Float64 = zero(Float64)
    SGN::Float64 = zero(Float64)
    STI::Float64 = zero(Float64)
    STR::Float64 = zero(Float64)
    TOL::Float64 = zero(Float64)
    UFL::Float64 = zero(Float64)
    ZNI::Float64 = zero(Float64)
    ZNR::Float64 = zero(Float64)
    ZTI::Float64 = zero(Float64)
    begin 
        HPI = 1.5707963267948966
    end
    IERR = Int32(0)
    NZ = Int32(0)
    if ZR == 0.0 && ZI == 0.0
        IERR = Int32(1)
    end
    if FNU < 0.0
        IERR = Int32(1)
    end
    if M < Int32(1) || M > Int32(2)
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
    FN = FNU + DBLE(FLOAT(NN - Int32(1)))
    MM = (Int32(3) - M) - M
    FMM = DBLE(FLOAT(MM))
    ZNR = FMM * ZI
    ZNI = -FMM * ZR
    AZ = abs(COMPLEX(ZR,ZI))
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
        @goto line230
    end
    if FNU > FNUL
        @goto line90
    end
    if FN <= 1.0
        @goto line70
    end
    if FN > 2.0
        @goto line60
    end
    if AZ > TOL
        @goto line70
    end
    ARG = 0.5AZ
    ALN = -FN * DLOG(ARG)
    if ALN > ELIM
        @goto line230
    end
    @goto line70
    @label line60
    (NUF,) = ZUOIK(ZNR,ZNI,FNU,KODE,Int32(2),NN,CYR,CYI,NUF,TOL,ELIM,ALIM)
    if NUF < Int32(0)
        @goto line230
    end
    NZ = NZ + NUF
    NN = NN - NUF
    if NN == Int32(0)
        @goto line140
    end
    @label line70
    if ZNR < 0.0 || ZNR == 0.0 && ZNI < 0.0 && M == Int32(2)
        @goto line80
    end
    (NZ,) = ZBKNU(ZNR,ZNI,FNU,KODE,NN,CYR,CYI,NZ,TOL,ELIM,ALIM)
    @goto line110
    @label line80
    MR = -MM
    (NW,) = ZACON(ZNR,ZNI,FNU,KODE,MR,NN,CYR,CYI,NW,RL,FNUL,TOL,ELIM,ALIM)
    if NW < Int32(0)
        @goto line240
    end
    NZ = NW
    @goto line110
    @label line90
    MR = Int32(0)
    if ZNR >= 0.0 && (ZNR != 0.0 || ZNI >= 0.0 || M != Int32(2))
        @goto line100
    end
    MR = -MM
    if ZNR != 0.0 || ZNI >= 0.0
        @goto line100
    end
    ZNR = -ZNR
    ZNI = -ZNI
    @label line100
    (NW,) = ZBUNK(ZNR,ZNI,FNU,KODE,MR,NN,CYR,CYI,NW,TOL,ELIM,ALIM)
    if NW < Int32(0)
        @goto line240
    end
    NZ = NZ + NW
    @label line110
    SGN = DSIGN(HPI,-FMM)
    INU = INT(SNGL(FNU))
    INUH = div(INU,Int32(2))
    IR = INU - Int32(2) * INUH
    ARG = (FNU - DBLE(FLOAT(INU - IR))) * SGN
    RHPI = 1.0 / SGN
    CSGNI = RHPI * DCOS(ARG)
    CSGNR = -RHPI * DSIN(ARG)
    if MOD(INUH,Int32(2)) == Int32(0)
        @goto line120
    end
    CSGNR = -CSGNR
    CSGNI = -CSGNI
    @label line120
    ZTI = -FMM
    RTOL = 1.0 / TOL
    ASCLE = UFL * RTOL
    for I = Int32(1):NN
        AA = CYR[I]
        BB = CYI[I]
        ATOL = 1.0
        if DMAX1(DABS(AA),DABS(BB)) > ASCLE
            @goto line135
        end
        AA = AA * RTOL
        BB = BB * RTOL
        ATOL = TOL
        @label line135
        STR = AA * CSGNR - BB * CSGNI
        STI = AA * CSGNI + BB * CSGNR
        CYR[I] = STR * ATOL
        CYI[I] = STI * ATOL
        STR = -CSGNI * ZTI
        CSGNI = CSGNR * ZTI
        CSGNR = STR
        @label line130
    end
    return (NZ,IERR)
    @label line140
    if ZNR < 0.0
        @goto line230
    end
    return (NZ,IERR)
    @label line230
    NZ = Int32(0)
    IERR = Int32(2)
    return (NZ,IERR)
    @label line240
    if NW == Int32(-1)
        @goto line230
    end
    NZ = Int32(0)
    IERR = Int32(5)
    return (NZ,IERR)
    @label line260
    NZ = Int32(0)
    IERR = Int32(4)
    return (NZ,IERR)
end
