const _ZUNK2_ZETA2R = Array{Float64}(2)
const _ZUNK2_ZETA2I = Array{Float64}(2)
const _ZUNK2_ZETA1R = Array{Float64}(2)
const _ZUNK2_ZETA1I = Array{Float64}(2)
const _ZUNK2_PHIR = Array{Float64}(2)
const _ZUNK2_PHII = Array{Float64}(2)
const _ZUNK2_CYR = Array{Float64}(2)
const _ZUNK2_CYI = Array{Float64}(2)
const _ZUNK2_CSSR = Array{Float64}(3)
const _ZUNK2_CSRR = Array{Float64}(3)
const _ZUNK2_CIPR = Array{Float64}(4)
const _ZUNK2_CIPI = Array{Float64}(4)
const _ZUNK2_BSUMR = Array{Float64}(2)
const _ZUNK2_BSUMI = Array{Float64}(2)
const _ZUNK2_BRY = Array{Float64}(3)
const _ZUNK2_ASUMR = Array{Float64}(2)
const _ZUNK2_ASUMI = Array{Float64}(2)
const _ZUNK2_ARGR = Array{Float64}(2)
const _ZUNK2_ARGI = Array{Float64}(2)
function ZUNK2(ZR::Float64,ZI::Float64,FNU::Float64,KODE::Int32,MR::Int32,N::Int32,YR::AbstractArray{Float64},YI::AbstractArray{Float64},NZ::Int32,TOL::Float64,ELIM::Float64,ALIM::Float64)
    AARG::Float64 = zero(Float64)
    AIC::Float64 = zero(Float64)
    AII::Float64 = zero(Float64)
    AIR::Float64 = zero(Float64)
    ANG::Float64 = zero(Float64)
    APHI::Float64 = zero(Float64)
    ARGDI::Float64 = zero(Float64)
    ARGDR::Float64 = zero(Float64)
    const ARGI = _ZUNK2_ARGI
    const ARGR = _ZUNK2_ARGR
    ASC::Float64 = zero(Float64)
    ASCLE::Float64 = zero(Float64)
    ASUMDI::Float64 = zero(Float64)
    ASUMDR::Float64 = zero(Float64)
    const ASUMI = _ZUNK2_ASUMI
    const ASUMR = _ZUNK2_ASUMR
    const BRY = _ZUNK2_BRY
    BSUMDI::Float64 = zero(Float64)
    BSUMDR::Float64 = zero(Float64)
    const BSUMI = _ZUNK2_BSUMI
    const BSUMR = _ZUNK2_BSUMR
    C1I::Float64 = zero(Float64)
    C1R::Float64 = zero(Float64)
    C2I::Float64 = zero(Float64)
    C2M::Float64 = zero(Float64)
    C2R::Float64 = zero(Float64)
    CAR::Float64 = zero(Float64)
    const CIPI = _ZUNK2_CIPI
    const CIPR = _ZUNK2_CIPR
    CKI::Float64 = zero(Float64)
    CKR::Float64 = zero(Float64)
    CONER::Float64 = zero(Float64)
    CR1I::Float64 = zero(Float64)
    CR1R::Float64 = zero(Float64)
    CR2I::Float64 = zero(Float64)
    CR2R::Float64 = zero(Float64)
    CRSC::Float64 = zero(Float64)
    CSCL::Float64 = zero(Float64)
    CSGNI::Float64 = zero(Float64)
    CSI::Float64 = zero(Float64)
    CSPNI::Float64 = zero(Float64)
    CSPNR::Float64 = zero(Float64)
    CSR::Float64 = zero(Float64)
    const CSRR = _ZUNK2_CSRR
    const CSSR = _ZUNK2_CSSR
    const CYI = _ZUNK2_CYI
    const CYR = _ZUNK2_CYR
    DAII::Float64 = zero(Float64)
    DAIR::Float64 = zero(Float64)
    FMR::Float64 = zero(Float64)
    FN::Float64 = zero(Float64)
    FNF::Float64 = zero(Float64)
    HPI::Float64 = zero(Float64)
    I::Int32 = zero(Int32)
    IB::Int32 = zero(Int32)
    IC::Int32 = zero(Int32)
    IDUM::Int32 = zero(Int32)
    IFLAG::Int32 = zero(Int32)
    IFN::Int32 = zero(Int32)
    IL::Int32 = zero(Int32)
    IN::Int32 = zero(Int32)
    INU::Int32 = zero(Int32)
    IPARD::Int32 = zero(Int32)
    IUF::Int32 = zero(Int32)
    J::Int32 = zero(Int32)
    K::Int32 = zero(Int32)
    KDFLG::Int32 = zero(Int32)
    KFLAG::Int32 = zero(Int32)
    KK::Int32 = zero(Int32)
    NAI::Int32 = zero(Int32)
    NDAI::Int32 = zero(Int32)
    NW::Int32 = zero(Int32)
    PHIDI::Float64 = zero(Float64)
    PHIDR::Float64 = zero(Float64)
    const PHII = _ZUNK2_PHII
    const PHIR = _ZUNK2_PHIR
    PI::Float64 = zero(Float64)
    PTI::Float64 = zero(Float64)
    PTR::Float64 = zero(Float64)
    RAST::Float64 = zero(Float64)
    RAZR::Float64 = zero(Float64)
    RS1::Float64 = zero(Float64)
    RZI::Float64 = zero(Float64)
    RZR::Float64 = zero(Float64)
    S1I::Float64 = zero(Float64)
    S1R::Float64 = zero(Float64)
    S2I::Float64 = zero(Float64)
    S2R::Float64 = zero(Float64)
    SAR::Float64 = zero(Float64)
    SGN::Float64 = zero(Float64)
    STI::Float64 = zero(Float64)
    STR::Float64 = zero(Float64)
    YY::Float64 = zero(Float64)
    ZBI::Float64 = zero(Float64)
    ZBR::Float64 = zero(Float64)
    ZEROI::Float64 = zero(Float64)
    ZEROR::Float64 = zero(Float64)
    ZET1DI::Float64 = zero(Float64)
    ZET1DR::Float64 = zero(Float64)
    ZET2DI::Float64 = zero(Float64)
    ZET2DR::Float64 = zero(Float64)
    const ZETA1I = _ZUNK2_ZETA1I
    const ZETA1R = _ZUNK2_ZETA1R
    const ZETA2I = _ZUNK2_ZETA2I
    const ZETA2R = _ZUNK2_ZETA2R
    ZNI::Float64 = zero(Float64)
    ZNR::Float64 = zero(Float64)
    ZRI::Float64 = zero(Float64)
    ZRR::Float64 = zero(Float64)
    begin 
        ZEROR = 0.0
        ZEROI = 0.0
        CONER = 1.0
        CR1R = 1.0
        CR1I = 1.7320508075688772
        CR2R = -0.5
        CR2I = -0.8660254037844386
    end
    begin 
        HPI = 1.5707963267948966
        PI = 3.141592653589793
        AIC = 1.2655121234846454
    end
    begin 
        CIPR[Int32(1)] = 1.0
        CIPI[Int32(1)] = 0.0
        CIPR[Int32(2)] = 0.0
        CIPI[Int32(2)] = -1.0
        CIPR[Int32(3)] = -1.0
        CIPI[Int32(3)] = 0.0
        CIPR[Int32(4)] = 0.0
        CIPI[Int32(4)] = 1.0
    end
    KDFLG = Int32(1)
    NZ = Int32(0)
    CSCL = 1.0 / TOL
    CRSC = TOL
    CSSR[Int32(1)] = CSCL
    CSSR[Int32(2)] = CONER
    CSSR[Int32(3)] = CRSC
    CSRR[Int32(1)] = CRSC
    CSRR[Int32(2)] = CONER
    CSRR[Int32(3)] = CSCL
    BRY[Int32(1)] = (1000.0D1MACH1) / TOL
    BRY[Int32(2)] = 1.0 / BRY[Int32(1)]
    BRY[Int32(3)] = D1MACH2
    ZRR = ZR
    ZRI = ZI
    if ZR >= 0.0
        @goto line10
    end
    ZRR = -ZR
    ZRI = -ZI
    @label line10
    YY = ZRI
    ZNR = ZRI
    ZNI = -ZRR
    ZBR = ZRR
    ZBI = ZRI
    INU = trunc(Int32,SNGL(FNU))
    FNF = FNU - DBLE(FLOAT(INU))
    ANG = -HPI * FNF
    CAR = cos(ANG)
    SAR = sin(ANG)
    C2R = HPI * SAR
    C2I = -HPI * CAR
    KK = mod(INU,Int32(4)) + Int32(1)
    STR = C2R * CIPR[KK] - C2I * CIPI[KK]
    STI = C2R * CIPI[KK] + C2I * CIPR[KK]
    CSR = CR1R * STR - CR1I * STI
    CSI = CR1R * STI + CR1I * STR
    if YY > 0.0
        @goto line20
    end
    ZNR = -ZNR
    ZBI = -ZBI
    @label line20
    J = Int32(2)
    for I = Int32(1):N
        J = Int32(3) - J
        FN = FNU + DBLE(FLOAT(I - Int32(1)))
        (PHIR[J],PHII[J],ARGR[J],ARGI[J],ZETA1R[J],ZETA1I[J],ZETA2R[J],ZETA2I[J],ASUMR[J],ASUMI[J],BSUMR[J],BSUMI[J]) = ZUNHJ(ZNR,ZNI,FN,Int32(0),TOL,PHIR[J],PHII[J],ARGR[J],ARGI[J],ZETA1R[J],ZETA1I[J],ZETA2R[J],ZETA2I[J],ASUMR[J],ASUMI[J],BSUMR[J],BSUMI[J])
        if KODE == Int32(1)
            @goto line30
        end
        STR = ZBR + ZETA2R[J]
        STI = ZBI + ZETA2I[J]
        RAST = FN / abs(complex(STR,STI))
        STR = STR * RAST * RAST
        STI = -STI * RAST * RAST
        S1R = ZETA1R[J] - STR
        S1I = ZETA1I[J] - STI
        @goto line40
        @label line30
        S1R = ZETA1R[J] - ZETA2R[J]
        S1I = ZETA1I[J] - ZETA2I[J]
        @label line40
        RS1 = S1R
        if abs(RS1) > ELIM
            @goto line70
        end
        if KDFLG == Int32(1)
            KFLAG = Int32(2)
        end
        if abs(RS1) < ALIM
            @goto line50
        end
        APHI = abs(complex(PHIR[J],PHII[J]))
        AARG = abs(complex(ARGR[J],ARGI[J]))
        RS1 = ((RS1 + log(APHI)) - 0.25 * log(AARG)) - AIC
        if abs(RS1) > ELIM
            @goto line70
        end
        if KDFLG == Int32(1)
            KFLAG = Int32(1)
        end
        if RS1 < 0.0
            @goto line50
        end
        if KDFLG == Int32(1)
            KFLAG = Int32(3)
        end
        @label line50
        C2R = ARGR[J] * CR2R - ARGI[J] * CR2I
        C2I = ARGR[J] * CR2I + ARGI[J] * CR2R
        (AIR,AII,NAI,IDUM) = ZAIRY(C2R,C2I,Int32(0),Int32(2),AIR,AII,NAI,IDUM)
        (DAIR,DAII,NDAI,IDUM) = ZAIRY(C2R,C2I,Int32(1),Int32(2),DAIR,DAII,NDAI,IDUM)
        STR = DAIR * BSUMR[J] - DAII * BSUMI[J]
        STI = DAIR * BSUMI[J] + DAII * BSUMR[J]
        PTR = STR * CR2R - STI * CR2I
        PTI = STR * CR2I + STI * CR2R
        STR = PTR + (AIR * ASUMR[J] - AII * ASUMI[J])
        STI = PTI + (AIR * ASUMI[J] + AII * ASUMR[J])
        PTR = STR * PHIR[J] - STI * PHII[J]
        PTI = STR * PHII[J] + STI * PHIR[J]
        S2R = PTR * CSR - PTI * CSI
        S2I = PTR * CSI + PTI * CSR
        STR = exp(S1R) * CSSR[KFLAG]
        S1R = STR * cos(S1I)
        S1I = STR * sin(S1I)
        STR = S2R * S1R - S2I * S1I
        S2I = S1R * S2I + S2R * S1I
        S2R = STR
        if KFLAG != Int32(1)
            @goto line60
        end
        if zuchk(complex(S2R, S2I), TOL)
            @goto line70
        end
        @label line60
        if YY <= 0.0
            S2I = -S2I
        end
        CYR[KDFLG] = S2R
        CYI[KDFLG] = S2I
        YR[I] = S2R * CSRR[KFLAG]
        YI[I] = S2I * CSRR[KFLAG]
        STR = CSI
        CSI = -CSR
        CSR = STR
        if KDFLG == Int32(2)
            @goto line85
        end
        KDFLG = Int32(2)
        @goto line80
        @label line70
        if RS1 > 0.0
            @goto line320
        end
        if ZR < 0.0
            @goto line320
        end
        KDFLG = Int32(1)
        YR[I] = ZEROR
        YI[I] = ZEROI
        NZ = NZ + Int32(1)
        STR = CSI
        CSI = -CSR
        CSR = STR
        if I == Int32(1)
            @goto line80
        end
        if YR[I - Int32(1)] == ZEROR && YI[I - Int32(1)] == ZEROI
            @goto line80
        end
        YR[I - Int32(1)] = ZEROR
        YI[I - Int32(1)] = ZEROI
        NZ = NZ + Int32(1)
        @label line80
    end
    I = N
    @label line85
    RAZR = 1.0 / abs(complex(ZRR,ZRI))
    STR = ZRR * RAZR
    STI = -ZRI * RAZR
    RZR = (STR + STR) * RAZR
    RZI = (STI + STI) * RAZR
    CKR = FN * RZR
    CKI = FN * RZI
    IB = I + Int32(1)
    if N < IB
        @goto line180
    end
    FN = FNU + DBLE(FLOAT(N - Int32(1)))
    IPARD = Int32(1)
    if MR != Int32(0)
        IPARD = Int32(0)
    end
    (PHIDR,PHIDI,ARGDR,ARGDI,ZET1DR,ZET1DI,ZET2DR,ZET2DI,ASUMDR,ASUMDI,BSUMDR,BSUMDI) = ZUNHJ(ZNR,ZNI,FN,IPARD,TOL,PHIDR,PHIDI,ARGDR,ARGDI,ZET1DR,ZET1DI,ZET2DR,ZET2DI,ASUMDR,ASUMDI,BSUMDR,BSUMDI)
    if KODE == Int32(1)
        @goto line90
    end
    STR = ZBR + ZET2DR
    STI = ZBI + ZET2DI
    RAST = FN / abs(complex(STR,STI))
    STR = STR * RAST * RAST
    STI = -STI * RAST * RAST
    S1R = ZET1DR - STR
    S1I = ZET1DI - STI
    @goto line100
    @label line90
    S1R = ZET1DR - ZET2DR
    S1I = ZET1DI - ZET2DI
    @label line100
    RS1 = S1R
    if abs(RS1) > ELIM
        @goto line105
    end
    if abs(RS1) < ALIM
        @goto line120
    end
    APHI = abs(complex(PHIDR,PHIDI))
    RS1 = RS1 + log(APHI)
    if abs(RS1) < ELIM
        @goto line120
    end
    @label line105
    if RS1 > 0.0
        @goto line320
    end
    if ZR < 0.0
        @goto line320
    end
    NZ = N
    for I = Int32(1):N
        YR[I] = ZEROR
        YI[I] = ZEROI
        @label line106
    end
    return NZ
    @label line120
    S1R = CYR[Int32(1)]
    S1I = CYI[Int32(1)]
    S2R = CYR[Int32(2)]
    S2I = CYI[Int32(2)]
    C1R = CSRR[KFLAG]
    ASCLE = BRY[KFLAG]
    for I = IB:N
        C2R = S2R
        C2I = S2I
        S2R = (CKR * C2R - CKI * C2I) + S1R
        S2I = CKR * C2I + CKI * C2R + S1I
        S1R = C2R
        S1I = C2I
        CKR = CKR + RZR
        CKI = CKI + RZI
        C2R = S2R * C1R
        C2I = S2I * C1R
        YR[I] = C2R
        YI[I] = C2I
        if KFLAG >= Int32(3)
            @goto line130
        end
        STR = abs(C2R)
        STI = abs(C2I)
        C2M = max(STR,STI)
        if C2M <= ASCLE
            @goto line130
        end
        KFLAG = KFLAG + Int32(1)
        ASCLE = BRY[KFLAG]
        S1R = S1R * C1R
        S1I = S1I * C1R
        S2R = C2R
        S2I = C2I
        S1R = S1R * CSSR[KFLAG]
        S1I = S1I * CSSR[KFLAG]
        S2R = S2R * CSSR[KFLAG]
        S2I = S2I * CSSR[KFLAG]
        C1R = CSRR[KFLAG]
        @label line130
    end
    @label line180
    if MR == Int32(0)
        return NZ
    end
    NZ = Int32(0)
    FMR = DBLE(FLOAT(MR))
    SGN = -(copysign(PI,FMR))
    CSGNI = SGN
    if YY <= 0.0
        CSGNI = -CSGNI
    end
    IFN = (INU + N) - Int32(1)
    ANG = FNF * SGN
    CSPNR = cos(ANG)
    CSPNI = sin(ANG)
    if mod(IFN,Int32(2)) == Int32(0)
        @goto line190
    end
    CSPNR = -CSPNR
    CSPNI = -CSPNI
    @label line190
    CSR = SAR * CSGNI
    CSI = CAR * CSGNI
    IN = mod(IFN,Int32(4)) + Int32(1)
    C2R = CIPR[IN]
    C2I = CIPI[IN]
    STR = CSR * C2R + CSI * C2I
    CSI = -CSR * C2I + CSI * C2R
    CSR = STR
    ASC = BRY[Int32(1)]
    IUF = Int32(0)
    KK = N
    KDFLG = Int32(1)
    IB = IB - Int32(1)
    IC = IB - Int32(1)
    for K = Int32(1):N
        FN = FNU + DBLE(FLOAT(KK - Int32(1)))
        if N > Int32(2)
            @goto line175
        end
        @label line172
        PHIDR = PHIR[J]
        PHIDI = PHII[J]
        ARGDR = ARGR[J]
        ARGDI = ARGI[J]
        ZET1DR = ZETA1R[J]
        ZET1DI = ZETA1I[J]
        ZET2DR = ZETA2R[J]
        ZET2DI = ZETA2I[J]
        ASUMDR = ASUMR[J]
        ASUMDI = ASUMI[J]
        BSUMDR = BSUMR[J]
        BSUMDI = BSUMI[J]
        J = Int32(3) - J
        @goto line210
        @label line175
        if KK == N && IB < N
            @goto line210
        end
        if KK == IB || KK == IC
            @goto line172
        end
        (PHIDR,PHIDI,ARGDR,ARGDI,ZET1DR,ZET1DI,ZET2DR,ZET2DI,ASUMDR,ASUMDI,BSUMDR,BSUMDI) = ZUNHJ(ZNR,ZNI,FN,Int32(0),TOL,PHIDR,PHIDI,ARGDR,ARGDI,ZET1DR,ZET1DI,ZET2DR,ZET2DI,ASUMDR,ASUMDI,BSUMDR,BSUMDI)
        @label line210
        if KODE == Int32(1)
            @goto line220
        end
        STR = ZBR + ZET2DR
        STI = ZBI + ZET2DI
        RAST = FN / abs(complex(STR,STI))
        STR = STR * RAST * RAST
        STI = -STI * RAST * RAST
        S1R = -ZET1DR + STR
        S1I = -ZET1DI + STI
        @goto line230
        @label line220
        S1R = -ZET1DR + ZET2DR
        S1I = -ZET1DI + ZET2DI
        @label line230
        RS1 = S1R
        if abs(RS1) > ELIM
            @goto line280
        end
        if KDFLG == Int32(1)
            IFLAG = Int32(2)
        end
        if abs(RS1) < ALIM
            @goto line240
        end
        APHI = abs(complex(PHIDR,PHIDI))
        AARG = abs(complex(ARGDR,ARGDI))
        RS1 = ((RS1 + log(APHI)) - 0.25 * log(AARG)) - AIC
        if abs(RS1) > ELIM
            @goto line280
        end
        if KDFLG == Int32(1)
            IFLAG = Int32(1)
        end
        if RS1 < 0.0
            @goto line240
        end
        if KDFLG == Int32(1)
            IFLAG = Int32(3)
        end
        @label line240
        (AIR,AII,NAI,IDUM) = ZAIRY(ARGDR,ARGDI,Int32(0),Int32(2),AIR,AII,NAI,IDUM)
        (DAIR,DAII,NDAI,IDUM) = ZAIRY(ARGDR,ARGDI,Int32(1),Int32(2),DAIR,DAII,NDAI,IDUM)
        STR = DAIR * BSUMDR - DAII * BSUMDI
        STI = DAIR * BSUMDI + DAII * BSUMDR
        STR = STR + (AIR * ASUMDR - AII * ASUMDI)
        STI = STI + (AIR * ASUMDI + AII * ASUMDR)
        PTR = STR * PHIDR - STI * PHIDI
        PTI = STR * PHIDI + STI * PHIDR
        S2R = PTR * CSR - PTI * CSI
        S2I = PTR * CSI + PTI * CSR
        STR = exp(S1R) * CSSR[IFLAG]
        S1R = STR * cos(S1I)
        S1I = STR * sin(S1I)
        STR = S2R * S1R - S2I * S1I
        S2I = S2R * S1I + S2I * S1R
        S2R = STR
        if IFLAG != Int32(1)
            @goto line250
        end
        if zuchk(complex(S2R, S2I), TOL)
            @goto line250
        end
        S2R = ZEROR
        S2I = ZEROI
        @label line250
        if YY <= 0.0
            S2I = -S2I
        end
        CYR[KDFLG] = S2R
        CYI[KDFLG] = S2I
        C2R = S2R
        C2I = S2I
        S2R = S2R * CSRR[IFLAG]
        S2I = S2I * CSRR[IFLAG]
        S1R = YR[KK]
        S1I = YI[KK]
        if KODE == Int32(1)
            @goto line270
        end
        (S1R,S1I,S2R,S2I,NW,IUF) = ZS1S2(ZRR,ZRI,S1R,S1I,S2R,S2I,NW,ASC,ALIM,IUF)
        NZ = NZ + NW
        @label line270
        YR[KK] = (S1R * CSPNR - S1I * CSPNI) + S2R
        YI[KK] = S1R * CSPNI + S1I * CSPNR + S2I
        KK = KK - Int32(1)
        CSPNR = -CSPNR
        CSPNI = -CSPNI
        STR = CSI
        CSI = -CSR
        CSR = STR
        if C2R != 0.0 || C2I != 0.0
            @goto line255
        end
        KDFLG = Int32(1)
        @goto line290
        @label line255
        if KDFLG == Int32(2)
            @goto line295
        end
        KDFLG = Int32(2)
        @goto line290
        @label line280
        if RS1 > 0.0
            @goto line320
        end
        S2R = ZEROR
        S2I = ZEROI
        @goto line250
        @label line290
    end
    K = N
    @label line295
    IL = N - K
    if IL == Int32(0)
        return NZ
    end
    S1R = CYR[Int32(1)]
    S1I = CYI[Int32(1)]
    S2R = CYR[Int32(2)]
    S2I = CYI[Int32(2)]
    CSR = CSRR[IFLAG]
    ASCLE = BRY[IFLAG]
    FN = DBLE(FLOAT(INU + IL))
    for I = Int32(1):IL
        C2R = S2R
        C2I = S2I
        S2R = S1R + (FN + FNF) * (RZR * C2R - RZI * C2I)
        S2I = S1I + (FN + FNF) * (RZR * C2I + RZI * C2R)
        S1R = C2R
        S1I = C2I
        FN = FN - 1.0
        C2R = S2R * CSR
        C2I = S2I * CSR
        CKR = C2R
        CKI = C2I
        C1R = YR[KK]
        C1I = YI[KK]
        if KODE == Int32(1)
            @goto line300
        end
        (C1R,C1I,C2R,C2I,NW,IUF) = ZS1S2(ZRR,ZRI,C1R,C1I,C2R,C2I,NW,ASC,ALIM,IUF)
        NZ = NZ + NW
        @label line300
        YR[KK] = (C1R * CSPNR - C1I * CSPNI) + C2R
        YI[KK] = C1R * CSPNI + C1I * CSPNR + C2I
        KK = KK - Int32(1)
        CSPNR = -CSPNR
        CSPNI = -CSPNI
        if IFLAG >= Int32(3)
            @goto line310
        end
        C2R = abs(CKR)
        C2I = abs(CKI)
        C2M = max(C2R,C2I)
        if C2M <= ASCLE
            @goto line310
        end
        IFLAG = IFLAG + Int32(1)
        ASCLE = BRY[IFLAG]
        S1R = S1R * CSR
        S1I = S1I * CSR
        S2R = CKR
        S2I = CKI
        S1R = S1R * CSSR[IFLAG]
        S1I = S1I * CSSR[IFLAG]
        S2R = S2R * CSSR[IFLAG]
        S2I = S2I * CSSR[IFLAG]
        CSR = CSRR[IFLAG]
        @label line310
    end
    return NZ
    @label line320
    NZ = Int32(-1)
    return NZ
end
