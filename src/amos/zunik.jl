const _ZUNIK_CON = Array(Float64,2)
const _ZUNIK_C = Array(Float64,120)
function ZUNIK(ZRR::Float64,ZRI::Float64,FNU::Float64,IKFLG::Int32,IPMTR::Int32,TOL::Float64,INIT::Int32,PHIR::Float64,PHII::Float64,ZETA1R::Float64,ZETA1I::Float64,ZETA2R::Float64,ZETA2I::Float64,SUMR::Float64,SUMI::Float64,CWRKR::AbstractArray{Float64},CWRKI::AbstractArray{Float64})
    AC::Float64 = zero(Float64)
    const C = _ZUNIK_C
    const CON = _ZUNIK_CON
    CONEI::Float64 = zero(Float64)
    CONER::Float64 = zero(Float64)
    CRFNI::Float64 = zero(Float64)
    CRFNR::Float64 = zero(Float64)
    I::Int32 = zero(Int32)
    IDUM::Int32 = zero(Int32)
    J::Int32 = zero(Int32)
    K::Int32 = zero(Int32)
    L::Int32 = zero(Int32)
    RFN::Float64 = zero(Float64)
    SI::Float64 = zero(Float64)
    SR::Float64 = zero(Float64)
    SRI::Float64 = zero(Float64)
    SRR::Float64 = zero(Float64)
    STI::Float64 = zero(Float64)
    STR::Float64 = zero(Float64)
    T2I::Float64 = zero(Float64)
    T2R::Float64 = zero(Float64)
    TEST::Float64 = zero(Float64)
    TI::Float64 = zero(Float64)
    TR::Float64 = zero(Float64)
    ZEROI::Float64 = zero(Float64)
    ZEROR::Float64 = zero(Float64)
    ZNI::Float64 = zero(Float64)
    ZNR::Float64 = zero(Float64)
    begin 
        ZEROR = 0.0
        ZEROI = 0.0
        CONER = 1.0
        CONEI = 0.0
    end
    begin 
        CON[int32(1)] = 0.3989422804014327
        CON[int32(2)] = 1.2533141373155003
    end
    begin 
        C[int32(1)] = 1.0
        C[int32(2)] = -0.20833333333333334
        C[int32(3)] = 0.125
        C[int32(4)] = 0.3342013888888889
        C[int32(5)] = -0.4010416666666667
        C[int32(6)] = 0.0703125
        C[int32(7)] = -1.0258125964506173
        C[int32(8)] = 1.8464626736111112
        C[int32(9)] = -0.8912109375
        C[int32(10)] = 0.0732421875
        C[int32(11)] = 4.669584423426247
        C[int32(12)] = -11.207002616222994
        C[int32(13)] = 8.78912353515625
        C[int32(14)] = -2.3640869140625
        C[int32(15)] = 0.112152099609375
        C[int32(16)] = -28.212072558200244
        C[int32(17)] = 84.63621767460073
        C[int32(18)] = -91.81824154324002
        C[int32(19)] = 42.53499874538846
        C[int32(20)] = -7.368794359479632
        C[int32(21)] = 0.22710800170898438
        C[int32(22)] = 212.57013003921713
        C[int32(23)] = -765.2524681411817
        C[int32(24)] = 1059.9904525279999
    end
    begin 
        C[int32(25)] = -699.5796273761325
        C[int32(26)] = 218.1905117442116
        C[int32(27)] = -26.491430486951554
        C[int32(28)] = 0.5725014209747314
        C[int32(29)] = -1919.457662318407
        C[int32(30)] = 8061.722181737309
        C[int32(31)] = -13586.550006434138
        C[int32(32)] = 11655.393336864534
        C[int32(33)] = -5305.646978613403
        C[int32(34)] = 1200.9029132163525
        C[int32(35)] = -108.09091978839466
        C[int32(36)] = 1.7277275025844574
        C[int32(37)] = 20204.29133096615
        C[int32(38)] = -96980.59838863752
        C[int32(39)] = 192547.00123253153
        C[int32(40)] = -203400.17728041555
        C[int32(41)] = 122200.46498301746
        C[int32(42)] = -41192.65496889755
        C[int32(43)] = 7109.514302489364
        C[int32(44)] = -493.915304773088
        C[int32(45)] = 6.074042001273483
        C[int32(46)] = -242919.18790055133
        C[int32(47)] = 1.3117636146629772e6
        C[int32(48)] = -2.9980159185381066e6
    end
    begin 
        C[int32(49)] = 3.763271297656404e6
        C[int32(50)] = -2.813563226586534e6
        C[int32(51)] = 1.2683652733216248e6
        C[int32(52)] = -331645.1724845636
        C[int32(53)] = 45218.76898136273
        C[int32(54)] = -2499.8304818112097
        C[int32(55)] = 24.380529699556064
        C[int32(56)] = 3.284469853072038e6
        C[int32(57)] = -1.9706819118432228e7
        C[int32(58)] = 5.095260249266464e7
        C[int32(59)] = -7.410514821153265e7
        C[int32(60)] = 6.634451227472903e7
        C[int32(61)] = -3.756717666076335e7
        C[int32(62)] = 1.3288767166421818e7
        C[int32(63)] = -2.7856181280864547e6
        C[int32(64)] = 308186.4046126624
        C[int32(65)] = -13886.08975371704
        C[int32(66)] = 110.01714026924674
        C[int32(67)] = -4.932925366450996e7
        C[int32(68)] = 3.2557307418576574e8
        C[int32(69)] = -9.394623596815784e8
        C[int32(70)] = 1.55359689957058e9
        C[int32(71)] = -1.6210805521083372e9
        C[int32(72)] = 1.1068428168230145e9
    end
    begin 
        C[int32(73)] = -4.958897842750303e8
        C[int32(74)] = 1.420629077975331e8
        C[int32(75)] = -2.447406272573873e7
        C[int32(76)] = 2.2437681779224495e6
        C[int32(77)] = -84005.43360302408
        C[int32(78)] = 551.3358961220206
        C[int32(79)] = 8.147890961183121e8
        C[int32(80)] = -5.866481492051847e9
        C[int32(81)] = 1.8688207509295826e10
        C[int32(82)] = -3.4632043388158775e10
        C[int32(83)] = 4.1280185579753975e10
        C[int32(84)] = -3.3026599749800724e10
        C[int32(85)] = 1.79542137311556e10
        C[int32(86)] = -6.563293792619285e9
        C[int32(87)] = 1.5592798648792574e9
        C[int32(88)] = -2.2510566188941526e8
        C[int32(89)] = 1.7395107553978164e7
        C[int32(90)] = -549842.3275722887
        C[int32(91)] = 3038.090510922384
        C[int32(92)] = -1.4679261247695616e10
        C[int32(93)] = 1.144982377320258e11
        C[int32(94)] = -3.990961752244665e11
        C[int32(95)] = 8.192186695485773e11
        C[int32(96)] = -1.0983751560812233e12
    end
    begin 
        C[int32(97)] = 1.0081581068653821e12
        C[int32(98)] = -6.453648692453765e11
        C[int32(99)] = 2.879006499061506e11
        C[int32(100)] = -8.786707217802327e10
        C[int32(101)] = 1.763473060683497e10
        C[int32(102)] = -2.167164983223795e9
        C[int32(103)] = 1.4315787671888897e8
        C[int32(104)] = -3.871833442572613e6
        C[int32(105)] = 18257.755474293175
        C[int32(106)] = 2.86464035717679e11
        C[int32(107)] = -2.406297900028504e12
        C[int32(108)] = 9.109341185239898e12
        C[int32(109)] = -2.0516899410934438e13
        C[int32(110)] = 3.056512551993532e13
        C[int32(111)] = -3.166708858478516e13
        C[int32(112)] = 2.334836404458184e13
        C[int32(113)] = -1.2320491305598287e13
        C[int32(114)] = 4.612725780849132e12
        C[int32(115)] = -1.1965528801961816e12
        C[int32(116)] = 2.0591450323241e11
        C[int32(117)] = -2.1822927757529224e10
        C[int32(118)] = 1.2470092935127103e9
    end
    begin 
        C[int32(119)] = -2.9188388122220814e7
        C[int32(120)] = 118838.42625678325
    end
    if INIT != int32(0)
        @goto line40
    end
    RFN = 1.0 / FNU
    TEST = D1MACH1 * 1000.0
    AC = FNU * TEST
    if DABS(ZRR) > AC || DABS(ZRI) > AC
        @goto line15
    end
    ZETA1R = 2.0 * DABS(DLOG(TEST)) + FNU
    ZETA1I = 0.0
    ZETA2R = FNU
    ZETA2I = 0.0
    PHIR = 1.0
    PHII = 0.0
    return (INIT,PHIR,PHII,ZETA1R,ZETA1I,ZETA2R,ZETA2I,SUMR,SUMI)
    @label line15
    TR = ZRR * RFN
    TI = ZRI * RFN
    SR = CONER + (TR * TR - TI * TI)
    SI = CONEI + (TR * TI + TI * TR)
    (SRR,SRI) = ZSQRT(SR,SI,SRR,SRI)
    STR = CONER + SRR
    STI = CONEI + SRI
    (ZNR,ZNI) = ZDIV(STR,STI,TR,TI,ZNR,ZNI)
    (STR,STI,IDUM) = ZLOG(ZNR,ZNI,STR,STI,IDUM)
    ZETA1R = FNU * STR
    ZETA1I = FNU * STI
    ZETA2R = FNU * SRR
    ZETA2I = FNU * SRI
    (TR,TI) = ZDIV(CONER,CONEI,SRR,SRI,TR,TI)
    SRR = TR * RFN
    SRI = TI * RFN
    (CWRKR[int32(16)],CWRKI[int32(16)]) = ZSQRT(SRR,SRI,CWRKR[int32(16)],CWRKI[int32(16)])
    PHIR = CWRKR[int32(16)] * CON[IKFLG]
    PHII = CWRKI[int32(16)] * CON[IKFLG]
    if IPMTR != int32(0)
        return (INIT,PHIR,PHII,ZETA1R,ZETA1I,ZETA2R,ZETA2I,SUMR,SUMI)
    end
    (T2R,T2I) = ZDIV(CONER,CONEI,SR,SI,T2R,T2I)
    CWRKR[int32(1)] = CONER
    CWRKI[int32(1)] = CONEI
    CRFNR = CONER
    CRFNI = CONEI
    AC = 1.0
    L = int32(1)
    for K = int32(2):int32(15)
        SR = ZEROR
        SI = ZEROI
        for J = int32(1):K
            L = L + int32(1)
            STR = (SR * T2R - SI * T2I) + C[L]
            SI = SR * T2I + SI * T2R
            SR = STR
            @label line10
        end
        STR = CRFNR * SRR - CRFNI * SRI
        CRFNI = CRFNR * SRI + CRFNI * SRR
        CRFNR = STR
        CWRKR[K] = CRFNR * SR - CRFNI * SI
        CWRKI[K] = CRFNR * SI + CRFNI * SR
        AC = AC * RFN
        TEST = DABS(CWRKR[K]) + DABS(CWRKI[K])
        if AC < TOL && TEST < TOL
            @goto line30
        end
        @label line20
    end
    K = int32(15)
    @label line30
    INIT = K
    @label line40
    if IKFLG == int32(2)
        @goto line60
    end
    SR = ZEROR
    SI = ZEROI
    for I = int32(1):INIT
        SR = SR + CWRKR[I]
        SI = SI + CWRKI[I]
        @label line50
    end
    SUMR = SR
    SUMI = SI
    PHIR = CWRKR[int32(16)] * CON[int32(1)]
    PHII = CWRKI[int32(16)] * CON[int32(1)]
    return (INIT,PHIR,PHII,ZETA1R,ZETA1I,ZETA2R,ZETA2I,SUMR,SUMI)
    @label line60
    SR = ZEROR
    SI = ZEROI
    TR = CONER
    for I = int32(1):INIT
        SR = SR + TR * CWRKR[I]
        SI = SI + TR * CWRKI[I]
        TR = -TR
        @label line70
    end
    SUMR = SR
    SUMI = SI
    PHIR = CWRKR[int32(16)] * CON[int32(2)]
    PHII = CWRKI[int32(16)] * CON[int32(2)]
    return (INIT,PHIR,PHII,ZETA1R,ZETA1I,ZETA2R,ZETA2I,SUMR,SUMI)
end
