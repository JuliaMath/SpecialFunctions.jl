const _ZUNIK_CON = Array{Float64}(2)
const _ZUNIK_C = Array{Float64}(120)
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
        CON[Int32(1)] = 0.3989422804014327
        CON[Int32(2)] = 1.2533141373155003
    end
    begin 
        C[Int32(1)] = 1.0
        C[Int32(2)] = -0.20833333333333334
        C[Int32(3)] = 0.125
        C[Int32(4)] = 0.3342013888888889
        C[Int32(5)] = -0.4010416666666667
        C[Int32(6)] = 0.0703125
        C[Int32(7)] = -1.0258125964506173
        C[Int32(8)] = 1.8464626736111112
        C[Int32(9)] = -0.8912109375
        C[Int32(10)] = 0.0732421875
        C[Int32(11)] = 4.669584423426247
        C[Int32(12)] = -11.207002616222994
        C[Int32(13)] = 8.78912353515625
        C[Int32(14)] = -2.3640869140625
        C[Int32(15)] = 0.112152099609375
        C[Int32(16)] = -28.212072558200244
        C[Int32(17)] = 84.63621767460073
        C[Int32(18)] = -91.81824154324002
        C[Int32(19)] = 42.53499874538846
        C[Int32(20)] = -7.368794359479632
        C[Int32(21)] = 0.22710800170898438
        C[Int32(22)] = 212.57013003921713
        C[Int32(23)] = -765.2524681411817
        C[Int32(24)] = 1059.9904525279999
    end
    begin 
        C[Int32(25)] = -699.5796273761325
        C[Int32(26)] = 218.1905117442116
        C[Int32(27)] = -26.491430486951554
        C[Int32(28)] = 0.5725014209747314
        C[Int32(29)] = -1919.457662318407
        C[Int32(30)] = 8061.722181737309
        C[Int32(31)] = -13586.550006434138
        C[Int32(32)] = 11655.393336864534
        C[Int32(33)] = -5305.646978613403
        C[Int32(34)] = 1200.9029132163525
        C[Int32(35)] = -108.09091978839466
        C[Int32(36)] = 1.7277275025844574
        C[Int32(37)] = 20204.29133096615
        C[Int32(38)] = -96980.59838863752
        C[Int32(39)] = 192547.00123253153
        C[Int32(40)] = -203400.17728041555
        C[Int32(41)] = 122200.46498301746
        C[Int32(42)] = -41192.65496889755
        C[Int32(43)] = 7109.514302489364
        C[Int32(44)] = -493.915304773088
        C[Int32(45)] = 6.074042001273483
        C[Int32(46)] = -242919.18790055133
        C[Int32(47)] = 1.3117636146629772e6
        C[Int32(48)] = -2.9980159185381066e6
    end
    begin 
        C[Int32(49)] = 3.763271297656404e6
        C[Int32(50)] = -2.813563226586534e6
        C[Int32(51)] = 1.2683652733216248e6
        C[Int32(52)] = -331645.1724845636
        C[Int32(53)] = 45218.76898136273
        C[Int32(54)] = -2499.8304818112097
        C[Int32(55)] = 24.380529699556064
        C[Int32(56)] = 3.284469853072038e6
        C[Int32(57)] = -1.9706819118432228e7
        C[Int32(58)] = 5.095260249266464e7
        C[Int32(59)] = -7.410514821153265e7
        C[Int32(60)] = 6.634451227472903e7
        C[Int32(61)] = -3.756717666076335e7
        C[Int32(62)] = 1.3288767166421818e7
        C[Int32(63)] = -2.7856181280864547e6
        C[Int32(64)] = 308186.4046126624
        C[Int32(65)] = -13886.08975371704
        C[Int32(66)] = 110.01714026924674
        C[Int32(67)] = -4.932925366450996e7
        C[Int32(68)] = 3.2557307418576574e8
        C[Int32(69)] = -9.394623596815784e8
        C[Int32(70)] = 1.55359689957058e9
        C[Int32(71)] = -1.6210805521083372e9
        C[Int32(72)] = 1.1068428168230145e9
    end
    begin 
        C[Int32(73)] = -4.958897842750303e8
        C[Int32(74)] = 1.420629077975331e8
        C[Int32(75)] = -2.447406272573873e7
        C[Int32(76)] = 2.2437681779224495e6
        C[Int32(77)] = -84005.43360302408
        C[Int32(78)] = 551.3358961220206
        C[Int32(79)] = 8.147890961183121e8
        C[Int32(80)] = -5.866481492051847e9
        C[Int32(81)] = 1.8688207509295826e10
        C[Int32(82)] = -3.4632043388158775e10
        C[Int32(83)] = 4.1280185579753975e10
        C[Int32(84)] = -3.3026599749800724e10
        C[Int32(85)] = 1.79542137311556e10
        C[Int32(86)] = -6.563293792619285e9
        C[Int32(87)] = 1.5592798648792574e9
        C[Int32(88)] = -2.2510566188941526e8
        C[Int32(89)] = 1.7395107553978164e7
        C[Int32(90)] = -549842.3275722887
        C[Int32(91)] = 3038.090510922384
        C[Int32(92)] = -1.4679261247695616e10
        C[Int32(93)] = 1.144982377320258e11
        C[Int32(94)] = -3.990961752244665e11
        C[Int32(95)] = 8.192186695485773e11
        C[Int32(96)] = -1.0983751560812233e12
    end
    begin 
        C[Int32(97)] = 1.0081581068653821e12
        C[Int32(98)] = -6.453648692453765e11
        C[Int32(99)] = 2.879006499061506e11
        C[Int32(100)] = -8.786707217802327e10
        C[Int32(101)] = 1.763473060683497e10
        C[Int32(102)] = -2.167164983223795e9
        C[Int32(103)] = 1.4315787671888897e8
        C[Int32(104)] = -3.871833442572613e6
        C[Int32(105)] = 18257.755474293175
        C[Int32(106)] = 2.86464035717679e11
        C[Int32(107)] = -2.406297900028504e12
        C[Int32(108)] = 9.109341185239898e12
        C[Int32(109)] = -2.0516899410934438e13
        C[Int32(110)] = 3.056512551993532e13
        C[Int32(111)] = -3.166708858478516e13
        C[Int32(112)] = 2.334836404458184e13
        C[Int32(113)] = -1.2320491305598287e13
        C[Int32(114)] = 4.612725780849132e12
        C[Int32(115)] = -1.1965528801961816e12
        C[Int32(116)] = 2.0591450323241e11
        C[Int32(117)] = -2.1822927757529224e10
        C[Int32(118)] = 1.2470092935127103e9
    end
    begin 
        C[Int32(119)] = -2.9188388122220814e7
        C[Int32(120)] = 118838.42625678325
    end
    if INIT != Int32(0)
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
    (SRR,SRI) = reim(sqrt(complex(SR,SI)))
    STR = CONER + SRR
    STI = CONEI + SRI
    (ZNR,ZNI) = reim(complex(STR,STI) / complex(TR,TI))
    (STR,STI) = reim(log(complex(ZNR,ZNI)))
    ZETA1R = FNU * STR
    ZETA1I = FNU * STI
    ZETA2R = FNU * SRR
    ZETA2I = FNU * SRI
    (TR,TI) = reim(complex(CONER,CONEI) / complex(SRR,SRI))
    SRR = TR * RFN
    SRI = TI * RFN
    (CWRKR[Int32(16)],CWRKI[Int32(16)]) = reim(sqrt(complex(SRR,SRI)))
    PHIR = CWRKR[Int32(16)] * CON[IKFLG]
    PHII = CWRKI[Int32(16)] * CON[IKFLG]
    if IPMTR != Int32(0)
        return (INIT,PHIR,PHII,ZETA1R,ZETA1I,ZETA2R,ZETA2I,SUMR,SUMI)
    end
    (T2R,T2I) = reim(complex(CONER,CONEI) / complex(SR,SI))
    CWRKR[Int32(1)] = CONER
    CWRKI[Int32(1)] = CONEI
    CRFNR = CONER
    CRFNI = CONEI
    AC = 1.0
    L = Int32(1)
    for K = Int32(2):Int32(15)
        SR = ZEROR
        SI = ZEROI
        for J = Int32(1):K
            L = L + Int32(1)
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
    K = Int32(15)
    @label line30
    INIT = K
    @label line40
    if IKFLG == Int32(2)
        @goto line60
    end
    SR = ZEROR
    SI = ZEROI
    for I = Int32(1):INIT
        SR = SR + CWRKR[I]
        SI = SI + CWRKI[I]
        @label line50
    end
    SUMR = SR
    SUMI = SI
    PHIR = CWRKR[Int32(16)] * CON[Int32(1)]
    PHII = CWRKI[Int32(16)] * CON[Int32(1)]
    return (INIT,PHIR,PHII,ZETA1R,ZETA1I,ZETA2R,ZETA2I,SUMR,SUMI)
    @label line60
    SR = ZEROR
    SI = ZEROI
    TR = CONER
    for I = Int32(1):INIT
        SR = SR + TR * CWRKR[I]
        SI = SI + TR * CWRKI[I]
        TR = -TR
        @label line70
    end
    SUMR = SR
    SUMI = SI
    PHIR = CWRKR[Int32(16)] * CON[Int32(2)]
    PHII = CWRKI[Int32(16)] * CON[Int32(2)]
    return (INIT,PHIR,PHII,ZETA1R,ZETA1I,ZETA2R,ZETA2I,SUMR,SUMI)
end
