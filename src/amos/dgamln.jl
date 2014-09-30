const _DGAMLN_GLN = Array(Float64,100)
const _DGAMLN_CF = Array(Float64,22)
function DGAMLN(Z::Float64,IERR::Int32)
    __DGAMLN__::Float64 = zero(Float64)
    const CF = _DGAMLN_CF
    CON::Float64 = zero(Float64)
    FLN::Float64 = zero(Float64)
    FZ::Float64 = zero(Float64)
    const GLN = _DGAMLN_GLN
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
    begin 
        GLN[int32(1)] = 0.0
        GLN[int32(2)] = 0.0
        GLN[int32(3)] = 0.6931471805599453
        GLN[int32(4)] = 1.791759469228055
        GLN[int32(5)] = 3.1780538303479458
        GLN[int32(6)] = 4.787491742782046
        GLN[int32(7)] = 6.579251212010101
        GLN[int32(8)] = 8.525161361065415
        GLN[int32(9)] = 10.60460290274525
        GLN[int32(10)] = 12.801827480081469
        GLN[int32(11)] = 15.104412573075516
        GLN[int32(12)] = 17.502307845873887
        GLN[int32(13)] = 19.987214495661885
        GLN[int32(14)] = 22.552163853123425
        GLN[int32(15)] = 25.19122118273868
        GLN[int32(16)] = 27.89927138384089
        GLN[int32(17)] = 30.671860106080672
        GLN[int32(18)] = 33.50507345013689
        GLN[int32(19)] = 36.39544520803305
        GLN[int32(20)] = 39.339884187199495
        GLN[int32(21)] = 42.335616460753485
        GLN[int32(22)] = 45.38013889847691
    end
    begin 
        GLN[int32(23)] = 48.47118135183523
        GLN[int32(24)] = 51.60667556776438
        GLN[int32(25)] = 54.78472939811232
        GLN[int32(26)] = 58.00360522298052
        GLN[int32(27)] = 61.261701761002
        GLN[int32(28)] = 64.55753862700634
        GLN[int32(29)] = 67.88974313718154
        GLN[int32(30)] = 71.25703896716801
        GLN[int32(31)] = 74.65823634883016
        GLN[int32(32)] = 78.0922235533153
        GLN[int32(33)] = 81.55795945611504
        GLN[int32(34)] = 85.05446701758152
        GLN[int32(35)] = 88.58082754219768
        GLN[int32(36)] = 92.1361756036871
        GLN[int32(37)] = 95.7196945421432
        GLN[int32(38)] = 99.33061245478743
        GLN[int32(39)] = 102.96819861451381
        GLN[int32(40)] = 106.63176026064346
        GLN[int32(41)] = 110.32063971475739
        GLN[int32(42)] = 114.0342117814617
        GLN[int32(43)] = 117.77188139974507
        GLN[int32(44)] = 121.53308151543864
    end
    begin 
        GLN[int32(45)] = 125.3172711493569
        GLN[int32(46)] = 129.12393363912722
        GLN[int32(47)] = 132.95257503561632
        GLN[int32(48)] = 136.80272263732635
        GLN[int32(49)] = 140.67392364823425
        GLN[int32(50)] = 144.5657439463449
        GLN[int32(51)] = 148.47776695177302
        GLN[int32(52)] = 152.40959258449735
        GLN[int32(53)] = 156.3608363030788
        GLN[int32(54)] = 160.3311282166309
        GLN[int32(55)] = 164.32011226319517
        GLN[int32(56)] = 168.32744544842765
        GLN[int32(57)] = 172.3527971391628
        GLN[int32(58)] = 176.39584840699735
        GLN[int32(59)] = 180.45629141754378
        GLN[int32(60)] = 184.53382886144948
        GLN[int32(61)] = 188.6281734236716
        GLN[int32(62)] = 192.7390472878449
        GLN[int32(63)] = 196.86618167289
        GLN[int32(64)] = 201.00931639928152
        GLN[int32(65)] = 205.1681994826412
        GLN[int32(66)] = 209.34258675253685
    end
    begin 
        GLN[int32(67)] = 213.53224149456327
        GLN[int32(68)] = 217.73693411395422
        GLN[int32(69)] = 221.95644181913033
        GLN[int32(70)] = 226.1905483237276
        GLN[int32(71)] = 230.43904356577696
        GLN[int32(72)] = 234.70172344281826
        GLN[int32(73)] = 238.97838956183432
        GLN[int32(74)] = 243.2688490029827
        GLN[int32(75)] = 247.57291409618688
        GLN[int32(76)] = 251.8904022097232
        GLN[int32(77)] = 256.22113555000954
        GLN[int32(78)] = 260.5649409718632
        GLN[int32(79)] = 264.9216497985528
        GLN[int32(80)] = 269.2910976510198
        GLN[int32(81)] = 273.6731242856937
        GLN[int32(82)] = 278.0675734403661
        GLN[int32(83)] = 282.4742926876304
        GLN[int32(84)] = 286.893133295427
        GLN[int32(85)] = 291.3239500942703
        GLN[int32(86)] = 295.76660135076065
        GLN[int32(87)] = 300.22094864701415
        GLN[int32(88)] = 304.6868567656687
    end
    begin 
        GLN[int32(89)] = 309.1641935801469
        GLN[int32(90)] = 313.65282994987905
        GLN[int32(91)] = 318.1526396202093
        GLN[int32(92)] = 322.66349912672615
        GLN[int32(93)] = 327.1852877037752
        GLN[int32(94)] = 331.7178871969285
        GLN[int32(95)] = 336.26118197919845
        GLN[int32(96)] = 340.815058870799
        GLN[int32(97)] = 345.37940706226686
        GLN[int32(98)] = 349.95411804077025
        GLN[int32(99)] = 354.5390855194408
        GLN[int32(100)] = 359.1342053695754
    end
    begin 
        CF[int32(1)] = 0.08333333333333333
        CF[int32(2)] = -0.002777777777777778
        CF[int32(3)] = 0.0007936507936507937
        CF[int32(4)] = -0.0005952380952380953
        CF[int32(5)] = 0.0008417508417508417
        CF[int32(6)] = -0.0019175269175269176
        CF[int32(7)] = 0.00641025641025641
        CF[int32(8)] = -0.029550653594771242
        CF[int32(9)] = 0.17964437236883057
        CF[int32(10)] = -1.3924322169059011
        CF[int32(11)] = 13.402864044168393
        CF[int32(12)] = -156.84828462600203
        CF[int32(13)] = 2193.1033333333335
        CF[int32(14)] = -36108.77125372499
        CF[int32(15)] = 691472.268851313
        CF[int32(16)] = -1.5238221539407415e7
        CF[int32(17)] = 3.8290075139141417e8
        CF[int32(18)] = -1.0882266035784391e10
        CF[int32(19)] = 3.4732028376500226e11
        CF[int32(20)] = -1.2369602142269275e13
        CF[int32(21)] = 4.887880647930793e14
        CF[int32(22)] = -2.1320333960919372e16
    end
    begin 
        CON = 1.8378770664093456
    end
    IERR = int32(0)
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
    if NZ > int32(100)
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
    MZ = INT(SNGL(ZM)) + int32(1)
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
    T1 = CF[int32(1)] * ZP
    S = T1
    if ZP < WDTOL
        @goto line40
    end
    ZSQ = ZP * ZP
    TST = T1 * WDTOL
    for K = int32(2):int32(22)
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
    for I = int32(1):NZ
        ZP = ZP * (Z + FLOAT(I - int32(1)))
        @label line60
    end
    TLG = DLOG(ZDMY)
    __DGAMLN__ = (ZDMY * (TLG - 1.0) - DLOG(ZP)) + 0.5 * (CON - TLG) + S
    return __DGAMLN__
    @label line70
    IERR = int32(1)
    return __DGAMLN__
end
