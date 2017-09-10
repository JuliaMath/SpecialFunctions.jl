const _ZUNHJ_UPR = Array{Float64}(14)
const _ZUNHJ_UPI = Array{Float64}(14)
const _ZUNHJ_PR = Array{Float64}(30)
const _ZUNHJ_PI = Array{Float64}(30)
const _ZUNHJ_GAMA = Array{Float64}(30)
const _ZUNHJ_DRR = Array{Float64}(14)
const _ZUNHJ_DRI = Array{Float64}(14)
const _ZUNHJ_CRR = Array{Float64}(14)
const _ZUNHJ_CRI = Array{Float64}(14)
const _ZUNHJ_C = Array{Float64}(105)
const _ZUNHJ_BR = Array{Float64}(14)
const _ZUNHJ_BETA = Array{Float64}(210)
const _ZUNHJ_AR = Array{Float64}(14)
const _ZUNHJ_AP = Array{Float64}(30)
const _ZUNHJ_ALFA = Array{Float64}(180)
function ZUNHJ(ZR::Float64,ZI::Float64,FNU::Float64,IPMTR::Int32,TOL::Float64,PHIR::Float64,PHII::Float64,ARGR::Float64,ARGI::Float64,ZETA1R::Float64,ZETA1I::Float64,ZETA2R::Float64,ZETA2I::Float64,ASUMR::Float64,ASUMI::Float64,BSUMR::Float64,BSUMI::Float64)
    AC::Float64 = zero(Float64)
    const ALFA = _ZUNHJ_ALFA
    ANG::Float64 = zero(Float64)
    const AP = _ZUNHJ_AP
    const AR = _ZUNHJ_AR
    ATOL::Float64 = zero(Float64)
    AW2::Float64 = zero(Float64)
    AZTH::Float64 = zero(Float64)
    const BETA = _ZUNHJ_BETA
    const BR = _ZUNHJ_BR
    BTOL::Float64 = zero(Float64)
    const C = _ZUNHJ_C
    CONEI::Float64 = zero(Float64)
    CONER::Float64 = zero(Float64)
    const CRI = _ZUNHJ_CRI
    const CRR = _ZUNHJ_CRR
    const DRI = _ZUNHJ_DRI
    const DRR = _ZUNHJ_DRR
    EX1::Float64 = zero(Float64)
    EX2::Float64 = zero(Float64)
    FN13::Float64 = zero(Float64)
    FN23::Float64 = zero(Float64)
    const GAMA = _ZUNHJ_GAMA
    GPI::Float64 = zero(Float64)
    HPI::Float64 = zero(Float64)
    IAS::Int32 = zero(Int32)
    IBS::Int32 = zero(Int32)
    IDUM::Int32 = zero(Int32)
    IS::Int32 = zero(Int32)
    J::Int32 = zero(Int32)
    JR::Int32 = zero(Int32)
    JU::Int32 = zero(Int32)
    K::Int32 = zero(Int32)
    KMAX::Int32 = zero(Int32)
    KP1::Int32 = zero(Int32)
    KS::Int32 = zero(Int32)
    L::Int32 = zero(Int32)
    L1::Int32 = zero(Int32)
    L2::Int32 = zero(Int32)
    LR::Int32 = zero(Int32)
    LRP1::Int32 = zero(Int32)
    M::Int32 = zero(Int32)
    const PI = _ZUNHJ_PI
    PP::Float64 = zero(Float64)
    const PR = _ZUNHJ_PR
    PRZTHI::Float64 = zero(Float64)
    PRZTHR::Float64 = zero(Float64)
    PTFNI::Float64 = zero(Float64)
    PTFNR::Float64 = zero(Float64)
    RAW::Float64 = zero(Float64)
    RAW2::Float64 = zero(Float64)
    RAZTH::Float64 = zero(Float64)
    RFN13::Float64 = zero(Float64)
    RFNU::Float64 = zero(Float64)
    RFNU2::Float64 = zero(Float64)
    RTZTI::Float64 = zero(Float64)
    RTZTR::Float64 = zero(Float64)
    RZTHI::Float64 = zero(Float64)
    RZTHR::Float64 = zero(Float64)
    STI::Float64 = zero(Float64)
    STR::Float64 = zero(Float64)
    SUMAI::Float64 = zero(Float64)
    SUMAR::Float64 = zero(Float64)
    SUMBI::Float64 = zero(Float64)
    SUMBR::Float64 = zero(Float64)
    T2I::Float64 = zero(Float64)
    T2R::Float64 = zero(Float64)
    TEST::Float64 = zero(Float64)
    TFNI::Float64 = zero(Float64)
    TFNR::Float64 = zero(Float64)
    THPI::Float64 = zero(Float64)
    TZAI::Float64 = zero(Float64)
    TZAR::Float64 = zero(Float64)
    const UPI = _ZUNHJ_UPI
    const UPR = _ZUNHJ_UPR
    W2I::Float64 = zero(Float64)
    W2R::Float64 = zero(Float64)
    WI::Float64 = zero(Float64)
    WR::Float64 = zero(Float64)
    ZAI::Float64 = zero(Float64)
    ZAR::Float64 = zero(Float64)
    ZBI::Float64 = zero(Float64)
    ZBR::Float64 = zero(Float64)
    ZCI::Float64 = zero(Float64)
    ZCR::Float64 = zero(Float64)
    ZEROI::Float64 = zero(Float64)
    ZEROR::Float64 = zero(Float64)
    ZETAI::Float64 = zero(Float64)
    ZETAR::Float64 = zero(Float64)
    ZTHI::Float64 = zero(Float64)
    ZTHR::Float64 = zero(Float64)
    begin 
        AR[Int32(1)] = 1.0
        AR[Int32(2)] = 0.10416666666666667
        AR[Int32(3)] = 0.08355034722222222
        AR[Int32(4)] = 0.12822657455632716
        AR[Int32(5)] = 0.29184902646414046
        AR[Int32(6)] = 0.8816272674437576
        AR[Int32(7)] = 3.3214082818627677
        AR[Int32(8)] = 14.995762986862555
        AR[Int32(9)] = 78.92301301158652
        AR[Int32(10)] = 474.4515388682643
        AR[Int32(11)] = 3207.490090890662
        AR[Int32(12)] = 24086.549640874004
        AR[Int32(13)] = 198923.1191695098
        AR[Int32(14)] = 1.7919020077753437e6
    end
    begin 
        BR[Int32(1)] = 1.0
        BR[Int32(2)] = -0.14583333333333334
        BR[Int32(3)] = -0.09874131944444445
        BR[Int32(4)] = -0.14331205391589505
        BR[Int32(5)] = -0.31722720267841353
        BR[Int32(6)] = -0.9424291479571203
        BR[Int32(7)] = -3.5112030408263544
        BR[Int32(8)] = -15.727263620368046
        BR[Int32(9)] = -82.28143909718594
        BR[Int32(10)] = -492.3553705236705
        BR[Int32(11)] = -3316.2185685479726
        BR[Int32(12)] = -24827.67424520859
        BR[Int32(13)] = -204526.5873151298
        BR[Int32(14)] = -1.83844491706821e6
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
    end
    begin 
        ALFA[Int32(1)] = -0.0044444444444444444
        ALFA[Int32(2)] = -0.000922077922077922
        ALFA[Int32(3)] = -8.848928848928849e-5
        ALFA[Int32(4)] = 0.00016592768783244973
        ALFA[Int32(5)] = 0.0002466913727417929
        ALFA[Int32(6)] = 0.0002659955893462548
        ALFA[Int32(7)] = 0.00026182429706150096
        ALFA[Int32(8)] = 0.0002487304373446556
        ALFA[Int32(9)] = 0.00023272104008323209
        ALFA[Int32(10)] = 0.00021636248571236508
        ALFA[Int32(11)] = 0.00020073885876275234
        ALFA[Int32(12)] = 0.00018626763663754517
        ALFA[Int32(13)] = 0.0001730607759178765
        ALFA[Int32(14)] = 0.00016109170592901574
        ALFA[Int32(15)] = 0.00015027477416090814
        ALFA[Int32(16)] = 0.0001405034973912698
        ALFA[Int32(17)] = 0.0001316688165459228
        ALFA[Int32(18)] = 0.00012366744559825325
        ALFA[Int32(19)] = 0.00011640527147473791
        ALFA[Int32(20)] = 0.00010979829837271337
        ALFA[Int32(21)] = 0.00010377241042299283
        ALFA[Int32(22)] = 9.826260783693634e-5
    end
    begin 
        ALFA[Int32(23)] = 9.321205172495032e-5
        ALFA[Int32(24)] = 8.857108524787117e-5
        ALFA[Int32(25)] = 8.429631057157003e-5
        ALFA[Int32(26)] = 8.034975484077912e-5
        ALFA[Int32(27)] = 7.669813453592074e-5
        ALFA[Int32(28)] = 7.331221574817778e-5
        ALFA[Int32(29)] = 7.016626251631414e-5
        ALFA[Int32(30)] = 6.723756337901603e-5
        ALFA[Int32(31)] = 0.000693735541354589
        ALFA[Int32(32)] = 0.00023224174518292166
        ALFA[Int32(33)] = -1.419862735566912e-5
        ALFA[Int32(34)] = -0.00011644493167204864
        ALFA[Int32(35)] = -0.00015080355805304876
        ALFA[Int32(36)] = -0.00015512192491809622
        ALFA[Int32(37)] = -0.00014680975664646556
        ALFA[Int32(38)] = -0.00013381550386749137
        ALFA[Int32(39)] = -0.00011974497568425405
        ALFA[Int32(40)] = -0.00010618431920797402
        ALFA[Int32(41)] = -9.376995498911944e-5
        ALFA[Int32(42)] = -8.269230455881933e-5
        ALFA[Int32(43)] = -7.293743481552213e-5
        ALFA[Int32(44)] = -6.440423577210163e-5
    end
    begin 
        ALFA[Int32(45)] = -5.69611566009369e-5
        ALFA[Int32(46)] = -5.0473104430356164e-5
        ALFA[Int32(47)] = -4.481348680088828e-5
        ALFA[Int32(48)] = -3.9868872771759884e-5
        ALFA[Int32(49)] = -3.554005329720425e-5
        ALFA[Int32(50)] = -3.174142566090225e-5
        ALFA[Int32(51)] = -2.839967939041748e-5
        ALFA[Int32(52)] = -2.5452272063487058e-5
        ALFA[Int32(53)] = -2.2845929716472455e-5
        ALFA[Int32(54)] = -2.053527531064806e-5
        ALFA[Int32(55)] = -1.848162176276661e-5
        ALFA[Int32(56)] = -1.665193300213938e-5
        ALFA[Int32(57)] = -1.5017941298011949e-5
        ALFA[Int32(58)] = -1.3555403137904052e-5
        ALFA[Int32(59)] = -1.2243474647385812e-5
        ALFA[Int32(60)] = -1.1064188481130817e-5
        ALFA[Int32(61)] = -0.00035421197145774384
        ALFA[Int32(62)] = -0.00015616126394515941
        ALFA[Int32(63)] = 3.044655035949364e-5
        ALFA[Int32(64)] = 0.0001301986557732427
        ALFA[Int32(65)] = 0.00016747110669971228
        ALFA[Int32(66)] = 0.00017022258768359256
    end
    begin 
        ALFA[Int32(67)] = 0.00015650142760859472
        ALFA[Int32(68)] = 0.00013633917097744512
        ALFA[Int32(69)] = 0.00011488669202982512
        ALFA[Int32(70)] = 9.458690930346882e-5
        ALFA[Int32(71)] = 7.644984192508983e-5
        ALFA[Int32(72)] = 6.0757033496519734e-5
        ALFA[Int32(73)] = 4.743942992905088e-5
        ALFA[Int32(74)] = 3.627575120053443e-5
        ALFA[Int32(75)] = 2.699397149792249e-5
        ALFA[Int32(76)] = 1.9321093824793926e-5
        ALFA[Int32(77)] = 1.3005667479396321e-5
        ALFA[Int32(78)] = 7.826208667444966e-6
        ALFA[Int32(79)] = 3.592574858193516e-6
        ALFA[Int32(80)] = 1.4404004981425182e-7
        ALFA[Int32(81)] = -2.653967696979391e-6
        ALFA[Int32(82)] = -4.913468670984859e-6
        ALFA[Int32(83)] = -6.727392960912483e-6
        ALFA[Int32(84)] = -8.17269379678658e-6
        ALFA[Int32(85)] = -9.313047150935612e-6
        ALFA[Int32(86)] = -1.0201141879801643e-5
        ALFA[Int32(87)] = -1.0880596251059288e-5
        ALFA[Int32(88)] = -1.1387548150960355e-5
    end
    begin 
        ALFA[Int32(89)] = -1.1751967567455642e-5
        ALFA[Int32(90)] = -1.1998736487094414e-5
        ALFA[Int32(91)] = 0.0003781941992017729
        ALFA[Int32(92)] = 0.00020247195276181616
        ALFA[Int32(93)] = -6.379385063188624e-5
        ALFA[Int32(94)] = -0.0002385982306030059
        ALFA[Int32(95)] = -0.0003109162560273616
        ALFA[Int32(96)] = -0.00031368011524757634
        ALFA[Int32(97)] = -0.0002789502737913234
        ALFA[Int32(98)] = -0.00022856408261914138
        ALFA[Int32(99)] = -0.00017524528034084676
        ALFA[Int32(100)] = -0.00012554406306069035
        ALFA[Int32(101)] = -8.229828728202083e-5
        ALFA[Int32(102)] = -4.628607305881165e-5
        ALFA[Int32(103)] = -1.7233430236696227e-5
        ALFA[Int32(104)] = 5.6069048230460226e-6
        ALFA[Int32(105)] = 2.313954431482868e-5
        ALFA[Int32(106)] = 3.626427458567939e-5
        ALFA[Int32(107)] = 4.5800612449018877e-5
        ALFA[Int32(108)] = 5.2459529495911405e-5
        ALFA[Int32(109)] = 5.683962085458153e-5
        ALFA[Int32(110)] = 5.9434982039310406e-5
    end
    begin 
        ALFA[Int32(111)] = 6.0647852757842175e-5
        ALFA[Int32(112)] = 6.080239077884365e-5
        ALFA[Int32(113)] = 6.0157789453946036e-5
        ALFA[Int32(114)] = 5.891996573446985e-5
        ALFA[Int32(115)] = 5.72515823777593e-5
        ALFA[Int32(116)] = 5.528043755858526e-5
        ALFA[Int32(117)] = 5.310637738028802e-5
        ALFA[Int32(118)] = 5.080693020123257e-5
        ALFA[Int32(119)] = 4.8441864762009484e-5
        ALFA[Int32(120)] = 4.6056858160747536e-5
        ALFA[Int32(121)] = -0.0006911413972882942
        ALFA[Int32(122)] = -0.0004299766330588719
        ALFA[Int32(123)] = 0.000183067735980039
        ALFA[Int32(124)] = 0.0006600881475420142
        ALFA[Int32(125)] = 0.0008759649699511859
        ALFA[Int32(126)] = 0.0008773352359582355
        ALFA[Int32(127)] = 0.0007493695853789907
        ALFA[Int32(128)] = 0.000563832329756981
        ALFA[Int32(129)] = 0.0003680593199714432
        ALFA[Int32(130)] = 0.0001884645355144556
    end
    begin 
        ALFA[Int32(131)] = 3.7066305766490415e-5
        ALFA[Int32(132)] = -8.28520220232137e-5
        ALFA[Int32(133)] = -0.000172751952869173
        ALFA[Int32(134)] = -0.00023631487360587297
        ALFA[Int32(135)] = -0.0002779661506949067
        ALFA[Int32(136)] = -0.00030207951415545694
        ALFA[Int32(137)] = -0.0003125947126438201
        ALFA[Int32(138)] = -0.00031287255875806717
        ALFA[Int32(139)] = -0.0003056780384663244
        ALFA[Int32(140)] = -0.0002932264706145573
        ALFA[Int32(141)] = -0.0002772556555829348
        ALFA[Int32(142)] = -0.0002591039284670317
        ALFA[Int32(143)] = -0.00023978401439648034
        ALFA[Int32(144)] = -0.00022004826004542284
        ALFA[Int32(145)] = -0.00020044391109497149
        ALFA[Int32(146)] = -0.00018135869221097068
        ALFA[Int32(147)] = -0.00016305767447865748
        ALFA[Int32(148)] = -0.00014571267217520584
        ALFA[Int32(149)] = -0.0001294254219839246
        ALFA[Int32(150)] = -0.00011424569194244596
    end
    begin 
        ALFA[Int32(151)] = 0.0019282196424877589
        ALFA[Int32(152)] = 0.0013559257630202223
        ALFA[Int32(153)] = -0.000717858090421303
        ALFA[Int32(154)] = -0.0025808480257527035
        ALFA[Int32(155)] = -0.0034927113082616847
        ALFA[Int32(156)] = -0.003469862993409606
        ALFA[Int32(157)] = -0.002822852333513102
        ALFA[Int32(158)] = -0.0018810307640489134
        ALFA[Int32(159)] = -0.0008895317183839476
        ALFA[Int32(160)] = 3.8791210263103525e-6
        ALFA[Int32(161)] = 0.0007286885401196914
        ALFA[Int32(162)] = 0.0012656637305345775
        ALFA[Int32(163)] = 0.0016251815837267443
        ALFA[Int32(164)] = 0.0018320315321637317
        ALFA[Int32(165)] = 0.0019158838899052792
        ALFA[Int32(166)] = 0.0019058884675554615
        ALFA[Int32(167)] = 0.0018279898242182574
        ALFA[Int32(168)] = 0.0017038950642112153
        ALFA[Int32(169)] = 0.0015509712717109768
        ALFA[Int32(170)] = 0.0013826142185227616
    end
    begin 
        ALFA[Int32(171)] = 0.0012088142423006478
        ALFA[Int32(172)] = 0.0010367653263834496
        ALFA[Int32(173)] = 0.0008714379180686191
        ALFA[Int32(174)] = 0.000716080155297701
        ALFA[Int32(175)] = 0.0005726370025581294
        ALFA[Int32(176)] = 0.0004420898194658023
        ALFA[Int32(177)] = 0.00032472494850309055
        ALFA[Int32(178)] = 0.00022034204273024659
        ALFA[Int32(179)] = 0.00012841289840135388
        ALFA[Int32(180)] = 4.8200592455209545e-5
    end
    begin 
        BETA[Int32(1)] = 0.01799887214135533
        BETA[Int32(2)] = 0.005599649110643881
        BETA[Int32(3)] = 0.0028850140223113277
        BETA[Int32(4)] = 0.0018009660676105393
        BETA[Int32(5)] = 0.001247531105891992
        BETA[Int32(6)] = 0.0009228788765729383
        BETA[Int32(7)] = 0.0007144304217272874
        BETA[Int32(8)] = 0.0005717872817897049
        BETA[Int32(9)] = 0.00046943100760648155
        BETA[Int32(10)] = 0.00039323283546291665
        BETA[Int32(11)] = 0.0003348188893182977
        BETA[Int32(12)] = 0.00028895214849575154
        BETA[Int32(13)] = 0.0002522116155495733
        BETA[Int32(14)] = 0.00022228058079888332
        BETA[Int32(15)] = 0.0001975418380330625
        BETA[Int32(16)] = 0.00017683685501971802
        BETA[Int32(17)] = 0.0001593168996618211
        BETA[Int32(18)] = 0.00014434793019733397
        BETA[Int32(19)] = 0.0001314480681199654
        BETA[Int32(20)] = 0.00012024544494930288
        BETA[Int32(21)] = 0.0001104491445045994
        BETA[Int32(22)] = 0.00010182877074056726
    end
    begin 
        BETA[Int32(23)] = 9.419982242042375e-5
        BETA[Int32(24)] = 8.741305457538345e-5
        BETA[Int32(25)] = 8.134662621628014e-5
        BETA[Int32(26)] = 7.590022696462193e-5
        BETA[Int32(27)] = 7.099063006341535e-5
        BETA[Int32(28)] = 6.654828748424682e-5
        BETA[Int32(29)] = 6.25146958969275e-5
        BETA[Int32(30)] = 5.884033944262518e-5
        BETA[Int32(31)] = -0.0014928295321342917
        BETA[Int32(32)] = -0.0008782047095463894
        BETA[Int32(33)] = -0.0005029165495720346
        BETA[Int32(34)] = -0.000294822138512746
        BETA[Int32(35)] = -0.00017546399697078284
        BETA[Int32(36)] = -0.00010400855046081644
        BETA[Int32(37)] = -5.961419530464579e-5
        BETA[Int32(38)] = -3.1203892907609836e-5
        BETA[Int32(39)] = -1.2608973598023005e-5
        BETA[Int32(40)] = -2.4289260857573037e-7
        BETA[Int32(41)] = 8.059961654142736e-6
        BETA[Int32(42)] = 1.3650700926214739e-5
        BETA[Int32(43)] = 1.7396412547292627e-5
        BETA[Int32(44)] = 1.9867297884213378e-5
    end
    begin 
        BETA[Int32(45)] = 2.1446326379082263e-5
        BETA[Int32(46)] = 2.2395465923245652e-5
        BETA[Int32(47)] = 2.2896778381471263e-5
        BETA[Int32(48)] = 2.307853898111778e-5
        BETA[Int32(49)] = 2.3032197608090914e-5
        BETA[Int32(50)] = 2.2823607372034874e-5
        BETA[Int32(51)] = 2.250058811052924e-5
        BETA[Int32(52)] = 2.2098101536199144e-5
        BETA[Int32(53)] = 2.164184274481039e-5
        BETA[Int32(54)] = 2.1150764925622083e-5
        BETA[Int32(55)] = 2.0638874978217072e-5
        BETA[Int32(56)] = 2.0116524199708165e-5
        BETA[Int32(57)] = 1.9591345014117925e-5
        BETA[Int32(58)] = 1.9068936791043675e-5
        BETA[Int32(59)] = 1.8553371964163667e-5
        BETA[Int32(60)] = 1.804757222596742e-5
        BETA[Int32(61)] = 0.0005522130767212928
        BETA[Int32(62)] = 0.00044793258155238465
        BETA[Int32(63)] = 0.0002795206539920206
        BETA[Int32(64)] = 0.0001524681561984466
        BETA[Int32(65)] = 6.932711056570436e-5
        BETA[Int32(66)] = 1.762586830699914e-5
    end
    begin 
        BETA[Int32(67)] = -1.3574499634326914e-5
        BETA[Int32(68)] = -3.179724133504272e-5
        BETA[Int32(69)] = -4.188618616966934e-5
        BETA[Int32(70)] = -4.6900488937914104e-5
        BETA[Int32(71)] = -4.8766544741378735e-5
        BETA[Int32(72)] = -4.8701003118673505e-5
        BETA[Int32(73)] = -4.747556208900866e-5
        BETA[Int32(74)] = -4.558130581386284e-5
        BETA[Int32(75)] = -4.33309644511266e-5
        BETA[Int32(76)] = -4.0923019315775034e-5
        BETA[Int32(77)] = -3.848226386032213e-5
        BETA[Int32(78)] = -3.608571675354105e-5
        BETA[Int32(79)] = -3.377933061233674e-5
        BETA[Int32(80)] = -3.158885607721096e-5
        BETA[Int32(81)] = -2.952695617508073e-5
        BETA[Int32(82)] = -2.7597891482833575e-5
        BETA[Int32(83)] = -2.5800617466688372e-5
        BETA[Int32(84)] = -2.413083567612802e-5
        BETA[Int32(85)] = -2.2582350951834605e-5
        BETA[Int32(86)] = -2.1147965676891298e-5
        BETA[Int32(87)] = -1.9820063888529493e-5
        BETA[Int32(88)] = -1.8590987080106508e-5
    end
    begin 
        BETA[Int32(89)] = -1.7453269984421023e-5
        BETA[Int32(90)] = -1.63997823854498e-5
        BETA[Int32(91)] = -0.0004746177965599598
        BETA[Int32(92)] = -0.0004778645671473215
        BETA[Int32(93)] = -0.00032039022806703763
        BETA[Int32(94)] = -0.00016110501611996228
        BETA[Int32(95)] = -4.257781012854352e-5
        BETA[Int32(96)] = 3.445712942949675e-5
        BETA[Int32(97)] = 7.97092684075675e-5
        BETA[Int32(98)] = 0.0001031382367082722
        BETA[Int32(99)] = 0.00011246677526220416
        BETA[Int32(100)] = 0.0001131036421084814
        BETA[Int32(101)] = 0.00010865163484877427
        BETA[Int32(102)] = 0.00010143795159766197
        BETA[Int32(103)] = 9.29298396593364e-5
        BETA[Int32(104)] = 8.4029313301609e-5
        BETA[Int32(105)] = 7.52727991349134e-5
        BETA[Int32(106)] = 6.696325219757309e-5
        BETA[Int32(107)] = 5.925645473231947e-5
        BETA[Int32(108)] = 5.2216930882697554e-5
        BETA[Int32(109)] = 4.585394851653606e-5
        BETA[Int32(110)] = 4.014455138914868e-5
    end
    begin 
        BETA[Int32(111)] = 3.504817300313281e-5
        BETA[Int32(112)] = 3.0515799503434667e-5
        BETA[Int32(113)] = 2.6495611995051603e-5
        BETA[Int32(114)] = 2.2936363369099816e-5
        BETA[Int32(115)] = 1.9789305666402162e-5
        BETA[Int32(116)] = 1.7009198463641262e-5
        BETA[Int32(117)] = 1.45547428261524e-5
        BETA[Int32(118)] = 1.238866409958784e-5
        BETA[Int32(119)] = 1.0477587607658323e-5
        BETA[Int32(120)] = 8.791799549784793e-6
        BETA[Int32(121)] = 0.0007364658105725784
        BETA[Int32(122)] = 0.000872790805146194
        BETA[Int32(123)] = 0.0006226148625731351
        BETA[Int32(124)] = 0.00028599815419430417
        BETA[Int32(125)] = 3.847376728793661e-6
        BETA[Int32(126)] = -0.00018790600363697156
        BETA[Int32(127)] = -0.00029760364659455455
        BETA[Int32(128)] = -0.00034599812683265633
        BETA[Int32(129)] = -0.00035338247091603773
        BETA[Int32(130)] = -0.00033571563577504876
    end
    begin 
        BETA[Int32(131)] = -0.0003043211247890398
        BETA[Int32(132)] = -0.00026672272304761283
        BETA[Int32(133)] = -0.00022765421412281953
        BETA[Int32(134)] = -0.00018992261185456235
        BETA[Int32(135)] = -0.00015505891859909386
        BETA[Int32(136)] = -0.00012377824076187363
        BETA[Int32(137)] = -9.629261477176441e-5
        BETA[Int32(138)] = -7.251783277144253e-5
        BETA[Int32(139)] = -5.220700288956338e-5
        BETA[Int32(140)] = -3.5034775051190054e-5
        BETA[Int32(141)] = -2.0648976103555174e-5
        BETA[Int32(142)] = -8.701060968497671e-6
        BETA[Int32(143)] = 1.136986866751003e-6
        BETA[Int32(144)] = 9.164264741227788e-6
        BETA[Int32(145)] = 1.564777854288726e-5
        BETA[Int32(146)] = 2.0822362948246685e-5
        BETA[Int32(147)] = 2.4892338100459516e-5
        BETA[Int32(148)] = 2.803405095741463e-5
        BETA[Int32(149)] = 3.039877746298619e-5
        BETA[Int32(150)] = 3.211567314067006e-5
    end
    begin 
        BETA[Int32(151)] = -0.0018018219196388571
        BETA[Int32(152)] = -0.0024340296293804253
        BETA[Int32(153)] = -0.001834226635498568
        BETA[Int32(154)] = -0.0007622045963540097
        BETA[Int32(155)] = 0.00023907947525692722
        BETA[Int32(156)] = 0.0009492661171768811
        BETA[Int32(157)] = 0.0013446744970154036
        BETA[Int32(158)] = 0.0014845749525944918
        BETA[Int32(159)] = 0.001447323398306176
        BETA[Int32(160)] = 0.0013026826128565718
        BETA[Int32(161)] = 0.0011035159737564268
        BETA[Int32(162)] = 0.0008860474404197917
        BETA[Int32(163)] = 0.0006730732081656654
        BETA[Int32(164)] = 0.00047760387285658237
        BETA[Int32(165)] = 0.00030599192635878935
        BETA[Int32(166)] = 0.00016031569459472162
        BETA[Int32(167)] = 4.007495552706133e-5
        BETA[Int32(168)] = -5.666074616352516e-5
        BETA[Int32(169)] = -0.00013250618677298264
        BETA[Int32(170)] = -0.00019029618798961406
    end
    begin 
        BETA[Int32(171)] = -0.0002328114503769374
        BETA[Int32(172)] = -0.00026262881146466884
        BETA[Int32(173)] = -0.00028205046986759866
        BETA[Int32(174)] = -0.00029308156319286116
        BETA[Int32(175)] = -0.0002974359621763166
        BETA[Int32(176)] = -0.0002965573342393481
        BETA[Int32(177)] = -0.0002916473633120909
        BETA[Int32(178)] = -0.0002836962038377342
        BETA[Int32(179)] = -0.00027351231709567335
        BETA[Int32(180)] = -0.0002617501558067686
        BETA[Int32(181)] = 0.006385858912120509
        BETA[Int32(182)] = 0.00962374215806378
        BETA[Int32(183)] = 0.0076187806120700105
        BETA[Int32(184)] = 0.0028321905554562804
        BETA[Int32(185)] = -0.002098413520127201
        BETA[Int32(186)] = -0.005738267642166265
        BETA[Int32(187)] = -0.0077080424449541465
        BETA[Int32(188)] = -0.008210116922648444
        BETA[Int32(189)] = -0.007658245203469054
        BETA[Int32(190)] = -0.006472097293910452
    end
    begin 
        BETA[Int32(191)] = -0.004991324120049665
        BETA[Int32(192)] = -0.0034561228971313326
        BETA[Int32(193)] = -0.002017855800141708
        BETA[Int32(194)] = -0.0007594306867819614
        BETA[Int32(195)] = 0.0002841736315238591
        BETA[Int32(196)] = 0.001108916675863374
        BETA[Int32(197)] = 0.0017290149387272878
        BETA[Int32(198)] = 0.002168125908026847
        BETA[Int32(199)] = 0.002453577104945397
        BETA[Int32(200)] = 0.0026128182105833488
        BETA[Int32(201)] = 0.002671410396562769
        BETA[Int32(202)] = 0.0026520307339598045
        BETA[Int32(203)] = 0.002574116528772873
        BETA[Int32(204)] = 0.0024538912623609443
        BETA[Int32(205)] = 0.002304600580717955
        BETA[Int32(206)] = 0.0021368483768671267
        BETA[Int32(207)] = 0.001958965284788709
        BETA[Int32(208)] = 0.0017773700867945441
        BETA[Int32(209)] = 0.0015969028076583906
        BETA[Int32(210)] = 0.0014211197566443854
    end
    begin 
        GAMA[Int32(1)] = 0.6299605249474366
        GAMA[Int32(2)] = 0.25198420997897464
        GAMA[Int32(3)] = 0.15479030041565583
        GAMA[Int32(4)] = 0.11071306241615901
        GAMA[Int32(5)] = 0.08573093955273949
        GAMA[Int32(6)] = 0.06971613169586843
        GAMA[Int32(7)] = 0.05860856718937136
        GAMA[Int32(8)] = 0.05046988735363107
        GAMA[Int32(9)] = 0.04426005806891548
        GAMA[Int32(10)] = 0.039372066154350994
        GAMA[Int32(11)] = 0.03542831959244554
        GAMA[Int32(12)] = 0.032181885750209825
        GAMA[Int32(13)] = 0.029464624079115768
        GAMA[Int32(14)] = 0.027158167711293448
        GAMA[Int32(15)] = 0.025176827297386177
        GAMA[Int32(16)] = 0.02345707553060789
        GAMA[Int32(17)] = 0.02195083901349072
        GAMA[Int32(18)] = 0.020621082823564625
        GAMA[Int32(19)] = 0.019438824089788084
        GAMA[Int32(20)] = 0.018381063380068317
        GAMA[Int32(21)] = 0.017429321323196318
        GAMA[Int32(22)] = 0.016568583778661234
    end
    begin 
        GAMA[Int32(23)] = 0.015786528598791844
        GAMA[Int32(24)] = 0.01507295014940956
        GAMA[Int32(25)] = 0.014419325083995464
        GAMA[Int32(26)] = 0.013818480573534178
        GAMA[Int32(27)] = 0.013264337899427657
        GAMA[Int32(28)] = 0.012751712197049864
        GAMA[Int32(29)] = 0.012276154531876277
        GAMA[Int32(30)] = 0.01183382623984824
    end
    begin 
        EX1 = 0.3333333333333333
        EX2 = 0.6666666666666666
        HPI = 1.5707963267948966
        GPI = 3.141592653589793
        THPI = 4.71238898038469
    end
    begin 
        ZEROR = 0.0
        ZEROI = 0.0
        CONER = 1.0
        CONEI = 0.0
    end
    RFNU = 1.0 / FNU
    TEST = D1MACH1 * 1000.0
    AC = FNU * TEST
    if DABS(ZR) > AC || DABS(ZI) > AC
        @goto line15
    end
    ZETA1R = 2.0 * DABS(DLOG(TEST)) + FNU
    ZETA1I = 0.0
    ZETA2R = FNU
    ZETA2I = 0.0
    PHIR = 1.0
    PHII = 0.0
    ARGR = 1.0
    ARGI = 0.0
    return (PHIR,PHII,ARGR,ARGI,ZETA1R,ZETA1I,ZETA2R,ZETA2I,ASUMR,ASUMI,BSUMR,BSUMI)
    @label line15
    ZBR = ZR * RFNU
    ZBI = ZI * RFNU
    RFNU2 = RFNU * RFNU
    FN13 = FNU^EX1
    FN23 = FN13 * FN13
    RFN13 = 1.0 / FN13
    W2R = (CONER - ZBR * ZBR) + ZBI * ZBI
    W2I = (CONEI - ZBR * ZBI) - ZBR * ZBI
    AW2 = ZABS(COMPLEX(W2R,W2I))
    if AW2 > 0.25
        @goto line130
    end
    K = Int32(1)
    PR[Int32(1)] = CONER
    PI[Int32(1)] = CONEI
    SUMAR = GAMA[Int32(1)]
    SUMAI = ZEROI
    AP[Int32(1)] = 1.0
    if AW2 < TOL
        @goto line20
    end
    for K = Int32(2):Int32(30)
        PR[K] = PR[K - Int32(1)] * W2R - PI[K - Int32(1)] * W2I
        PI[K] = PR[K - Int32(1)] * W2I + PI[K - Int32(1)] * W2R
        SUMAR = SUMAR + PR[K] * GAMA[K]
        SUMAI = SUMAI + PI[K] * GAMA[K]
        AP[K] = AP[K - Int32(1)] * AW2
        if AP[K] < TOL
            @goto line20
        end
        @label line10
    end
    K = Int32(30)
    @label line20
    KMAX = K
    ZETAR = W2R * SUMAR - W2I * SUMAI
    ZETAI = W2R * SUMAI + W2I * SUMAR
    ARGR = ZETAR * FN23
    ARGI = ZETAI * FN23
    (ZAR,ZAI) = ZSQRT(SUMAR,SUMAI,ZAR,ZAI)
    (STR,STI) = ZSQRT(W2R,W2I,STR,STI)
    ZETA2R = STR * FNU
    ZETA2I = STI * FNU
    STR = CONER + EX2 * (ZETAR * ZAR - ZETAI * ZAI)
    STI = CONEI + EX2 * (ZETAR * ZAI + ZETAI * ZAR)
    ZETA1R = STR * ZETA2R - STI * ZETA2I
    ZETA1I = STR * ZETA2I + STI * ZETA2R
    ZAR = ZAR + ZAR
    ZAI = ZAI + ZAI
    (STR,STI) = ZSQRT(ZAR,ZAI,STR,STI)
    PHIR = STR * RFN13
    PHII = STI * RFN13
    if IPMTR == Int32(1)
        @goto line120
    end
    SUMBR = ZEROR
    SUMBI = ZEROI
    for K = Int32(1):KMAX
        SUMBR = SUMBR + PR[K] * BETA[K]
        SUMBI = SUMBI + PI[K] * BETA[K]
        @label line30
    end
    ASUMR = ZEROR
    ASUMI = ZEROI
    BSUMR = SUMBR
    BSUMI = SUMBI
    L1 = Int32(0)
    L2 = Int32(30)
    BTOL = TOL * (DABS(BSUMR) + DABS(BSUMI))
    ATOL = TOL
    PP = 1.0
    IAS = Int32(0)
    IBS = Int32(0)
    if RFNU2 < TOL
        @goto line110
    end
    for IS = Int32(2):Int32(7)
        ATOL = ATOL / RFNU2
        PP = PP * RFNU2
        if IAS == Int32(1)
            @goto line60
        end
        SUMAR = ZEROR
        SUMAI = ZEROI
        for K = Int32(1):KMAX
            M = L1 + K
            SUMAR = SUMAR + PR[K] * ALFA[M]
            SUMAI = SUMAI + PI[K] * ALFA[M]
            if AP[K] < ATOL
                @goto line50
            end
            @label line40
        end
        @label line50
        ASUMR = ASUMR + SUMAR * PP
        ASUMI = ASUMI + SUMAI * PP
        if PP < TOL
            IAS = Int32(1)
        end
        @label line60
        if IBS == Int32(1)
            @goto line90
        end
        SUMBR = ZEROR
        SUMBI = ZEROI
        for K = Int32(1):KMAX
            M = L2 + K
            SUMBR = SUMBR + PR[K] * BETA[M]
            SUMBI = SUMBI + PI[K] * BETA[M]
            if AP[K] < ATOL
                @goto line80
            end
            @label line70
        end
        @label line80
        BSUMR = BSUMR + SUMBR * PP
        BSUMI = BSUMI + SUMBI * PP
        if PP < BTOL
            IBS = Int32(1)
        end
        @label line90
        if IAS == Int32(1) && IBS == Int32(1)
            @goto line110
        end
        L1 = L1 + Int32(30)
        L2 = L2 + Int32(30)
        @label line100
    end
    @label line110
    ASUMR = ASUMR + CONER
    PP = RFNU * RFN13
    BSUMR = BSUMR * PP
    BSUMI = BSUMI * PP
    @label line120
    return (PHIR,PHII,ARGR,ARGI,ZETA1R,ZETA1I,ZETA2R,ZETA2I,ASUMR,ASUMI,BSUMR,BSUMI)
    @label line130
    (WR,WI) = ZSQRT(W2R,W2I,WR,WI)
    if WR < 0.0
        WR = 0.0
    end
    if WI < 0.0
        WI = 0.0
    end
    STR = CONER + WR
    STI = WI
    (ZAR,ZAI) = ZDIV(STR,STI,ZBR,ZBI,ZAR,ZAI)
    (ZCR,ZCI,IDUM) = ZLOG(ZAR,ZAI,ZCR,ZCI,IDUM)
    if ZCI < 0.0
        ZCI = 0.0
    end
    if ZCI > HPI
        ZCI = HPI
    end
    if ZCR < 0.0
        ZCR = 0.0
    end
    ZTHR = (ZCR - WR) * 1.5
    ZTHI = (ZCI - WI) * 1.5
    ZETA1R = ZCR * FNU
    ZETA1I = ZCI * FNU
    ZETA2R = WR * FNU
    ZETA2I = WI * FNU
    AZTH = ZABS(COMPLEX(ZTHR,ZTHI))
    ANG = THPI
    if ZTHR >= 0.0 && ZTHI < 0.0
        @goto line140
    end
    ANG = HPI
    if ZTHR == 0.0
        @goto line140
    end
    ANG = DATAN(ZTHI / ZTHR)
    if ZTHR < 0.0
        ANG = ANG + GPI
    end
    @label line140
    PP = AZTH^EX2
    ANG = ANG * EX2
    ZETAR = PP * DCOS(ANG)
    ZETAI = PP * DSIN(ANG)
    if ZETAI < 0.0
        ZETAI = 0.0
    end
    ARGR = ZETAR * FN23
    ARGI = ZETAI * FN23
    (RTZTR,RTZTI) = ZDIV(ZTHR,ZTHI,ZETAR,ZETAI,RTZTR,RTZTI)
    (ZAR,ZAI) = ZDIV(RTZTR,RTZTI,WR,WI,ZAR,ZAI)
    TZAR = ZAR + ZAR
    TZAI = ZAI + ZAI
    (STR,STI) = ZSQRT(TZAR,TZAI,STR,STI)
    PHIR = STR * RFN13
    PHII = STI * RFN13
    if IPMTR == Int32(1)
        @goto line120
    end
    RAW = 1.0 / DSQRT(AW2)
    STR = WR * RAW
    STI = -WI * RAW
    TFNR = STR * RFNU * RAW
    TFNI = STI * RFNU * RAW
    RAZTH = 1.0 / AZTH
    STR = ZTHR * RAZTH
    STI = -ZTHI * RAZTH
    RZTHR = STR * RAZTH * RFNU
    RZTHI = STI * RAZTH * RFNU
    ZCR = RZTHR * AR[Int32(2)]
    ZCI = RZTHI * AR[Int32(2)]
    RAW2 = 1.0 / AW2
    STR = W2R * RAW2
    STI = -W2I * RAW2
    T2R = STR * RAW2
    T2I = STI * RAW2
    STR = T2R * C[Int32(2)] + C[Int32(3)]
    STI = T2I * C[Int32(2)]
    UPR[Int32(2)] = STR * TFNR - STI * TFNI
    UPI[Int32(2)] = STR * TFNI + STI * TFNR
    BSUMR = UPR[Int32(2)] + ZCR
    BSUMI = UPI[Int32(2)] + ZCI
    ASUMR = ZEROR
    ASUMI = ZEROI
    if RFNU < TOL
        @goto line220
    end
    PRZTHR = RZTHR
    PRZTHI = RZTHI
    PTFNR = TFNR
    PTFNI = TFNI
    UPR[Int32(1)] = CONER
    UPI[Int32(1)] = CONEI
    PP = 1.0
    BTOL = TOL * (DABS(BSUMR) + DABS(BSUMI))
    KS = Int32(0)
    KP1 = Int32(2)
    L = Int32(3)
    IAS = Int32(0)
    IBS = Int32(0)
    for LR = Int32(2):Int32(2):Int32(12)
        LRP1 = LR + Int32(1)
        for K = LR:LRP1
            KS = KS + Int32(1)
            KP1 = KP1 + Int32(1)
            L = L + Int32(1)
            ZAR = C[L]
            ZAI = ZEROI
            for J = Int32(2):KP1
                L = L + Int32(1)
                STR = (ZAR * T2R - T2I * ZAI) + C[L]
                ZAI = ZAR * T2I + ZAI * T2R
                ZAR = STR
                @label line150
            end
            STR = PTFNR * TFNR - PTFNI * TFNI
            PTFNI = PTFNR * TFNI + PTFNI * TFNR
            PTFNR = STR
            UPR[KP1] = PTFNR * ZAR - PTFNI * ZAI
            UPI[KP1] = PTFNI * ZAR + PTFNR * ZAI
            CRR[KS] = PRZTHR * BR[KS + Int32(1)]
            CRI[KS] = PRZTHI * BR[KS + Int32(1)]
            STR = PRZTHR * RZTHR - PRZTHI * RZTHI
            PRZTHI = PRZTHR * RZTHI + PRZTHI * RZTHR
            PRZTHR = STR
            DRR[KS] = PRZTHR * AR[KS + Int32(2)]
            DRI[KS] = PRZTHI * AR[KS + Int32(2)]
            @label line160
        end
        PP = PP * RFNU2
        if IAS == Int32(1)
            @goto line180
        end
        SUMAR = UPR[LRP1]
        SUMAI = UPI[LRP1]
        JU = LRP1
        for JR = Int32(1):LR
            JU = JU - Int32(1)
            SUMAR = (SUMAR + CRR[JR] * UPR[JU]) - CRI[JR] * UPI[JU]
            SUMAI = SUMAI + CRR[JR] * UPI[JU] + CRI[JR] * UPR[JU]
            @label line170
        end
        ASUMR = ASUMR + SUMAR
        ASUMI = ASUMI + SUMAI
        TEST = DABS(SUMAR) + DABS(SUMAI)
        if PP < TOL && TEST < TOL
            IAS = Int32(1)
        end
        @label line180
        if IBS == Int32(1)
            @goto line200
        end
        SUMBR = (UPR[LR + Int32(2)] + UPR[LRP1] * ZCR) - UPI[LRP1] * ZCI
        SUMBI = UPI[LR + Int32(2)] + UPR[LRP1] * ZCI + UPI[LRP1] * ZCR
        JU = LRP1
        for JR = Int32(1):LR
            JU = JU - Int32(1)
            SUMBR = (SUMBR + DRR[JR] * UPR[JU]) - DRI[JR] * UPI[JU]
            SUMBI = SUMBI + DRR[JR] * UPI[JU] + DRI[JR] * UPR[JU]
            @label line190
        end
        BSUMR = BSUMR + SUMBR
        BSUMI = BSUMI + SUMBI
        TEST = DABS(SUMBR) + DABS(SUMBI)
        if PP < BTOL && TEST < BTOL
            IBS = Int32(1)
        end
        @label line200
        if IAS == Int32(1) && IBS == Int32(1)
            @goto line220
        end
        @label line210
    end
    @label line220
    ASUMR = ASUMR + CONER
    STR = -BSUMR * RFN13
    STI = -BSUMI * RFN13
    (BSUMR,BSUMI) = ZDIV(STR,STI,RTZTR,RTZTI,BSUMR,BSUMI)
    @goto line120
end
