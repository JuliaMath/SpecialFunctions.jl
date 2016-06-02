# airy
@test_approx_eq SpecialFunctions.airy(1.8) 0.0470362168668458052247
@test_approx_eq SpecialFunctions.airyprime(1.8) -0.0685247801186109345638
@test_approx_eq SpecialFunctions.airybi(1.8) 2.595869356743906290060
@test_approx_eq SpecialFunctions.airybiprime(1.8) 2.98554005084659907283
@test_throws SpecialFunctions.AmosException SpecialFunctions.airy(200im)
@test_throws SpecialFunctions.AmosException SpecialFunctions.airybi(200)
z = 1.8 + 1.0im
@test_approx_eq SpecialFunctions.airyx(0, z) SpecialFunctions.airy(0, z) * exp(2/3 * z * sqrt(z))
@test_approx_eq SpecialFunctions.airyx(1, z) SpecialFunctions.airy(1, z) * exp(2/3 * z * sqrt(z))
@test_approx_eq SpecialFunctions.airyx(2, z) SpecialFunctions.airy(2, z) * exp(-abs(real(2/3 * z * sqrt(z))))
@test_approx_eq SpecialFunctions.airyx(3, z) SpecialFunctions.airy(3, z) * exp(-abs(real(2/3 * z * sqrt(z))))

# besselh
true_h133 = 0.30906272225525164362 - 0.53854161610503161800im
@test_approx_eq SpecialFunctions.besselh(3,1,3) true_h133
@test_approx_eq SpecialFunctions.besselh(-3,1,3) -true_h133
@test_approx_eq SpecialFunctions.besselh(3,2,3) conj(true_h133)
@test_approx_eq SpecialFunctions.besselh(-3,2,3) -conj(true_h133)
@test_throws SpecialFunctions.AmosException SpecialFunctions.besselh(1,0)

# besseli
true_i33 = 0.95975362949600785698
@test_approx_eq SpecialFunctions.besseli(3,3) true_i33
@test_approx_eq SpecialFunctions.besseli(-3,3) true_i33
@test_approx_eq SpecialFunctions.besseli(3,-3) -true_i33
@test_approx_eq SpecialFunctions.besseli(-3,-3) -true_i33
@test_throws SpecialFunctions.AmosException SpecialFunctions.besseli(1,1000)

# besselj
@test SpecialFunctions.besselj(0,0) == 1
for i = 1:5
    @test SpecialFunctions.besselj(i,0) == 0
    @test SpecialFunctions.besselj(-i,0) == 0
end

j33 = SpecialFunctions.besselj(3,3.)
@test SpecialFunctions.besselj(3,3) == j33
@test SpecialFunctions.besselj(-3,-3) == j33
@test SpecialFunctions.besselj(-3,3) == -j33
@test SpecialFunctions.besselj(3,-3) == -j33

j43 = SpecialFunctions.besselj(4,3.)
@test SpecialFunctions.besselj(4,3) == j43
@test SpecialFunctions.besselj(-4,-3) == j43
@test SpecialFunctions.besselj(-4,3) == j43
@test SpecialFunctions.besselj(4,-3) == j43

@test_approx_eq j33 0.30906272225525164362
@test_approx_eq j43 0.13203418392461221033
@test_throws DomainError    SpecialFunctions.besselj(0.1, -0.4)
@test_approx_eq SpecialFunctions.besselj(0.1, complex(-0.4)) 0.820421842809028916 + 0.266571215948350899im
@test_approx_eq SpecialFunctions.besselj(3.2, 1.3+0.6im) 0.01135309305831220201 + 0.03927719044393515275im
@test_approx_eq SpecialFunctions.besselj(1, 3im) 3.953370217402609396im
@test_throws SpecialFunctions.AmosException SpecialFunctions.besselj(20,1000im)

# besselk
true_k33 = 0.12217037575718356792
@test_approx_eq SpecialFunctions.besselk(3,3) true_k33
@test_approx_eq SpecialFunctions.besselk(-3,3) true_k33
true_k3m3 = -0.1221703757571835679 - 3.0151549516807985776im
@test_throws DomainError SpecialFunctions.besselk(3,-3)
@test_approx_eq SpecialFunctions.besselk(3,complex(-3)) true_k3m3
@test_approx_eq SpecialFunctions.besselk(-3,complex(-3)) true_k3m3
@test_throws SpecialFunctions.AmosException SpecialFunctions.besselk(200,0.01)
# issue #6564
@test SpecialFunctions.besselk(1.0,0.0) == Inf

# bessely
y33 = SpecialFunctions.bessely(3,3.)
@test SpecialFunctions.bessely(3,3) == y33
@test_approx_eq SpecialFunctions.bessely(-3,3) -y33
@test_approx_eq y33 -0.53854161610503161800
@test_throws DomainError SpecialFunctions.bessely(3,-3)
@test_approx_eq SpecialFunctions.bessely(3,complex(-3)) 0.53854161610503161800 - 0.61812544451050328724im
@test_throws SpecialFunctions.AmosException SpecialFunctions.bessely(200.5,0.1)

# issue #6653
for f in (SpecialFunctions.besselj,SpecialFunctions.bessely,SpecialFunctions.besseli,SpecialFunctions.besselk,SpecialFunctions.hankelh1,SpecialFunctions.hankelh2)
    @test_approx_eq f(0,1) f(0,Complex128(1))
    @test_approx_eq f(0,1) f(0,Complex64(1))
end

# scaled bessel[ijky] and hankelh[12]
for x in (1.0, 0.0, -1.0), y in (1.0, 0.0, -1.0), nu in (1.0, 0.0, -1.0)
    z = Complex128(x + y * im)
    z == zero(z) || @test_approx_eq SpecialFunctions.hankelh1x(nu, z) SpecialFunctions.hankelh1(nu, z) * exp(-z * im)
    z == zero(z) || @test_approx_eq SpecialFunctions.hankelh2x(nu, z) SpecialFunctions.hankelh2(nu, z) * exp(z * im)
    (nu < 0 && z == zero(z)) || @test_approx_eq SpecialFunctions.besselix(nu, z) SpecialFunctions.besseli(nu, z) * exp(-abs(real(z)))
    (nu < 0 && z == zero(z)) || @test_approx_eq SpecialFunctions.besseljx(nu, z) SpecialFunctions.besselj(nu, z) * exp(-abs(imag(z)))
    z == zero(z) || @test_approx_eq SpecialFunctions.besselkx(nu, z) SpecialFunctions.besselk(nu, z) * exp(z)
    z == zero(z) || @test_approx_eq SpecialFunctions.besselyx(nu, z) SpecialFunctions.bessely(nu, z) * exp(-abs(imag(z)))
end
@test_throws SpecialFunctions.AmosException SpecialFunctions.hankelh1x(1, 0)
@test_throws SpecialFunctions.AmosException SpecialFunctions.hankelh2x(1, 0)
@test_throws SpecialFunctions.AmosException SpecialFunctions.besselix(-1, 0)
@test_throws SpecialFunctions.AmosException SpecialFunctions.besseljx(-1, 0)
@test SpecialFunctions.besselkx(1, 0) == Inf
@test_throws SpecialFunctions.AmosException SpecialFunctions.besselyx(1, 0)

# values from Abramowitz & Stegun, Table 10.11 (p475)
#   x     Ai(x)         Ai'(x)       Bi(x)        Bi'(x)       x     Ai(x)         Ai'(x)       Bi(x)        Bi'(x)
table10p11 = [
  0.00  0.35502_805  -0.25881_940  0.61492_663  0.44828_836  0.50  0.23169_361  -0.22491_053  0.85427_704  0.54457_256;
  0.01  0.35243_992  -0.25880_174  0.61940_962  0.44831_926  0.51  0.22945_031  -0.22374_617  0.85974_431  0.54890_049;
  0.02  0.34985_214  -0.25874_909  0.62389_322  0.44841_254  0.52  0.22721_872  -0.22257_027  0.86525_543  0.55334_239;
  0.03  0.34726_505  -0.25866_197  0.62837_808  0.44856_911  0.53  0.22499_894  -0.22138_322  0.87081_154  0.55789_959;
  0.04  0.34467_901  -0.25854_090  0.63286_482  0.44878_987  0.54  0.22279_109  -0.22018_541  0.87641_381  0.56257_345;
  0.05  0.34209_435  -0.25838_640  0.63735_409  0.44907_570  0.55  0.22059_527  -0.21897_720  0.88206_341  0.56736_532;
  0.06  0.33951_139  -0.25819_898  0.64184_655  0.44942_752  0.56  0.21841_158  -0.21775_898  0.88776_152  0.57227_662;
  0.07  0.33693_047  -0.25797_916  0.64634_286  0.44984_622  0.57  0.21624_012  -0.21653_112  0.89350_934  0.57730_873;
  0.08  0.33435_191  -0.25772_745  0.65084_370  0.45033_270  0.58  0.21408_099  -0.21529_397  0.89930_810  0.58246_311;
  0.09  0.33177_603  -0.25744_437  0.65534_975  0.45088_787  0.59  0.21193_427  -0.21404_790  0.90515_902  0.58774_120;
  0.10  0.32920_313  -0.25713_042  0.65986_169  0.45151_263  0.60  0.20980_006  -0.21279_326  0.91106_334  0.59314_448;
  0.11  0.32663_352  -0.25678_613  0.66438_023  0.45220_789  0.61  0.20767_844  -0.21153_041  0.91702_233  0.59867_447;
  0.12  0.32406_751  -0.25641_200  0.66890_609  0.45297_457  0.62  0.20556_948  -0.21025_970  0.92303_726  0.60433_267;
  0.13  0.32150_538  -0.25600_854  0.67343_997  0.45381_357  0.63  0.20347_327  -0.20898_146  0.92910_941  0.61012_064;
  0.14  0.31894_743  -0.25557_625  0.67798_260  0.45472_582  0.64  0.20138_987  -0.20769_605  0.93524_011  0.61603_997;
  0.15  0.31639_395  -0.25511_565  0.68253_473  0.45571_223  0.65  0.19931_937  -0.20640_378  0.94143_066  0.62209_226;
  0.16  0.31384_521  -0.25462_724  0.68709_709  0.45677_373  0.66  0.19726_182  -0.20510_500  0.94768_241  0.62827_912;
  0.17  0.31130_150  -0.25411_151  0.69167_046  0.45791_125  0.67  0.19521_729  -0.20380_004  0.95399_670  0.63460_222;
  0.18  0.30876_307  -0.25356_898  0.69625_558  0.45912_572  0.68  0.19318_584  -0.20248_920  0.96037_491  0.64106_324;
  0.19  0.30623_020  -0.25300_013  0.70085_323  0.46041_808  0.69  0.19116_752  -0.20117_281  0.96681_843  0.64766_389;
  0.20  0.30370_315  -0.25240_547  0.70546_420  0.46178_928  0.70  0.18916_240  -0.19985_119  0.97332_866  0.65440_592;
  0.21  0.30118_218  -0.25178_548  0.71008_928  0.46324_026  0.71  0.18717_052  -0.19852_464  0.97990_703  0.66129_109;
  0.22  0.29866_753  -0.25114_067  0.71472_927  0.46477_197  0.72  0.18519_192  -0.19719_347  0.98655_496  0.66832_121;
  0.23  0.29615_945  -0.25047_151  0.71938_499  0.46638_539  0.73  0.18322_666  -0.19585_798  0.99327_394  0.67549_810;
  0.24  0.29365_818  -0.24977_850  0.72405_726  0.46808_147  0.74  0.18127_478  -0.19451_846  1.00006_542  0.68282_363;
  0.25  0.29116_395  -0.24906_211  0.72874_690  0.46986_119  0.75  0.17933_631  -0.19317_521  1.00693_091  0.69029_970;
  0.26  0.28867_701  -0.24832_284  0.73345_477  0.47172_554  0.76  0.17741_128  -0.19182_851  1.01387_192  0.69792_824;
  0.27  0.28619_757  -0.24756_115  0.73818_170  0.47367_549  0.77  0.17549_975  -0.19047_865  1.02088_999  0.70571_121;
  0.28  0.28372_586  -0.24677_753  0.74292_857  0.47571_205  0.78  0.17360_172  -0.18912_591  1.02798_667  0.71365_062;
  0.29  0.28126_209  -0.24597_244  0.74769_624  0.47783_623  0.79  0.17171_724  -0.18777_055  1.03516_353  0.72174_849;
  0.30  0.27880_648  -0.24514_636  0.75248_559  0.48004_903  0.80  0.16984_632  -0.18641_286  1.04242_217  0.73000_690;
  0.31  0.27635_923  -0.24429_976  0.75729_752  0.48235_148  0.81  0.16798_899  -0.18505_310  1.04976_421  0.73842_795;
  0.32  0.27392_055  -0.24343_309  0.76213_292  0.48474_462  0.82  0.16614_526  -0.18369_153  1.05719_128  0.74701_380;
  0.33  0.27149_064  -0.24254_682  0.76699_272  0.48722_948  0.83  0.16431_516  -0.18232_840  1.06470_504  0.75576_663;
  0.34  0.26906_968  -0.24164_140  0.77187_782  0.48980_713  0.84  0.16249_870  -0.18096_398  1.07230_717  0.76468_865;
  0.35  0.26665_787  -0.24071_730  0.77678_917  0.49247_861  0.85  0.16069_588  -0.17959_851  1.07999_939  0.77378_215;
  0.36  0.26425_540  -0.23977_495  0.78172_770  0.49524_501  0.86  0.15890_673  -0.17823_223  1.08778_340  0.78304_942;
  0.37  0.26186_243  -0.23881_481  0.78669_439  0.49810_741  0.87  0.15713_124  -0.17686_539  1.09566_096  0.79249_282;
  0.38  0.25947_916  -0.23783_731  0.79169_018  0.50106_692  0.88  0.15536_942  -0.17549_823  1.10363_385  0.80211_473;
  0.39  0.25710_574  -0.23684_291  0.79671_605  0.50412_463  0.89  0.15362_128  -0.17413_097  1.11170_386  0.81191_759;
  0.40  0.25474_235  -0.23583_203  0.80177_300  0.50728_168  0.90  0.15188_680  -0.17276_384  1.11987_281  0.82190_389;
  0.41  0.25238_916  -0.23480_512  0.80686_202  0.51053_920  0.91  0.15016_600  -0.17139_708  1.12814_255  0.83207_615;
  0.42  0.25004_630  -0.23376_259  0.81198_412  0.51389_833  0.92  0.14845_886  -0.17003_090  1.13651_496  0.84243_695;
  0.43  0.24771_395  -0.23270_487  0.81714_033  0.51736_025  0.93  0.14676_538  -0.16866_551  1.14499_193  0.85298_891;
  0.44  0.24539_226  -0.23163_239  0.82233_167  0.52092_614  0.94  0.14508_555  -0.16730_113  1.15357_539  0.86373_470;
  0.45  0.24308_135  -0.23054_556  0.82755_920  0.52459_717  0.95  0.14341_935  -0.16593_797  1.16226_728  0.87467_704;
  0.46  0.24078_139  -0.22944_479  0.83282_397  0.52837_457  0.96  0.14176_678  -0.16457_623  1.17106_959  0.88581_871;
  0.47  0.23849_250  -0.22833_050  0.83812_705  0.53225_956  0.97  0.14012_782  -0.16321_611  1.17998_433  0.89716_253;
  0.48  0.23621_482  -0.22720_310  0.84346_952  0.53625_338  0.98  0.13850_245  -0.16185_781  1.18901_352  0.90871_137;
  0.49  0.23394_848  -0.22606_297  0.84885_248  0.54035_729  0.99  0.13689_066  -0.16050_153  1.19815_925  0.92046_818;
  0.50  0.23169_361  -0.22491_053  0.85427_704  0.54457_256  1.00  0.13529_242  -0.15914_744  1.20742_359  0.93243_593;
]

tol = 11e-9
for i = 1:size(table10p11,1)
    x   = table10p11[i,1]
    ai  = table10p11[i,2]
    aip = table10p11[i,3]
    bi  = table10p11[i,4]
    bip = table10p11[i,5]
    @test_approx_eq_eps ai SpecialFunctions.airyai(x) tol
    @test_approx_eq_eps aip SpecialFunctions.airyaiprime(x) tol
    @test_approx_eq_eps bi SpecialFunctions.airybi(x) tol
    @test_approx_eq_eps bip SpecialFunctions.airybiprime(x) tol
    x   = table10p11[i,6]
    ai  = table10p11[i,7]
    aip = table10p11[i,8]
    bi  = table10p11[i,9]
    bip = table10p11[i,10]
    @test_approx_eq_eps ai SpecialFunctions.airyai(x) tol
    @test_approx_eq_eps aip SpecialFunctions.airyaiprime(x) tol
    @test_approx_eq_eps bi SpecialFunctions.airybi(x) tol
    @test_approx_eq_eps bip SpecialFunctions.airybiprime(x) tol
end
