@testset "incomplete beta" begin
    # Compared with Wolfram Engine's BetaRegularized
    @testset "a=0.5, b=0.5, x=$x" for (x, val) in zip(
        0.0:0.01:1.0,
        [0.e-25,0.06376856085851985,0.0903344706017331,0.1108246866044594,0.1281884336979499,
         0.1435662931287063,0.1575424254240443,0.170463428383199,0.1825548909938284,0.1939733680413566,
         0.2048327646991335,0.2152190255436425,0.2251989006184253,0.2348254690531804,0.2441417753461602,
         0.253183311106635,0.2619797608689092,0.2705562622683007,0.2789343361135709,0.2871325862574125,
         0.2951672353008665,0.3030525408934699,0.3108011236524005,0.3184242286319003,0.3259319361244763,
         0.3333333333333333,0.3406366554771171,0.3478494027630338,0.3549784381258894,0.3620300695016579,
         0.3690101195655454,0.375923985234415,0.3827766887550388,0.3895729218328954,0.3963170839725418,
         0.4030133159793217,0.4096655293982669,0.4162774325283954,0.42285255354095,0.4293942611422476,
         0.4359057831510251,0.4423902233029032,0.4488505765480783,0.4552897430704916,0.4617105412258534,
         0.4681157195707401,0.4745079681345231,0.48088992906936,0.4872642068002451,0.49363337778673,0.5,
         0.50636662221327,0.5127357931997549,0.51911007093064,0.5254920318654769,0.5318842804292599,
         0.5382894587741466,0.5447102569295084,0.5511494234519217,0.5576097766970968,0.5640942168489749,
         0.5706057388577524,0.57714744645905,0.5837225674716046,0.5903344706017331,0.5969866840206783,
         0.6036829160274582,0.6104270781671046,0.6172233112449612,0.624076014765585,0.6309898804344546,
         0.6379699304983421,0.6450215618741106,0.6521505972369662,0.6593633445228829,0.6666666666666667,
         0.6740680638755237,0.6815757713680997,0.6891988763475995,0.6969474591065301,0.7048327646991335,
         0.7128674137425875,0.7210656638864291,0.7294437377316993,0.7380202391310908,0.746816688893365,
         0.7558582246538398,0.7651745309468196,0.7748010993815747,0.7847809744563575,0.7951672353008665,
         0.8060266319586434,0.8174451090061716,0.829536571616801,0.8424575745759557,0.8564337068712937,
         0.8718115663020501,0.8891753133955406,0.9096655293982669,0.9362314391414802,1.])
        @test beta_inc(0.5, 0.5, x, 1.0 - x)[1] ≈ val rtol=1e-14
    end
    @testset "a=0.5, b=0.8, x=$x" for (x, val) in zip(
        0.0:0.01:1.0,
        [0.e-24,0.08704165448867688,0.1231783937242339,0.1509643905099393,0.1744376469592185,
         0.1951612417150757,0.2139363526723944,0.2312386644356416,0.247377962648835,0.2625693322193375,
         0.2769693157845956,0.2906960523171838,0.3038413001517716,0.3164780157848126,0.3286653449719896,
         0.3404520303383471,0.3518788093124823,0.3629801455598847,0.3737855072250243,0.3843203290282974,
         0.394606748829225,0.4046641800910133,0.4145097628356969,0.424158723204309,0.4336246632967352,
         0.4429197971428956,0.4520551445675527,0.4610406917926526,0.4698855255076289,0.4785979455868333,
         0.4871855604805254,0.4956553684394569,0.5040138270750072,0.5122669132520822,0.5204201749213043,
         0.5284787761920692,0.5364475367081124,0.5443309661970611,0.5521332949136762,0.5598585005745779,
         0.5675103322837338,0.5750923318679227,0.5826078529759856,0.5900600782419893,0.5974520347681598,
         0.6047866081468104,0.612066555210059,0.6192945156707836,0.6264730227971188,0.6336045132451182,
         0.6406913361594411,0.6477357616396057,0.6547399886591094,0.6617061525162714,0.6686363318887601,
         0.6755325555582559,0.6823968088674286,0.6892310399682792,0.6960371659188544,0.7028170786843475,
         0.7095726510986684,0.7163057428437263,0.7230182065060113,0.7297118937736915,0.7363886618425424,
         0.7430503801058253,0.7496989372120396,0.7563362485856976,0.7629642645204542,0.7695849789717461,
         0.7762004391984992,0.7828127564316306,0.7894241177826232,0.7960367996505097,0.8026531829430502,
         0.8092757705016214,0.8159072072147228,0.8225503034294581,0.8292080624343213,0.8358837130049195,
         0.8425807482983071,0.8493029727827216,0.8560545594443592,0.8628401202921954,0.869664794295266,
         0.8765343585075641,0.8834553705447895,0.8904353542429415,0.8974830460562918,0.9046087289748337,
         0.9118246961145984,0.9191459128120601,0.9265909945948355,0.9341837118759302,0.941955425134812,
         0.9499492886199519,0.9582281548957946,0.9668913290074657,0.9761171938969928,0.9863137214671616,1.])
        @test beta_inc(0.5, 0.8, x, 1.0 - x)[1] ≈ val rtol=1e-14
    end
    @testset "a=0.9, b=0.8, x=$x" for (x, val) in zip(
        0.0:0.01:1.0,
        [0.e-45,0.01287348375077395,0.0240457784957195,0.0346688685639659,0.04495791979284836,
         0.0550110492927223,0.06488448192684902,0.07461461051118363,0.08422688488805261,0.0937400884197129,
         0.1031686494864254,0.1125239993592263,0.1218154216032231,0.1310506102439909,0.1402360514251891,
         0.1493772930740891,0.1584791407604112,0.1675458033326384,0.1765810034326174,0.1855880628612014,
         0.1945699695603456,0.2035294309106093,0.2124689166764504,0.2213906940068579,0.2302968562599333,
         0.2391893469703827,0.2480699799571045,0.2569404563342797,0.2658023790171647,0.2746572651853306,
         0.2835065570691575,0.2923516313514671,0.3011938074192502,0.3100343546562188,0.3188744989322686,
         0.3277154284186073,0.3365582988355958,0.3454042382230052,0.3542543513084796,0.3631097235387736,
         0.3719714248292748,0.3808405130799951,0.3897180375002914,0.3986050417798133,0.4075025671393664,
         0.4164116552923807,0.4253333513453621,0.4342687066639933,0.4432187817303744,0.452184649016201,
         0.4611673958964416,0.4701681276282771,0.4791879704207046,0.4882280746212869,0.4972896180480853,
         0.5063738094968672,0.5154818924563003,0.5246151490670925,0.5337749043650042,0.5429625308524718,
         0.5521794534493672,0.5614271548803832,0.5707071815648901,0.5800211500851545,0.5893707543209051,
         0.5987577733528309,0.608184080255283,0.6176516519199679,0.6271625800787011,0.6367190837255744,
         0.6463235231787504,0.6559784160716227,0.6656864556250031,0.6754505316299672,0.6852737546699122,
         0.695159484236884,0.7051113615604331,0.7151333481797695,0.7252297715695021,0.7354053795009442,
         0.74566540532109,0.7560156470130347,0.7664625638438512,0.7770133957285906,0.7876763123275041,
         0.7984606016409664,0.8093769119481182,0.8204375671426882,0.8316569852132864,0.8430522452167303,
         0.8546438740859578,0.8664569696979065,0.8785228586162501,0.8908816447296843,0.9035863305210574,
         0.9167099250794181,0.9303587974356043,0.9447009522476523,0.9600378941965647,0.977058098680342,1.])
        @test beta_inc(0.9, 0.8, x, 1.0 - x)[1] ≈ val rtol=1e-14
    end
    @testset "a=80.9, b=0.8, x=$x" for (x, val) in zip(
        0.0:0.01:1.0,
        [0.e-4037,5.659850188208874e-163,1.279388616124076e-138,2.257693896707542e-124,
         2.897994605366003e-114,2.009164686416781e-106,5.124902812310911e-100,1.338436546592981e-94,
         6.592846323627736e-90,9.083404801612933e-86,4.581171164611029e-82,1.024712535445972e-78,
         1.171290435847841e-75,7.618998306182543e-73,3.066399996029232e-70,8.159340151736502e-68,
         1.514230271101736e-65,2.047520095589173e-63,2.091671353695881e-61,1.664040903808287e-59,
         1.057759838815721e-57,5.49132667443713e-56,2.372578688736847e-54,8.671833527840728e-53,
         2.719793317701389e-51,7.411825180712226e-50,1.774472557976592e-48,3.768914466846032e-47,
         7.163917757923168e-46,1.228169152371666e-44,1.912408289817481e-43,2.721825167444e-42,
         3.561050588818447e-41,4.305131971329379e-40,4.832096799800394e-39,5.057055939273883e-38,
         4.954351648811546e-37,4.560116807292394e-36,3.956512105211734e-35,3.245869269592197e-34,
         2.525020123177697e-33,1.867480042455247e-32,1.316327433960537e-31,8.862853181709876e-31,
         5.712179225668737e-30,3.531051129157696e-29,2.097389394183402e-28,1.199152837592699e-27,
         6.609838550976534e-27,3.517914083771492e-26,1.810402708594319e-25,9.020730533511399e-25,
         4.357446464530653e-24,2.042971920816074e-23,9.307222587505059e-23,4.124448596487034e-22,
         1.779654898083663e-21,7.4841565165423e-21,3.070285970414923e-20,1.22974643568586e-19,
         4.812895987534055e-19,1.841991715270897e-18,6.898889614425955e-18,2.530381730672556e-17,
         9.094925675449364e-17,3.205503996194959e-16,1.108518744159741e-15,3.763511868841169e-15,
         1.255134520239131e-14,4.114026783167532e-14,1.326012138884187e-13,4.204809280909274e-13,
         1.312411918693827e-12,4.033834438688928e-12,1.221467371969832e-11,3.645413830024844e-11,
         1.072734311570093e-10,3.113805924559142e-10,8.918947491417688e-10,2.521869395260185e-9,
         7.041721770169991e-9,1.942404070962948e-8,5.294920235920981e-8,1.426893004118986e-7,
         3.802663384681951e-7,1.002539998693497e-6,2.615699030509946e-6,6.756237677059152e-6,
         0.00001728297704388556,0.00004380303345697832,0.0001100404885846187,0.0002741427949192357,
         0.0006776797945177599,0.001663395040736671,0.004057633817143372,0.009848744173397283,
         0.02382855932322502,0.05763883015096727,0.1402012176136536,0.348133320555652,1.])
        @test beta_inc(80.9, 0.8, x, 1.0 - x)[1] ≈ val rtol=1e-13
    end
    @testset "a=1.7, b=10.5, x=$x" for (x, val) in zip(
        0.0:0.01:1.0,
        [0.e-83,0.01396204974701841,0.04273844739996358,0.0802355082453466,0.1233212707222279,
         0.169883883601507,0.2183876645881966,0.2676795794438451,0.3168824713253145,0.3653277197583937,
         0.4125089764419372,0.4580484843629353,0.5016715394987723,0.5431865527277862,0.5824691495850261,
         0.6194492913831362,0.654100724564198,0.6864322667088425,0.7164805688276856,0.7443040822000412,
         0.7699780198879489,0.7935901474806426,0.8152372703488115,0.8350223093424488,0.8530518758192943,
         0.8694342717267099,0.8842778522619746,0.8976896981615276,0.9097745524576825,0.9206339829823968,
         0.9303657372794607,0.9390632611240216,0.9468153557047174,0.9537059518240998,0.9598139823148841,
         0.9652133363295071,0.9699728812994702,0.974156540228006,0.9778234136145069,0.9810279367445075,
         0.9838200643417606,0.9862454756915709,0.9883457943259302,0.9901588172271333,0.9917187492710794,
         0.993056439306095,0.9941996148579756,0.9951731129758233,0.9959991051938595,0.9966973149884605,
         0.997285226463178,0.997778283302743,0.9981900773047402,0.9985325260289722,0.998816039303301,
         0.9990496744943506,0.9992412805949365,0.9993976313002439,0.9995245473440961,0.9996270084474357,
         0.9997092552954444,0.9997748820094408,0.9998269196165406,0.9998679110455995,0.9998999781936217,
         0.9999248816138942,0.9999440733768028,0.9999587436476631,0.9999698615139476,0.9999782105779107,
         0.9999844198105916,0.9999889901402758,0.9999923172233462,0.9999947108186769,0.9999964111588281,
         0.9999976026827751,0.9999984254661612,0.9999989846564765,0.9999993581924477,0.9999996030595673,
         0.9999997603073096,0.9999998590283993,0.9999999194766476,0.9999999554775051,0.9999999762646713,
         0.9999999878569335,0.9999999940719059,0.9999999972575227,0.9999999988080026,0.999999999518498,
         0.999999999821743,0.999999999940621,0.99999999998264,0.999999999995699,0.999999999999142,
         0.999999999999873,0.999999999999988,0.999999999999999,1.,1.,1.])
        @test beta_inc(1.7, 10.5, x, 1.0 - x)[1] ≈ val rtol=1e-14
    end
    @testset "a=100.5, b=100.5, x=$x" for (x, val) in zip(
        0.0:0.01:1.0,
        [0.e-4966,3.356014382599029e-143,2.21351837164392e-113,4.014552683228967e-96,
         5.21091312517625e-84,1.020389199029434e-74,3.267737964647878e-67,6.101573725809464e-61,
         1.418580292320946e-55,6.698584764694289e-51,8.971972641791594e-47,4.325143541605987e-43,
         8.94664782857427e-40,9.074256619313614e-37,5.005701167291761e-34,1.630336905229558e-31,
         3.349388657777907e-29,4.581338772370011e-27,4.362832768445309e-25,3.002816845923933e-23,
         1.541642496405598e-21,6.064729213175483e-20,1.870837362242701e-18,4.616584448897111e-17,
         9.27255359484341e-16,1.539075197361062e-14,2.139403432844863e-13,2.520059271944931e-12,
         2.541836676624415e-11,2.215842537778625e-10,1.683408890736215e-9,1.122874450435998e-8,
         6.620173165982372e-8,3.470763841511403e-7,1.626934118513085e-6,6.852626959254337e-6,
         0.00002605241980849797,0.00008977130313754974,0.0002814332703651759,0.0008055416902949821,
         0.002112019781014849,0.00508797075711991,0.01129547597461094,0.0231746884205069,
         0.04406489065169599,0.07787036008215303,0.1282703308247337,0.1975636232169813,0.2854862534447546,
         0.3885022505229126,0.5,0.6114977494770874,0.7145137465552454,0.8024363767830187,0.8717296691752663,
         0.922129639917847,0.955935109348304,0.9768253115794931,0.9887045240253891,0.9949120292428801,
         0.9978879802189852,0.999194458309705,0.9997185667296348,0.9999102286968625,0.9999739475801915,
         0.9999931473730407,0.9999983730658815,0.9999996529236158,0.9999999337982683,0.9999999887712555,
         0.9999999983165911,0.999999999778416,0.999999999974582,0.99999999999748,0.999999999999786,
         0.999999999999985,0.999999999999999,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,
         1.,1.,1.])
        @test beta_inc(100.5, 100.5, x, 1.0 - x)[1] ≈ val rtol=1e-13
    end
    @testset "a=1.5, b=1e-20, x=$x" for (x, val) in zip(
        0.0:0.01:1.0,
        [0.e-95,6.706954621511613e-24,1.908573860989784e-23,3.527823559800854e-23,
         5.465108108164382e-23,7.685476511624814e-23,1.016654428420256e-22,1.289283267454361e-22,
         1.585315312526521e-22,1.903920840622343e-22,2.244476844084102e-22,2.606514996691386e-22,
         2.989687056036453e-22,3.393740920821702e-22,3.818503532521727e-22,4.26386835711828e-22,
         4.729786038720361e-22,5.21625731553983e-22,5.72332759176175e-22,6.251082749838933e-22,
         6.799645911929102e-22,7.369174942192775e-22,7.95986053851056e-22,8.571924801968895e-22,
         9.205620200883313e-22,9.861228866810969e-22,1.053906217532698e-21,1.12394605758909e-21,
         1.196279364399263e-21,1.270946033569195e-21,1.347988943018207e-21,1.427454015049574e-21,
         1.50939029562088e-21,1.593850050518841e-21,1.680888878423608e-21,1.77056584110116e-21,
         1.862943611198906e-21,1.958088638346952e-21,2.056071334492481e-21,2.156966279623312e-21,
         2.260852449274566e-21,2.367813465464866e-21,2.477937872981194e-21,2.591319443230187e-21,
         2.708057508204485e-21,2.828257327482878e-21,2.952030491600198e-21,3.079495365596476e-21,
         3.210777577095268e-21,3.346010553881064e-21,3.48533611665991e-21,3.628905133513816e-21,
         3.776878243519303e-21,3.929426658119426e-21,4.086733050147924e-21,4.248992541941733e-21,
         4.416413805790527e-21,4.589220292116253e-21,4.767651603322342e-21,4.951965034288911e-21,
         5.142437304125938e-21,5.339366508167657e-21,5.543074324470629e-21,5.753908515483088e-21,
         5.972245773362194e-21,6.19849496698491e-21,6.433100860486935e-21,6.67654838776915e-21,
         6.929367585613005e-21,7.192139310873532e-21,7.465501896037407e-21,7.75015893407293e-21,
         8.046888430431106e-21,8.356553620653897e-21,8.680115830952776e-21,9.018649862807561e-21,
         9.373362520215628e-21,9.745615082692113e-21,1.013695077726479e-20,1.05491286460553e-20,
         1.098416568357789e-20,1.14443897916644e-20,1.193250706559674e-20,1.245168833304345e-20,
         1.300568195962486e-20,1.359896310974523e-20,1.423693458276194e-20,1.492620220780857e-20,
         1.567496069508049e-20,1.649354777878116e-20,1.739526322363106e-20,1.839762082220938e-20,
         1.952434090962006e-20,2.080868037972059e-20,2.229934708861127e-20,2.407185551639959e-20,
         2.625271544895813e-20,2.907909078985462e-20,3.308342534794925e-20,3.996470818039522e-20,1.])
        @test beta_inc(1.5, 1e-20, x, 1.0 - x)[1] ≈ val rtol=1e-14
    end
    @testset "a=1e-20, b=0.99, x=$x" for (x, val) in zip(
        0.0:0.01:1.0,
        [1.0,4.621640571507914e-20,3.928392625456373e-20,3.522826232996826e-20,3.235042350168417e-20,
         3.011796455118398e-20,2.82937201372135e-20,2.675117900735502e-20,2.541482518522368e-20,
         2.423594928780016e-20,2.318129286271936e-20,2.222713398379447e-20,2.135595723375987e-20,
         2.05544611885328e-20,1.981230641875294e-20,1.912129648210598e-20,1.847482377915803e-20,
         1.786748370079641e-20,1.729479923207514e-20,1.675302011465343e-20,1.623897358450671e-20,
         1.574995156483639e-20,1.528362412553572e-20,1.483797219538297e-20,1.441123460536074e-20,
         1.400186594951488e-20,1.360850271543896e-20,1.322993581024732e-20,1.286508808543862e-20,
         1.251299580739944e-20,1.217279327043997e-20,1.184369993372369e-20,1.152500960103253e-20,
         1.121608126599496e-20,1.091633132431086e-20,1.062522691510309e-20,1.034228020045673e-20,
         1.006704342884576e-20,9.799104656964201e-21,9.538084027304864e-21,9.283630517029258e-21,
         9.035419088275798e-21,8.793148181840232e-21,8.556537505729393e-21,8.325326077896178e-21,
         8.099270488866332e-21,7.878143355244009e-21,7.661731939451283e-21,7.449836914688898e-21,
         7.242271257138218e-21,7.03885924996661e-21,6.839435585837694e-21,6.643844556434574e-21,
         6.451939319035191e-21,6.263581231480807e-21,6.078639247988984e-21,5.896989369212464e-21,
         5.718514140760424e-21,5.543102195099824e-21,5.370647832359387e-21,5.201050636081688e-21,
         5.034215120421851e-21,4.870050405684712e-21,4.708469919434276e-21,4.549391120707029e-21,
         4.392735245120113e-21,4.238427068891405e-21,4.086394689985318e-21,3.936569324769163e-21,
         3.788885118712737e-21,3.643278969790933e-21,3.499690363356917e-21,3.358061217343224e-21,
         3.21833573672026e-21,3.080460276196264e-21,2.944383210178798e-21,2.810054809033594e-21,
         2.677427120668968e-21,2.546453856438284e-21,2.41709028028191e-21,2.289293099913057e-21,
         2.163020358672793e-21,2.038231326414587e-21,1.914886387391528e-21,1.79294692255393e-21,
         1.67237518283336e-21,1.553134148749398e-21,1.435187369793683e-21,1.318498774123754e-21,
         1.20303243443583e-21,1.088752268197901e-21,9.756216372555088e-22,8.63602788212576e-22,
         7.526560302762232e-22,6.427384567004279e-22,5.338018164784628e-22,4.25788652801385e-22,
         3.186244273953916e-22,2.121983888512568e-22,1.063002842170579e-22,0])
        @test beta_inc(1e-20, 0.99, x, 1.0 - x)[2] ≈ val rtol=1e-14
    end

    @test beta_inc(1.5, 200.5,  0.07,0.93)[1] ≈ 0.99999790408564
    @test beta_inc(1e-20, 0.000001, 0.2)[2] ≈ 1.0000013862929421e-14

    # test promotions and return types
    for T in (Float16, Float32, Float64)
        for x in (T(0.1), rand(T))
            for a in (1, randexp(T)), b in (1, randexp(T))
                @test beta_inc(a, b, x) isa Tuple{T,T}
                @test T.(beta_inc(a, b, x, 1 - Float64(x))) === beta_inc(a, b, x)
            end
        end
    end
    a = randexp()
    b = randexp(Float32)
    x = rand(Float32)
    @test beta_inc(a, b, x) isa Tuple{Float64,Float64}
    @test beta_inc(a, b, x) == beta_inc(a, Float64(b), Float64(x))
    @test beta_inc(a, b, x, 1 - x) isa Tuple{Float64,Float64}
    @test beta_inc(a, b, x, 1 - x) == beta_inc(a, Float64(b), Float64(x), 1 - x)

    @test SpecialFunctions.loggammadiv(13.89, 21.0001) ≈ log(gamma(big(21.0001))/gamma(big(21.0001)+big(13.89)))
    @test SpecialFunctions.loggammadiv(big(13.89), big(21.0001)) ≈ log(gamma(big(21.0001))/gamma(big(21.0001)+big(13.89)))
    @test SpecialFunctions.stirling_corr(11.99, 100.1) ≈ SpecialFunctions.stirling_error(11.99) + SpecialFunctions.stirling_error(100.1) - SpecialFunctions.stirling_error(11.99 + 100.1)

    @testset "Issue 334. Underflow without erroring" begin
        @test beta_inc(0.1, 4000, 0.2) == (1.0, 0.0)
        @test beta_inc(4000, 0.1, 0.2) == (0.0, 1.0)
    end
end

@testset "inverse of incomplete beta" begin
    f(a,b,p) = beta_inc_inv(a,b,p)[1]
    @test f(.5, .5, 0.6376856085851985E-01)    ≈ 0.01   rtol=8eps()
    @test f(.5, .5, 0.20483276469913355)       ≈ 0.1    rtol=8eps()
    @test f(.5, .5, 1.0000)                    ≈ 1.0000 rtol=8eps()
    @test f(1.0, .5, 0.0)                     == 0.0
    @test f(1.0, .5, 0.5012562893380045E-02)   ≈ 0.01   rtol=8eps()
    @test f(1.0, .5, 0.5131670194948620E-01)   ≈ 0.1    rtol=8eps()
    @test f(1.0, .5, 0.2928932188134525)       ≈ 0.5    rtol=8eps()
    @test f(1.0, 1.0, .5)                      ≈ 0.5    rtol=8eps()
    @test f(2.0, 2.0, .028)                    ≈ 0.1    rtol=8eps()
    @test f(2.0, 2.0, 0.104)                   ≈ 0.2    rtol=8eps()
    @test f(2.0, 2.0, .216)                    ≈ 0.3    rtol=8eps()
    @test f(2.0, 2.0, .352)                    ≈ 0.4    rtol=8eps()
    @test f(2.0, 2.0, .5)                      ≈ 0.5    rtol=8eps()
    @test f(2.0, 2.0, 0.648)                   ≈ 0.6    rtol=8eps()
    @test f(2.0, 2.0, 0.784)                   ≈ 0.7    rtol=8eps()
    @test f(2.0, 2.0, 0.896)                   ≈ 0.8    rtol=8eps()
    @test f(2.0, 2.0, .972)                    ≈ 0.9    rtol=8eps()
    @test f(5.5, 5.0, 0.4361908850559777)      ≈ 0.5    rtol=8eps()
    @test f(10.0, .5, 0.1516409096347099)      ≈ 0.9    rtol=8eps()
    @test f(10.0, 5.0, 0.8978271484375000E-01) ≈ 0.5    rtol=8eps()
    @test f(10.0, 5.0, 1.00)                   ≈ 1.0    rtol=8eps()
    @test f(10.0, 10.0, .5)                    ≈ 0.5    rtol=8eps()
    @test f(20.0, 5.0, 0.4598773297575791)     ≈ 0.8    rtol=8eps()
    @test f(20.0, 10.0, 0.2146816102371739)    ≈ 0.6    rtol=8eps()
    @test f(20.0, 10.0, 0.9507364826957875)    ≈ 0.8    rtol=8eps()
    @test f(20.0, 20.0, .5)                    ≈ 0.5    rtol=8eps()
    @test f(20.0, 20.0, 0.8979413687105918)    ≈ 0.6    rtol=8eps()
    @test f(30.0, 10.0, 0.2241297491808366)    ≈ 0.7    rtol=8eps()
    @test f(30.0, 10.0, 0.7586405487192086)    ≈ 0.8    rtol=8eps()
    @test f(40.0, 20.0, 0.7001783247477069)    ≈ 0.7    rtol=8eps()
    @test f(1.0, 0.5, 0.5131670194948620E-01)  ≈ 0.1    rtol=8eps()
    @test f(1.0, 0.5, 0.1055728090000841)      ≈ 0.2    rtol=8eps()
    @test f(1.0, 0.5, 0.1633399734659245)      ≈ 0.3    rtol=8eps()
    @test f(1.0, 0.5, 0.2254033307585166)      ≈ 0.4    rtol=8eps()
    @test f(1.0, 2.0, .36)                     ≈ 0.2    rtol=8eps()
    @test f(1.0, 3.0, .488)                    ≈ 0.2    rtol=8eps()
    @test f(1.0, 4.0, .5904)                   ≈ 0.2    rtol=8eps()
    @test f(1.0, 5.0, .67232)                  ≈ 0.2    rtol=8eps()
    @test f(2.0, 2.0, .216)                    ≈ 0.3    rtol=8eps()
    @test f(3.0, 2.0, 0.837e-1)                ≈ 0.3    rtol=8eps()
    @test f(4.0, 2.0, 0.3078e-1)               ≈ 0.3    rtol=8eps()
    @test f(5.0, 2.0, 0.10935e-1)              ≈ 0.3    rtol=8eps()
    @test f(1.30625000, 11.75620000, 0.9188846846205182)  ≈ 0.225609 rtol=8eps()
    @test f(1.30625000, 11.75620000, 0.21053116418502513) ≈ 0.033557 rtol=8eps()
    @test f(1.30625000, 11.75620000, 0.18241165418408148) ≈ 0.029522 rtol=8eps()
    @test f(1000.0, 2.0, 9.0797754e-317) ≈ 0.48 # This one is a bit slow (but also hard)
    @test f(1.5, 5.0, beta_inc(1.5, 5.0, 0.2142857142857142)[1])[1] ≈ 0.2142857142857142 rtol=8eps()

    for T in (Float16, Float32, Float64)
        p = rand(T)
        for a in (1, randexp(T)), b in (1, randexp(T))
            @test beta_inc_inv(a, b, p) isa Tuple{T,T}
            @test beta_inc_inv(a, b, p, 1 - p) === beta_inc_inv(a, b, p)
        end
    end

    @testset "Avoid NaN when p=q=1" begin
        @test first(beta_inc_inv(1.0, 1.0, 1e-21)) ≈ 1e-21
        @test beta_inc_inv(1.0e30, 1.0, 0.49) == (1.0, 0.0)
    end

    @testset "Avoid infinite loops" begin
        # See https://github.com/JuliaStats/StatsFuns.jl/issues/133#issuecomment-1069602721
        y = 2.0e-280
        @test beta_inc(2.0, 1.0, beta_inc_inv(2.0, 1.0, y, 1.0)[1])[1] ≈ y
        # See https://github.com/JuliaStats/GLM.jl/issues/538#issuecomment-1603447448
        @test isequal(beta_inc(6.0, 112.5, NaN), (NaN, NaN))
    end

    @testset "StatsFuns#145" begin
        y = 0.92
        @test beta_inc_inv(0.01, 0.1, y)[1] ≈ 0.7803014210919872
        @test beta_inc(0.01, 0.1, beta_inc_inv(0.01, 0.1, y)[1])[1] ≈ y
    end
end
