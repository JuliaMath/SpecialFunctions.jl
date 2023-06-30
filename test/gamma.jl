@testset "gamma and friends" begin
    @testset "digamma" begin
        @testset "$elty" for elty in (Float32, Float64)
            @test digamma(convert(elty, 9)) ≈ convert(elty, 2.140641477955609996536345)
            @test digamma(convert(elty, 2.5)) ≈ convert(elty, 0.7031566406452431872257)
            @test digamma(convert(elty, 0.1)) ≈ convert(elty, -10.42375494041107679516822)
            @test digamma(convert(elty, 7e-4)) ≈ convert(elty, -1429.147493371120205005198)
            @test digamma(convert(elty, 7e-5)) ≈ convert(elty, -14286.29138623969227538398)
            @test digamma(convert(elty, 7e-6)) ≈ convert(elty, -142857.7200612932791081972)
            @test digamma(convert(elty, 2e-6)) ≈ convert(elty, -500000.5772123750382073831)
            @test digamma(convert(elty, 1e-6)) ≈ convert(elty, -1000000.577214019968668068)
            @test digamma(convert(elty, 7e-7)) ≈ convert(elty, -1428572.005785942019703646)
            @test digamma(convert(elty, -0.5)) ≈ convert(elty, .03648997397857652055902367)
            @test digamma(convert(elty, -1.1)) ≈ convert(elty,  10.15416395914385769902271)

            @test digamma(convert(elty, 0.1)) ≈ convert(elty, -10.42375494041108)
            @test digamma(convert(elty, 1/2)) ≈ convert(elty, -γ - log(4))
            @test digamma(convert(elty, 1)) ≈ convert(elty, -γ)
            @test digamma(convert(elty, 2)) ≈ convert(elty, 1 - γ)
            @test digamma(convert(elty, 3)) ≈ convert(elty, 3/2 - γ)
            @test digamma(convert(elty, 4)) ≈ convert(elty, 11/6 - γ)
            @test digamma(convert(elty, 5)) ≈ convert(elty, 25/12 - γ)
            @test digamma(convert(elty, 10)) ≈ convert(elty, 7129/2520 - γ)
        end
    end

    @testset "trigamma" begin
        @testset "$elty" for elty in (Float32, Float64)
            @test trigamma(convert(elty, 0.1)) ≈ convert(elty, 101.433299150792758817)
            @test trigamma(convert(elty, 0.1)) ≈ convert(elty, 101.433299150792758817)
            @test trigamma(convert(elty, 1/2)) ≈ convert(elty, π^2/2)
            @test trigamma(convert(elty, 1)) ≈ convert(elty, π^2/6)
            @test trigamma(convert(elty, 2)) ≈ convert(elty, π^2/6 - 1)
            @test trigamma(convert(elty, 3)) ≈ convert(elty, π^2/6 - 5/4)
            @test trigamma(convert(elty, 4)) ≈ convert(elty, π^2/6 - 49/36)
            @test trigamma(convert(elty, 5)) ≈ convert(elty, π^2/6 - 205/144)
            @test trigamma(convert(elty, 10)) ≈ convert(elty, π^2/6 - 9778141/6350400)
        end
    end

    @testset "invdigamma" begin
        @testset "$elty" for elty in (Float32, Float64)
            for val in [0.001, 0.01, 0.1, 1.0, 10.0]
                @test abs(invdigamma(digamma(convert(elty, val))) - convert(elty, val)) < 1e-8
            end
        end
        @test abs(invdigamma(2)) == abs(invdigamma(2.))
    end

    @testset "polygamma" begin
        @test polygamma(20, 7.) ≈ -4.644616027240543262561198814998587152547
        @test polygamma(20, Float16(7.)) ≈ -4.644616027240543262561198814998587152547
    end

    @testset "eta" begin
        @test eta(1) ≈ log(2)
        @test eta(2) ≈ pi^2/12
        @test eta(Float32(2)) ≈ eta(2)
        @test eta(Complex{Float32}(2)) ≈ eta(2)
    end
    @testset "gamma, loggamma, logabsgamma (complex argument)" begin
        @test gamma.(Float64[1:25;]) ≈ gamma.(1:25)
        for elty in (Float32, Float64)
            @test gamma(convert(elty,1/2)) ≈ convert(elty,sqrt(π))
            @test gamma(convert(elty,-1/2)) ≈ convert(elty,-2sqrt(π))
            @test logabsgamma(convert(elty,-1/2))[1] ≈ convert(elty,log(abs(gamma(-1/2))))
        end
        @test loggamma(1.4+3.7im) ≈ -3.7094025330996841898 + 2.4568090502768651184im
        @test loggamma(1.4+3.7im) ≈ log(gamma(1.4+3.7im))
        @test loggamma(-4.2+0im) ≈ logabsgamma(-4.2)[1] - 5pi*im
        @test logfactorial(0) == logfactorial(1) == 0
        @test logfactorial(2) == loggamma(3)
        # Ensure that the domain of logfactorial matches that of factorial (issue #21318)
        @test_throws DomainError logfactorial(-3)
        @test_throws DomainError loggamma(-4.2)
        @test_throws MethodError logfactorial(1.0)
    end

    # values taken from Wolfram Alpha
    @testset "loggamma & logabsgamma test cases" begin
        @test loggamma(-300im) ≅ -473.17185074259241355733179182866544204963885920016823743 - 1410.3490664555822107569308046418321236643870840962522425im
        @test loggamma(3.099) ≅ loggamma(3.099+0im) ≅ 0.786413746900558058720665860178923603134125854451168869796
        @test loggamma(1.15) ≅ loggamma(1.15+0im) ≅ -0.06930620867104688224241731415650307100375642207340564554
        @test logabsgamma(0.89)[1] ≅ loggamma(0.89+0im) ≅ 0.074022173958081423702265889979810658434235008344573396963
        @test loggamma(0.91) ≅ loggamma(0.91+0im) ≅ 0.058922567623832379298241751183907077883592982094770449167
        @test loggamma(0.01) ≅ loggamma(0.01+0im) ≅ 4.599479878042021722513945411008748087261001413385289652419
        @test loggamma(-3.4-0.1im) ≅ -1.1733353322064779481049088558918957440847715003659143454 + 12.331465501247826842875586104415980094316268974671819281im
        @test loggamma(-13.4-0.1im) ≅ -22.457344044212827625152500315875095825738672314550695161 + 43.620560075982291551250251193743725687019009911713182478im
        @test loggamma(-13.4+0.0im) ≅ conj(loggamma(-13.4-0.0im)) ≅ -22.404285036964892794140985332811433245813398559439824988 - 43.982297150257105338477007365913040378760371591251481493im
        @test loggamma(-13.4+8im) ≅ -44.705388949497032519400131077242200763386790107166126534 - 22.208139404160647265446701539526205774669649081807864194im
        @test logabsgamma(1+exp2(-20))[1] ≅ loggamma(1+exp2(-20)+0im) ≅ -5.504750066148866790922434423491111098144565651836914e-7
        @test loggamma(1+exp2(-20)+exp2(-19)*im) ≅ -5.5047799872835333673947171235997541985495018556426e-7 - 1.1009485171695646421931605642091915847546979851020e-6im
        @test loggamma(-300+2im) ≅ -1419.3444991797240659656205813341478289311980525970715668 - 932.63768120761873747896802932133229201676713644684614785im
        @test loggamma(300+2im) ≅ 1409.19538972991765122115558155209493891138852121159064304 + 11.4042446282102624499071633666567192538600478241492492652im
        @test loggamma(1-6im) ≅ -7.6099596929506794519956058191621517065972094186427056304 - 5.5220531255147242228831899544009162055434670861483084103im
        @test loggamma(1-8im) ≅ -10.607711310314582247944321662794330955531402815576140186 - 9.4105083803116077524365029286332222345505790217656796587im
        @test loggamma(1+6.5im) ≅ conj(loggamma(1-6.5im)) ≅ -8.3553365025113595689887497963634069303427790125048113307 + 6.4392816159759833948112929018407660263228036491479825744im
        @test loggamma(1+1im) ≅ conj(loggamma(1-1im)) ≅ -0.6509231993018563388852168315039476650655087571397225919 - 0.3016403204675331978875316577968965406598997739437652369im
        @test loggamma(-pi*1e7 + 6im) ≅ -5.10911758892505772903279926621085326635236850347591e8 - 9.86959420047365966439199219724905597399295814979993e7im
        @test loggamma(-pi*1e7 + 8im) ≅ -5.10911765175690634449032797392631749405282045412624e8 - 9.86959074790854911974415722927761900209557190058925e7im
        @test loggamma(-pi*1e14 + 6im) ≅ -1.0172766411995621854526383224252727000270225301426e16 - 9.8696044010873714715264929863618267642124589569347e14im
        @test loggamma(-pi*1e14 + 8im) ≅ -1.0172766411995628137711690403794640541491261237341e16 - 9.8696044010867038531027376655349878694397362250037e14im
        @test loggamma(2.05 + 0.03im) ≅ conj(loggamma(2.05 - 0.03im)) ≅ 0.02165570938532611215664861849215838847758074239924127515 + 0.01363779084533034509857648574107935425251657080676603919im
        @test loggamma(2+exp2(-20)+exp2(-19)*im) ≅ 4.03197681916768997727833554471414212058404726357753e-7 + 8.06398296652953575754782349984315518297283664869951e-7im
    end

    @testset "loggamma for non-finite arguments" begin
        @test loggamma(Inf + 0im) === Inf + 0im
        @test loggamma(Inf - 0.0im) === Inf - 0.0im
        @test loggamma(Inf + 1im) === Inf + Inf*im
        @test loggamma(Inf - 1im) === Inf - Inf*im
        @test loggamma(-Inf + 0.0im) === -Inf - Inf*im
        @test loggamma(-Inf - 0.0im) === -Inf + Inf*im
        @test loggamma(Inf*im) === -Inf + Inf*im
        @test loggamma(-Inf*im) === -Inf - Inf*im
        @test loggamma(Inf + Inf*im) === loggamma(NaN + 0im) === loggamma(NaN*im) === NaN + NaN*im
    end

    @testset "BigFloat" begin
        # test cases (taken from WolframAlpha, computed to 78 digits ≈ 256 bits)
        @test loggamma(big"3.099") ≈ big"0.78641374690055805872066586017892360313412585445116886979672329071050823224651" rtol=1e-75
        @test loggamma(big"1.15") ≈ big"-0.06930620867104688224241731415650307100375642207340564554412494594148673455871" rtol=1e-75
        @test logabsgamma(big"0.89")[1] ≈ big"0.0740221739580814237022658899798106584342350083445733969634566129726762260738245" rtol=1e-75
        @test loggamma(big"0.91") ≈ big"0.0589225676238323792982417511839070778835929820947704491677379048793029707373113" rtol=1e-75
        @test loggamma(big"0.01") ≈ big"4.59947987804202172251394541100874808726100141338528965241917138771477998920321" rtol=1e-75
        @test loggamma(1 + exp2(big"-20.0")) ≈ big"-5.50475006614886679092243442349111109814456565183691425527816079744208067935466e-7" rtol=1e-75

        # consistency
        @test loggamma(big(3124.0)) == log(gamma(big(3124.0)))
        @test loggamma(big(3124.0)) ≈ loggamma(3124.0)
        @test logabsgamma(big(3124.0)) == (loggamma(big(3124.0)), 1)
        @test logabsgamma(big(3124.0))[1] ≈ logabsgamma(3124.0)[1]

        # promotions
        @test loggamma(big(3124)) == loggamma(big(3124.0))
        @test loggamma(big(3//2)) == loggamma(big(1.5))
        @test logabsgamma(big(3124)) == logabsgamma(big(3124.0))
        @test logabsgamma(big(3//2)) == logabsgamma(big(1.5))

        # negative values
        @test loggamma(big(-3.0)) == big(Inf)
        @test_throws DomainError loggamma(big(-2.76))

        # non-finite values
        @test isnan(loggamma(big(NaN)))
        @test isnan(logabsgamma(big(NaN))[1])
        @test loggamma(big(Inf)) == big(Inf)
        @test logabsgamma(big(Inf))[1] == big(Inf)
    end

    @testset "Other float types" begin
        struct NotAFloat <: AbstractFloat
        end
        let x = one(Float16)
            @test gamma(x) ≈ one(Float16)
            @test gamma(x) isa Float16
            @test loggamma(x) ≈ zero(Float16)
            @test loggamma(x) isa Float16
        end
        @test_throws MethodError gamma(NotAFloat())
        @test_throws MethodError logabsgamma(NotAFloat())
        @test_throws MethodError loggamma(NotAFloat())

        # https://github.com/JuliaMath/SpecialFunctions.jl/issues/339
        # https://github.com/JuliaMath/SpecialFunctions.jl/issues/233
        @test_throws MethodError gamma(complex(big(1.0)))
        @test_throws MethodError loggamma(complex(big(1.0)))
    end
end

@testset "zeta" begin
    @test zeta(0) ≈ -0.5
    @test zeta(2) ≈ pi^2/6
    @test zeta(Complex{Float32}(2)) ≈ zeta(2)
    @test zeta(4) ≈ pi^4/90
    @test zeta(1,Float16(2.)) ≈ zeta(1,2.)
    @test zeta(1.,Float16(2.)) ≈ zeta(1,2.)
    @test zeta(Float16(1.),Float16(2.)) ≈ zeta(1,2.)
    @test isnan(zeta(NaN))
    @test isnan(zeta(1.0e0))
    @test isnan(zeta(1.0f0))
    @test isnan(zeta(complex(0,Inf)))
    @test isnan(zeta(complex(-Inf,0)))
end

#(compared to Wolfram Alpha)
@testset "digamma, trigamma, polygamma & zeta" begin
    for x in -10.2:0.3456:50
        @test 1e-12 > relerr(digamma(x+0im), digamma(x))
    end
    @test digamma(7+0im) ≅ 1.872784335098467139393487909917597568957840664060076401194232
    @test digamma(7im) ≅ 1.94761433458434866917623737015561385331974500663251349960124 + 1.642224898223468048051567761191050945700191089100087841536im
    @test digamma(-3.2+0.1im) ≅ 4.65022505497781398615943030397508454861261537905047116427511+2.32676364843128349629415011622322040021960602904363963042380im
    @test trigamma(8+0im) ≅ 0.133137014694031425134546685920401606452509991909746283540546
    @test trigamma(8im) ≅ -0.0078125000000000000029194973110119898029284994355721719150 - 0.12467345030312762782439017882063360876391046513966063947im
    @test trigamma(-3.2+0.1im) ≅ 15.2073506449733631753218003030676132587307964766963426965699+15.7081038855113567966903832015076316497656334265029416039199im
    @test polygamma(2, 8.1+0im) ≅ -0.01723882695611191078960494454602091934457319791968308929600
    @test polygamma(30, 8.1+2im) ≅ -2722.8895150799704384107961215752996280795801958784600407589+6935.8508929338093162407666304759101854270641674671634631058im
    @test polygamma(3, 2.1+1im) ≅ 0.00083328137020421819513475400319288216246978855356531898998-0.27776110819632285785222411186352713789967528250214937861im
    @test 1e-11 > relerr(polygamma(3, -4.2 + 2im),-0.0037752884324358856340054736472407163991189965406070325067-0.018937868838708874282432870292420046797798431078848805822im)
    @test polygamma(13, 5.2 - 2im) ≅ 0.08087519202975913804697004241042171828113370070289754772448-0.2300264043021038366901951197725318713469156789541415899307im
    @test 1e-11 > relerr(polygamma(123, -47.2 + 0im), 5.7111648667225422758966364116222590509254011308116701029e291)
    @test zeta(4.1+0.3im, -3.2+0.1im) ≅ -281.34474134962502296077659347175501181994490498591796647 + 286.55601240093672668066037366170168712249413003222992205im
    @test zeta(4.1+0.3im, 3.2+0.1im) ≅ 0.0121197525131633219465301571139288562254218365173899270675-0.00687228692565614267981577154948499247518236888933925740902im
    @test zeta(4.1, 3.2+0.1im) ≅ 0.0137637451187986846516125754047084829556100290057521276517-0.00152194599531628234517456529686769063828217532350810111482im
    @test 1e-12 > relerr(zeta(1.0001, -4.5e2+3.2im), 10003.765660925877260544923069342257387254716966134385170 - 0.31956240712464746491659767831985629577542514145649468090im)
    @test zeta(3.1,-4.2) ≅ zeta(3.1,-4.2+0im) ≅ 149.7591329008219102939352965761913107060718455168339040295
    @test 1e-15 > relerr(zeta(3.1+0im,-4.2), zeta(3.1,-4.2+0im))
    @test zeta(3.1,4.2) ≅ 0.029938344862645948405021260567725078588893266227472565010234
    @test zeta(27, 3.1) ≅ 5.413318813037879056337862215066960774064332961282599376e-14
    @test zeta(27, 2) ≅ 7.4507117898354294919810041706041194547190318825658299932e-9
    @test 1e-12 > relerr(zeta(27, -105.3), 1.3113726525492708826840989036205762823329453315093955e14)
    @test polygamma(4, -3.1+Inf*im) == polygamma(4, 3.1+Inf*im) == 0
    @test polygamma(4, -0.0) == Inf == -polygamma(4, +0.0)
    @test zeta(4, +0.0) == zeta(4, -0.0) ≅ pi^4 / 90
    @test zeta(5, +0.0) == zeta(5, -0.0) ≅ 1.036927755143369926331365486457034168057080919501912811974
    @test zeta(Inf, 1.) == 1
    @test zeta(Inf, 2.) == 0
    @test isnan(zeta(NaN, 1.))
    @test isa([digamma(x) for x in [1.0]], Vector{Float64})
    @test isa([trigamma(x) for x in [1.0]], Vector{Float64})
    @test isa([polygamma(3,x) for x in [1.0]], Vector{Float64})
    @test zeta(2 + 1im, -1.1) ≅ zeta(2 + 1im, -1.1+0im) ≅ -64.580137707692178058665068045847533319237536295165484548 + 73.992688148809018073371913557697318846844796582012921247im
    @test polygamma(3,5) ≈ polygamma(3,5.)

    @test zeta(-3.0, 7.0) ≅ -52919/120
    @test zeta(-3.0, -7.0) ≅ 94081/120
    @test zeta(-3.1, 7.2) ≅ -587.457736596403704429103489502662574345388906240906317350719
    @test zeta(-3.1, -7.2) ≅ 1042.167459863862249173444363794330893294733001752715542569576
    @test zeta(-3.1, 7.0) ≅ -518.431785723446831868686653718848680989961581500352503093748
    @test zeta(-3.1, -7.0) ≅ 935.1284612957581823462429983411337864448020149908884596048161
    @test zeta(-3.1-0.1im, 7.2) ≅ -579.29752287650299181119859268736614017824853865655709516268 - 96.551907752211554484321948972741033127192063648337407683877im
    @test zeta(-3.1-0.1im, -7.2) ≅ 1025.17607931184231774568797674684390615895201417983173984531 + 185.732454778663400767583204948796029540252923367115805842138im
    @test zeta(-3.1-0.1im, 7.2 + 0.1im) ≅ -571.66133526455569807299410569274606007165253039948889085762 - 131.86744836357808604785415199791875369679879576524477540653im
    @test zeta(-3.1-0.1im, -7.2 + 0.1im) ≅ 1035.35760409421020754141207226034191979220047873089445768189 + 130.905870774271320475424492384335798304480814695778053731045im
    @test zeta(-3.1-0.1im, -7.0 + 0.1im) ≅ 929.546530292101383210555114424269079830017210969572819344670 + 113.646687807533854478778193456684618838875194573742062527301im
    @test zeta(-3.1, 7.2 + 0.1im) ≅ -586.61801005507638387063781112254388285799318636946559637115 - 36.148831292706044180986261734913443701649622026758378669700im
    @test zeta(-3.1, -7.2 + 0.1im) ≅ 1041.04241628770682295952302478199641560576378326778432301623 - 55.7154858634145071137760301929537184886497572057171143541058im
    @test zeta(-13.4, 4.1) ≅ -3.860040842156185186414774125656116135638705768861917e6
    @test zeta(3.2, -4) ≅ 2.317164896026427640718298719837102378116771112128525719078
    @test zeta(3.2, 0) ≅ 1.166773370984467020452550350896512026772734054324169010977
    @test zeta(-3.2+0.1im, 0.0) ≅ zeta(-3.2+0.1im, 0.0+0im) ≅ 0.0070547946138977701155565365569392198424378109226519905493 + 0.00076891821792430587745535285452496914239014050562476729610im
    @test zeta(-3.2, 0.0) ≅ zeta(-3.2, 0.0+0im) ≅ 0.007011972077091051091698102914884052994997144629191121056378

    @test 1e-14 > relerr(eta(1+1e-9), 0.693147180719814213126976796937244130533478392539154928250926)
    @test 1e-14 > relerr(eta(1+5e-3), 0.693945708117842473436705502427198307157819636785324430166786)
    @test 1e-13 > relerr(eta(1+7.1e-3), 0.694280602623782381522315484518617968911346216413679911124758)
    @test 1e-13 > relerr(eta(1+8.1e-3), 0.694439974969407464789106040237272613286958025383030083792151)
    @test 1e-13 > relerr(eta(1 - 2.1e-3 + 2e-3 * im), 0.69281144248566007063525513903467244218447562492555491581+0.00032001240133205689782368277733081683574922990400416791019im)
    @test 1e-13 > relerr(eta(1 + 5e-3 + 5e-3 * im), 0.69394652468453741050544512825906295778565788963009705146+0.00079771059614865948716292388790427833787298296229354721960im)
    @test 1e-12 > relerr(zeta(1e-3+1e-3im), -0.5009189365276307665899456585255302329444338284981610162-0.0009209468912269622649423786878087494828441941303691216750im)
    @test 1e-13 > relerr(zeta(1e-4 + 2e-4im), -0.5000918637469642920007659467492165281457662206388959645-0.0001838278317660822408234942825686513084009527096442173056im)

    # Issue #7169:
    @test 1e-13  > relerr(zeta(0 + 99.69im), 4.67192766128949471267133846066040655597942700322077493021802+3.89448062985266025394674304029984849370377607524207984092848im)
    @test 1e-12 > relerr(zeta(3 + 99.69im), 1.09996958148566565003471336713642736202442134876588828500-0.00948220959478852115901654819402390826992494044787958181148im)
    @test 1e-13  > relerr(zeta(-3 + 99.69im), 10332.6267578711852982128675093428012860119184786399673520976+13212.8641740351391796168658602382583730208014957452167440726im)
    @test 1e-13 > relerr(zeta(2 + 99.69im, 1.3), 0.41617652544777996034143623540420694985469543821307918291931-0.74199610821536326325073784018327392143031681111201859489991im)

    # issue #128
    @test 1e-13 > relerr(zeta(.4 + 453.0im), 5.595631794716693 - 4.994584420588448im)
    @test 1e-10 > relerr(zeta(.4 + 4053.0im), -0.1248993234383550+0.9195498409364987im)
    @test 1e-13 > relerr(zeta(.4 + 12.01im), 1.0233184799021265846512208845-0.8008078492939259287905322251im)
    @test zeta(.4 + 12.01im) == conj(zeta(.4 - 12.01im))

    # issue #420
    @test zeta(-2+13im) ≅ conj(zeta(-2-13im)) ≅ -0.30019019877262619754737023564024299182018857012958761814433485-5.5583626885487917197617298283836431070419020764882132809770386im
    @test zeta(-6+13im) ≅ conj(zeta(-6-13im)) ≅ 133.4764526350263089084083707864441932569167866714712267139316498-54.15465727586582149098585229287107039070546786014930791081909684im
    @test 1e-12 > relerr(zeta(-2+13im, 3), 2.3621038290867825837364823054-3.9497600485207119519185591345im)
    @test 1e-12 > relerr(zeta(-2-13im, 3), 2.3621038290867825837364823054+3.9497600485207119519185591345im)
end

@testset "logabsbinomial" begin
    @test logabsbinomial(10, -1) == (-Inf, 0.0)
    @test logabsbinomial(10, 11) == (-Inf, 0.0)
    @test logabsbinomial(10,  0) == ( 0.0, 1.0)
    @test logabsbinomial(10, 10) == ( 0.0, 1.0)

    @test logabsbinomial(10,  1)[1]   ≈ log(10)
    @test logabsbinomial(10,  1)[2]   == 1.0
    @test logabsbinomial(-6, 10)[1]   ≈ log(binomial(-6, 10))
    @test logabsbinomial(-6, 10)[2]   == 1.0
    @test logabsbinomial(-6, 11)[1]   ≈ log(abs(binomial(-6, 11)))
    @test logabsbinomial(-6, 11)[2]   == -1.0
    @test first.(logabsbinomial.(200, 0:200)) ≈ log.(binomial.(BigInt(200), (0:200)))
end

@testset "beta, logbeta, and logabsbeta" begin
    @test beta(3/2, 7/2) ≈ 5π/128
    @test beta(3, 5)     ≈ 1/105
    @test logbeta(5, 4)  ≈ log(beta(5,4))
    @test beta(5, 4)     ≈ beta(4,5)

    @testset "negative integer argument" begin
        @test beta(-2.0, 1.0)          ≈ -1/2    rtol=1e-14
        @test beta(-2.0, 2.0)          ≈  1/2    rtol=1e-14
        @test beta(-5.0, 2.0)          ≈  1/20   rtol=1e-14
        @test logabsbeta(-2.0, 2.0)[1] ≈ -log(2) rtol=1e-14
        @test beta(-2.0, -2.0)         == Inf
        @test logbeta(-2.0, -2.0)      == Inf
        @test beta(-2.0, 1.9)          == Inf
        @test logbeta(-2.0, 1.9)       == Inf
        @test beta(-2.0, 2.1)          == Inf
        @test logbeta(-2.0, 2.1)       == Inf
    end

    @testset "large difference in magnitude" begin
        @test beta(1e4, 1.5)          ≈     8.86193693673874630607029e-7  rtol=1e-14
        @test logabsbeta(1e4, 1.5)[1] ≈ log(8.86193693673874630607029e-7) rtol=1e-14
        @test beta(1e8, 0.5)          ≈     0.00017724538531210809        rtol=1e-14
        @test logabsbeta(1e8, 0.5)[1] ≈ log(0.00017724538531210809)       rtol=1e-14
    end

    @testset "BigFloat" begin
        @test beta(big(3), big(5)) ≈ inv(big(105))
        @test beta(big(3//2), big(7//2)) ≈ 5 * big(π) / 128
        @test beta(big(1e4), big(3//2)) ≈ 8.86193693673874630607029e-7 rtol=1e-14
        @test beta(big(1e8), big(1//2)) ≈ 0.00017724538531210809 rtol=1e-14

        @test logbeta(big(5), big(4)) ≈ log(beta(big(5), big(4)))
        @test logbeta(big(5.0), big(4)) ≈ log(beta(big(5), big(4)))
        @test logbeta(big(1e4), big(3//2)) ≈ log(8.86193693673874630607029e-7) rtol=1e-14
        @test logbeta(big(1e8), big(1//2)) ≈ log(0.00017724538531210809) rtol=1e-14

        @test logabsbeta(-big(2), big(2))[1] ≈ -log(big(2))
        @test logabsbeta(-big(2.0), big(2))[1] ≈ -log(big(2))
        @test logabsbeta(-big(2), big(2//1))[1] ≈ -log(big(2))
        @test logabsbeta(big(1e4), big(3//2))[1] ≈ log(8.86193693673874630607029e-7) rtol=1e-14
        @test logabsbeta(big(1e8), big(1//2))[1] ≈ log(0.00017724538531210809) rtol=1e-14
    end

    @test beta(-1/2, 3) ≈ beta(-1/2 + 0im, 3 + 0im) ≈    -16/3
    @test logabsbeta(-1/2, 3)[1]                    ≈ log(16/3)
    @test beta(Float32(5), Float32(4))             == beta(Float32(4), Float32(5))
    @test beta(3, 5)                                ≈ beta(3 + 0im, 5 + 0im)
    @test(beta(3.2 + 0.1im, 5.3 + 0.3im)            ≈ exp(logbeta(3.2 + 0.1im, 5.3 + 0.3im)) ≈
          0.00634645247782269506319336871208405439180447035257028310080 -
          0.00169495384841964531409376316336552555952269360134349446910im)
    @test beta(big(1.0), big(1.2))                  ≈ beta(1.0, 1.2) rtol=4*eps()
end
