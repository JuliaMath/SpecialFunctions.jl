@testset "airy" begin
    @test_throws AmosException airyai(200im)
    @test_throws AmosException airybi(200)

    for T in [Float16, Float32, Float64,Complex{Float16}, Complex{Float32},Complex{Float64}]
        @test airyai(T(1.8)) ≈ 0.0470362168668458052247
        @test airyaiprime(T(1.8)) ≈ -0.0685247801186109345638
        @test airybi(T(1.8)) ≈ 2.595869356743906290060
        @test airybiprime(T(1.8)) ≈ 2.98554005084659907283
    end
    for T in [Complex{Float16}, Complex{Float32}, Complex{Float64}]
        z = convert(T,1.8 + 1.0im)
        @test airyaix(z) ≈ airyai(z) * exp(2/3 * z * sqrt(z))
        @test airyaiprimex(z) ≈ airyaiprime(z) * exp(2/3 * z * sqrt(z))
        @test airybix(z) ≈ airybi(z) * exp(-abs(real(2/3 * z * sqrt(z))))
        @test airybiprimex(z) ≈ airybiprime(z) * exp(-abs(real(2/3 * z * sqrt(z))))
    end
    @test_throws MethodError airyai(complex(big(1.0)))

    for x = -3:3
        @test airyai(x) ≈ airyai(complex(x))
        @test airyaiprime(x) ≈ airyaiprime(complex(x))
        @test airybi(x) ≈ airybi(complex(x))
        @test airybiprime(x) ≈ airybiprime(complex(x))
        if x >= 0
            @test airyaix(x) ≈ airyaix(complex(x))
            @test airyaiprimex(x) ≈ airyaiprimex(complex(x))
        else
            @test_throws DomainError airyaix(x)
            @test_throws DomainError airyaiprimex(x)
        end
        @test airybix(x) ≈ airybix(complex(x))
        @test airybiprimex(x) ≈ airybiprimex(complex(x))
    end
end

@testset "bessel functions" begin
    bessel_funcs = [(bessely0, bessely1, bessely), (besselj0, besselj1, besselj)]
    @testset "$z, $o" for (z, o, f) in bessel_funcs
        @test z(Float32(2.0)) ≈ z(Float64(2.0))
        @test o(Float32(2.0)) ≈ o(Float64(2.0))
        @test z(Float16(2.0)) ≈ z(Float64(2.0))
        @test o(Float16(2.0)) ≈ o(Float64(2.0))
        @test z(2) ≈ z(2.0)
        @test o(2) ≈ o(2.0)
        @test z(2.0 + im) ≈ f(0, 2.0 + im)
        @test o(2.0 + im) ≈ f(1, 2.0 + im)
    end
    @testset "besselj error throwing" begin
        @test_throws MethodError besselj(1.2,big(1.0))
        @test_throws MethodError besselj(1,complex(big(1.0)))
        @test_throws MethodError besseljx(1,big(1.0))
        @test_throws MethodError besseljx(1,complex(big(1.0)))
    end
    @testset "besselh" begin
        true_h133 = 0.30906272225525164362 - 0.53854161610503161800im
        @test besselh(3,1,3) ≈ true_h133
        @test besselh(-3,1,3) ≈ -true_h133
        @test besselh(Float32(3),1,Float32(3)) ≈ true_h133
        @test besselh(Float16(3),1,Float16(3)) ≈ true_h133
        @test besselh(3,2,3) ≈ conj(true_h133)
        @test besselh(-3,2,3) ≈ -conj(true_h133)
        @test besselh(1,1,0) ≈ complex(0,-Inf)
        @test besselh(1,2,0) ≈ complex(0,Inf)
        @testset "Error throwing" begin
            @test_throws AmosException besselh(1,5,0)
            @test_throws MethodError besselh(1,big(1.0))
            @test_throws MethodError besselh(1,complex(big(1.0)))
            @test_throws MethodError besselhx(1,big(1.0))
            @test_throws MethodError besselhx(1,complex(big(1.0)))
        end
    end
    @testset "besseli" begin
        true_i33 = 0.95975362949600785698
        @test besseli(3,3) ≈ true_i33
        @test besseli(-3,3) ≈ true_i33
        @test besseli(3,-3) ≈ -true_i33
        @test besseli(-3,-3) ≈ -true_i33
        @test besseli(Float32(-3),Complex{Float32}(-3,0)) ≈ -true_i33
        @test besseli(Float16(-3),Complex{Float16}(-3,0)) ≈ -true_i33
        true_im3p1_3 = 0.84371226532586351965
        @test besseli(-3.1,3) ≈ true_im3p1_3
        for i in [-5 -3 -1 1 3 5]
            @test besseli(i,0) == 0.0
            @test besseli(i,Float32(0)) == 0
            @test besseli(i,Complex{Float32}(0)) == 0
            @test besseli(i,Float16(0)) == 0
            @test besseli(i,Complex{Float16}(0)) == 0
        end
        @testset "Error throwing" begin
            @test_throws AmosException besseli(1,1000)
            @test_throws DomainError besseli(0.4,-1.0)
            @test_throws MethodError besseli(1,big(1.0))
            @test_throws MethodError besseli(1,complex(big(1.0)))
            @test_throws MethodError besselix(1,big(1.0))
            @test_throws MethodError besselix(1,complex(big(1.0)))
        end
    end
    @testset "besselj" begin
        @test besselj(0,0) == 1
        for i in [-5 -3 -1 1 3 5]
            @test besselj(i,0) == 0
            @test besselj(i,Float32(0)) == 0
            @test besselj(i,Complex{Float32}(0)) == 0.0
        end

        j33 = besselj(3,3.)
        @test besselj(3,3) == j33
        @test besselj(-3,-3) == j33
        @test besselj(-3,3) == -j33
        @test besselj(3,-3) == -j33
        @test besselj(3,3f0) ≈ j33
        @test besselj(3,complex(3.)) ≈ j33
        @test besselj(3,complex(3f0)) ≈ j33
        @test besselj(3,complex(3)) ≈ j33

        j43 = besselj(4,3.)
        @test besselj(4,3) == j43
        @test besselj(-4,-3) == j43
        @test besselj(-4,3) == j43
        @test besselj(4,-3) == j43
        @test besselj(4,3f0) ≈ j43
        @test besselj(4,complex(3.)) ≈ j43
        @test besselj(4,complex(3f0)) ≈ j43
        @test besselj(4,complex(3)) ≈ j43

        @test j33 ≈ 0.30906272225525164362
        @test j43 ≈ 0.13203418392461221033
        @test besselj(0.1, complex(-0.4)) ≈ 0.820421842809028916 + 0.266571215948350899im
        @test besselj(3.2, 1.3+0.6im) ≈ 0.01135309305831220201 + 0.03927719044393515275im
        @test besselj(1, 3im) ≈ 3.953370217402609396im
        @test besselj(1.0,3im) ≈ besselj(1,3im)

        true_jm3p1_3 = -0.45024252862270713882
        @test besselj(-3.1,3) ≈ true_jm3p1_3
        @test besselj(Float16(-3.1),Float16(3)) ≈ true_jm3p1_3

        @testset "Error throwing" begin
            @test_throws DomainError    besselj(0.1, -0.4)
            @test_throws AmosException besselj(20,1000im)
            @test_throws MethodError besselj(big(1.0),3im)
        end
    end

    @testset "besselk" begin
        true_k33 = 0.12217037575718356792
        @test besselk(3,3) ≈ true_k33
        @test besselk(Float32(3),Float32(3)) ≈ true_k33
        @test besselk(Float16(3),Float16(3)) ≈ true_k33
        @test besselk(-3,3) ≈ true_k33
        true_k3m3 = -0.1221703757571835679 - 3.0151549516807985776im
        @test besselk(3,complex(-3)) ≈ true_k3m3
        @test besselk(-3,complex(-3)) ≈ true_k3m3
        # issue #6564
        @test besselk(1.0,0.0) == Inf
        @testset "Error throwing" begin
            @test_throws AmosException besselk(200,0.01)
            @test_throws DomainError besselk(3,-3)
            @test_throws MethodError besselk(1,big(1.0))
            @test_throws MethodError besselk(1,complex(big(1.0)))
            @test_throws MethodError besselkx(1,big(1.0))
            @test_throws MethodError besselkx(1,complex(big(1.0)))
        end
    end

    @testset "bessely" begin
        y33 = bessely(3,3.)
        @test bessely(3,3) == y33
        @test bessely(3.,3.) == y33
        @test bessely(3,Float32(3.)) ≈ y33
        @test bessely(-3,3) ≈ -y33
        @test y33 ≈ -0.53854161610503161800
        @test bessely(3,complex(-3)) ≈ 0.53854161610503161800 - 0.61812544451050328724im
        @testset "Error throwing" begin
            @test_throws AmosException bessely(200.5,0.1)
            @test_throws DomainError bessely(3,-3)
            @test_throws DomainError bessely(0.4,-1.0)
            @test_throws DomainError bessely(0.4,Float32(-1.0))
            @test_throws DomainError bessely(1,Float32(-1.0))
            @test_throws DomainError bessely(0.4,BigFloat(-1.0))
            @test_throws DomainError bessely(1,BigFloat(-1.0))
            @test_throws DomainError bessely(Cint(3),Float32(-3.))
            @test_throws DomainError bessely(Cint(3),Float64(-3.))

            @test_throws MethodError bessely(1.2,big(1.0))
            @test_throws MethodError bessely(1,complex(big(1.0)))
            @test_throws MethodError besselyx(1,big(1.0))
            @test_throws MethodError besselyx(1,complex(big(1.0)))
        end
    end

    @testset "sphericalbesselj" begin
        @test sphericalbesselj(1, 1)      ≈ 0.3011686789397568
        @test sphericalbesselj(10, 5.5)   ≈ 0.0009369210263385842
        @test sphericalbesselj(1.25, 5.5) ≈ -0.1123558799930763
        @test sphericalbesselj(1.25, -5.5+0im) ≈ 0.079447604649286 + 0.079447604649286im

        @test sphericalbesselj(0, 0.01) ≈ 0.999983333416666
        @test sphericalbesselj(0, 0)    == 1.0
        @test sphericalbesselj(1, 0)    == 0.0
        @test sphericalbesselj(1, 0.01) ≈ 0.003333300000119047
        
        @test_throws DomainError sphericalbesselj(1.25, -5.5)
    end

    @testset "sphericalbessely" begin
        @test sphericalbessely(1, 1)      ≈ -1.381773290676036
        @test sphericalbessely(10, 5.5)   ≈ -10.89087037026398
        @test sphericalbessely(1.25, 5.5) ≈ 0.148322390312228
        @test sphericalbessely(1.25, -5.5+0im) ≈ -0.054015441306998 - 0.104879767991574im

        @test sphericalbessely(0, 1e-5) ≈ -99999.9999950000000
        @test sphericalbessely(1, 1e-5) ≈ -1e10
        
        @test_throws DomainError sphericalbessely(1.25, -5.5)
        @test_throws AmosException sphericalbessely(1, 0)
    end

    @testset "besselhx" begin
        for elty in [Complex{Float16},Complex{Float32},Complex{Float64}]
            z = convert(elty, 1.0 + 1.9im)
            @test besselhx(1.0, 1, z) ≈ convert(elty,-0.5949634147786144 - 0.18451272807835967im)
            @test besselhx(Float32(1.0), 1, z) ≈ convert(elty,-0.5949634147786144 - 0.18451272807835967im)
        end
        @testset "Error throwing" begin
            @test_throws MethodError besselh(1,1,big(1.0))
            @test_throws MethodError besselh(1,1,complex(big(1.0)))
            @test_throws MethodError besselhx(1,1,big(1.0))
            @test_throws MethodError besselhx(1,1,complex(big(1.0)))
        end
    end
    @testset "scaled bessel[ijky] and hankelh[12]" begin
        for x in (1.0, 0.0, -1.0), y in (1.0, 0.0, -1.0), nu in (1.0, 0.0, -1.0)
            z = Complex{Float64}(x + y * im)
            z == zero(z) || @test hankelh1x(nu, z) ≈ hankelh1(nu, z) * exp(-z * im)
            z == zero(z) || @test hankelh2x(nu, z) ≈ hankelh2(nu, z) * exp(z * im)
            (nu < 0 && z == zero(z)) || @test besselix(nu, z) ≈ besseli(nu, z) * exp(-abs(real(z)))
            (nu < 0 && z == zero(z)) || @test besseljx(nu, z) ≈ besselj(nu, z) * exp(-abs(imag(z)))
            z == zero(z) || @test besselkx(nu, z) ≈ besselk(nu, z) * exp(z)
            z == zero(z) || @test besselyx(nu, z) ≈ bessely(nu, z) * exp(-abs(imag(z)))
        end
        @test besselkx(1, 0) == Inf
        for i = [-5 -3 -1 1 3 5]
            @test besseljx(i,0) == 0
            @test besselix(i,0) == 0
            @test besseljx(i,Float32(0)) == 0
            @test besselix(i,Float32(0)) == 0
            @test besseljx(i,Complex{Float32}(0)) == 0
            @test besselix(i,Complex{Float32}(0)) == 0
            @test besseljx(i,Float16(0)) == 0
            @test besselix(i,Float16(0)) == 0
            @test besseljx(i,Complex{Float16}(0)) == 0
            @test besselix(i,Complex{Float16}(0)) == 0
        end
        @testset "Error throwing" begin
            @test_throws AmosException hankelh1x(1, 0)
            @test_throws AmosException hankelh2x(1, 0)
            @test_throws AmosException besselix(-1.01, 0)
            @test_throws AmosException besseljx(-1.01, 0)
            @test_throws AmosException besselyx(1, 0)
            @test_throws DomainError besselix(0.4,-1.0)
            @test_throws DomainError besseljx(0.4, -1.0)
            @test_throws DomainError besselkx(0.4,-1.0)
            @test_throws DomainError besselyx(0.4,-1.0)
        end
    end
    @testset "issue #6653" begin
        @testset "$f" for f in (besselj,bessely,besseli,besselk,hankelh1,hankelh2)
            @test f(0,1) ≈ f(0,Complex{Float64}(1))
            @test f(0,1) ≈ f(0,Complex{Float32}(1))
            @test f(0,1) ≈ f(0,Complex{Float16}(1))
        end
    end

    @testset "jinc function" begin
        
        x1 = big"1.8"
        res1 = big"-0.116401402395234261196060221416100814277048989838340270527643894539781608507904800871249915340068849141209639297793338757544621191431028495364540423087719279432811811"
        x2 = big"0.00001"
        res2 = big"0.9999999998766299449914564074224855053704639002961349883799411022244488435839201448375732357069707105107902868318079100001712780899305533775708875568993566765287247356"
        x3 = big"0.02"
        res3 = big"0.9995066009475120760121935544369072906921634576001177613216005027724298234311657036029970452969291362640681823594159438475599412434802474217121139032145806264666152465837468052082455102325311416327533"
        
        root1 = big"1.21966989126650445492653884746525517787935933077511212945638126557694328028076014425087191879391333148570410142424433731699"
        for T in [Float16, Float32, Float64]
            @test jinc(T(x1)) ≅ T(res1)
            @test jinc(T(0)) ≅ T(1)
            x = root1 
            # check the first root for negative and positive
            @test jinc(T(x)) + T(1) ≅ T(1)
            @test jinc(T(-x)) + T(1) ≅ T(1)

            # check for symmetry
            for x = -5.0:1.0:3.0
                @test jinc(T(10^x)) ≅ jinc(T(- 10^x))

            end

            x = T(x2)
            res = T(res2)
            @test jinc(x) ≅ res

            x = T(x3)
            res = T(res3)
            @test jinc(x) ≅ res

            # type stability for 0
            @test typeof(jinc(T(0))) === T
        end

        # test big float
        begin
            T = BigFloat
            @test T(jinc(T(x1))) ≈  T(res1)
            #= @test jinc(T(0)) ≈ T(1) =#
            x = root1
            # check the first root for negative and positive
            @test jinc(T(x)) + T(1) ≈ T(1)
            @test jinc(T(-x)) + T(1) ≈ T(1)

            # check for symmetry
            for x = -5.0:1.0:3.0
                @test jinc(T(10^x)) ≈ jinc(T(- 10^x))

            end

            x = T(x2)
            res = T(res2)
            @test T(jinc(x)) ≈ res

            x = T(x3)
            res = T(res3)
            @test T(jinc(x)) ≈ res

            # type stability for 0
            @test typeof(jinc(T(0))) === T
        end


        x1 = big"1.0" + big"0.0001"*1im
        res1 = big"0.1811917493658519438199345187518953631147425189362440669941643986059307908067224704314102961998438798566229144099942611910610729246007355473812330821120666591291397178705295331963573655500130149361351" - big"0.0000970867870757251899047609572999879575662488432474075980361828257989869563037501527523384396635449088435640500216646530541742421182815911923875286404692082399771721882711636718746219759416544606077"*1im

        x2 = big"0.00001" + big"10"*1im
        res2 = big"1.971105539196620210542904558564976891448851717386016064466205835417372529649732623778146798009632930245313240967710133352055225553347549825667896762467087506209685330050397030027682762563950964962993e11" - big"5.89917589196246805438832769814196828580119495400655210987269285701113079909408122244782544425726443894463237062342568390611353700647229983489787030738912900107228049806443988000825445873073676204e6"*1im

        x3 = big"0.00002"*1im
        res3 = big"1.0000000004934802201356421734767362282021102443681623763651414603818274528461323973841772525693725599437763705213335500965076723534424132751557185487068200219257886657409530057472519875335446749379510"

        # check also some complex numbers
        # besselj1 is not type stable and returns Complex{Float64}
        # the stricter ≅ doesn't hold here
        # besselj1 seems to produce inaccurate results for some Float32 cases
        for T in [Complex{Float16}, Complex{Float32}, Complex{Float64}]
            x = T(x1)
            res = T(res1)
            @test T(jinc(x)) ≈ res

            x = T(x2)
            res = T(res2)
            @test T(jinc(x)) ≈ res

            x = T(x3)
            res = T(res3)
            @test T(jinc(x)) + T(1im) ≅ res + T(1im) 

        end

    end
end
