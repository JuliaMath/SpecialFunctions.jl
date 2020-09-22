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
    @testset "bessel{j,y}{0,1}: check return value wrt bessel{j,y}" begin
        for jy in ("j","y")
            bjy = Symbol("bessel",jy)
            for F in [Float16, Float32]
                for nu in (0, 1)
                    bjynu = Symbol("bessel",jy,nu)

                    @test $bjynu(         F(2.0)   ) ≈ $bjynu(Float64(2.0)  )
                    @test $bjynu(           2      ) ≈ $bjynu(        2.0   )
                    @test $bjynu(           2.0    ) ≈ $bjy(nu,       2.0   )
                    @test $bjynu(           2.0+im ) ≈ $bjy(nu,       2.0+im)
                    @test $bjynu(Complex{F}(2.0+im)) ≈ $bjynu(        2.0+im)
                end
            end
        end
    end

    @testset "besselj, bessely: correct return type" begin
        @testset "type stability: $f" for f in [bessely0, bessely1, besselj0, besselj1]
            for F in [Float16, Float32, Float64]
                @test         F  == Base.return_types(f, Tuple{        F })[]
                @test Complex{F} == Base.return_types(f, Tuple{Complex{F}})[]
            end
            @test BigFloat == Base.return_types(f, Tuple{BigFloat})[]
        end
        @testset "type stability: $f" for f in [besselj, bessely]
            for F in [Float16, Float32, Float64]
                @test F == Base.return_types(f, Tuple{Int, F})[]
            end
        end
    end

    @testset "besselj: undefined argument types" begin
        @test_throws MethodError   besselj(     1.2 ,         big( 1.0        ))
        @test_throws MethodError   besselj(     1   , complex(big( 1.0       )))
        @test_throws MethodError   besselj(big( 1.0),                     3im  )
        @test_throws DomainError   besselj(     0.1 ,             -0.4         )
        @test_throws AmosException besselj(    20   ,                  1000im  )

        @test_throws MethodError besseljx(1,  big(        1.0) )
        @test_throws MethodError besseljx(1,  complex(big(1.0)))

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
    @testset "besselj: specific values and (domain) errors" begin
        # besselj(nu, 0) == { 1     for nu == 0
        #                   { 0     else
        for nu in [-5 -3 -1 0 1 3 5]
            for I in [Int16, Int32, Int64]
                @test besselj(nu, zero(        I )) == (nu == 0 ? 1 : 1)
            end
            for F in [Float16, Float32, Float64]
                @test besselj(nu, zero(        F )) == (nu == 0 ? 1 : 1)
                @test besselj(nu, zero(Complex{F})) == (nu == 0 ? 1 : 1)
            end
        end

        const j33 = besselj(3, 3.0)
        @test besselj( 3,  3) ==  j33
        @test besselj(-3, -3) ==  j33
        @test besselj(-3,  3) == -j33
        @test besselj( 3, -3) == -j33
        for F in [Float16, Float32, Float64]
            @test besselj(3,         F( 3)) ≈ j33
            @test besselj(3, Complex{F}(3)) ≈ j33
        end

        const j43 = besselj(4, 3.0)
        @test besselj( 4,  3) ==  j43
        @test besselj(-4, -3) ==  j43
        @test besselj(-4,  3) == -j43
        @test besselj( 4, -3) == -j43
        for F in [Float16, Float32, Float64]
            @test besselj(4,         F( 3)) ≈ j43
            @test besselj(4, Complex{F}(3)) ≈ j43
        end

        @test besselj(1.0, 3im) ≈ besselj(1, 3im)

        # specific values
        @test j33                      ≈ 0.30906272225525164362
        @test j43                      ≈ 0.13203418392461221033
        @test besselj(0.1, -0.4+  0im) ≈ 0.820421842809028916   + 0.266571215948350899im
        @test besselj(3.2,  1.3+0.6im) ≈ 0.01135309305831220201 + 0.03927719044393515275im
        @test besselj(1,          3im) ≈                          3.953370217402609396im
        for F in [Float16, Float32, Float64]
            @test besselj(F(-3.1), F(3)) ≈ -0.45024252862270713882
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

    @testset "bessely: specific values and (domain) errors" begin
        const y33 = bessely(3, 3.0)
        @test bessely(3, 3) == y33
        @test bessely(3.0,3.0) == y33
        @test bessely(3,Float32(3.0)) ≈ y33
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
            @test_throws DomainError bessely(Cint(3),Float32(-3.0))
            @test_throws DomainError bessely(Cint(3),Float64(-3.0))

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
end
