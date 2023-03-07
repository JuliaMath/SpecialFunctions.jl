@testset "vectorization of 2-arg functions" begin
    binary_math_functions = [
        besselh, hankelh1, hankelh2, hankelh1x, hankelh2x,
        besseli, besselix, besselj, besseljx, besselk, besselkx, bessely, besselyx,
        polygamma, zeta
    ]
    @testset "$f" for f in binary_math_functions
        x = y = 2
        v = [f(x,y)]
        @test f.([x],y) == v
        @test f.(x,[y]) == v
        @test f.([x],[y]) == v
    end
end

@testset "MPFR" begin
    @testset "bessel functions" begin
        setprecision(53) do
            @test besselj(4, BigFloat(2)) ≈ besselj(4, 2.)
            @test besselj0(BigFloat(2)) ≈ besselj0(2.)
            @test besselj1(BigFloat(2)) ≈ besselj1(2.)
            @test bessely(4, BigFloat(2)) ≈ bessely(4, 2.)
            @test bessely0(BigFloat(2)) ≈ bessely0(2.)
            @test bessely1(BigFloat(2)) ≈ bessely1(2.)
        end
    end

    let err(z, x) = abs(z - x) / abs(x)
        @test 1e-60 > err(eta(parse(BigFloat,"1.005")), parse(BigFloat,"0.693945708117842473436705502427198307157819636785324430166786"))
        @test 1e-60 > err(exp(eta(big(1.0))), 2.0)
    end

    let a = parse(BigInt, "315135")
        @test typeof(erf(a)) == BigFloat
        @test typeof(erfc(a)) == BigFloat
    end

    # issue #101
    for i in 0:5
        @test gamma(big(i)) == gamma(i)
    end
end

@testset "Base Julia issue #17474" begin
    @test f64(complex(1f0,1f0)) === complex(1.0, 1.0)
    @test f64(1f0) === 1.0

    @test typeof(eta(big"2")) == BigFloat
    @test typeof(zeta(big"2")) == BigFloat
    @test typeof(digamma(big"2")) == BigFloat

    @test_throws MethodError trigamma(big"2")
    @test_throws MethodError trigamma(big"2.0")
    @test_throws MethodError invdigamma(big"2")
    @test_throws MethodError invdigamma(big"2.0")

    @test_throws MethodError eta(Complex(big"2"))
    @test_throws MethodError eta(Complex(big"2.0"))
    @test_throws MethodError zeta(Complex(big"2"))
    @test_throws MethodError zeta(Complex(big"2.0"))
    @test_throws MethodError zeta(1.0,big"2")
    @test_throws MethodError zeta(1.0,big"2.0")
    @test_throws MethodError zeta(big"1.0",2.0)
    @test_throws MethodError zeta(big"1",2.0)

    @test polygamma(3, 0x2) isa Float64
    @test polygamma(big"3", 2f0) isa Float32
    @test zeta(1, 2.0) isa Float64
    @test zeta(1, 2f0) isa Float32
    @test zeta(2f0, complex(2f0,0f0)) isa Complex{Float32}
    @test zeta(complex(1,1), 2f0) isa Complex{Float32}
    @test zeta(complex(1), 2.0) isa Complex{Float64}
end

@test sprint(showerror, AmosException(1)) == "AmosException with id 1: input error."

@testset "missing data" begin
    for f in (digamma, erf, erfc, erfcinv, erfcx, erfi, erfinv, eta, gamma,
              invdigamma, logfactorial, trigamma)
        @test f(missing) === missing
    end
    @test beta(1.0, missing) === missing
    @test beta(missing, 1.0) === missing
    @test beta(missing, missing) === missing
    @test polygamma(4, missing) === missing
end

@testset "Test fastabs" begin
    numbers = [1, 2, 0, 1e3]
    for n in numbers
        @test abs(n) == SpecialFunctions.fastabs(n)
    end
            
    numbers = [1im, 2 + 2im, 0 + 100im, 1e3 + 1e-10im]
    for n in numbers
        @test abs(real(n)) + abs(imag(n)) == SpecialFunctions.fastabs(n)
    end
end
