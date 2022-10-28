# From openlibm/test/libm-test-ulps.h, openlibm/test/libm-test.c

using SpecialFunctions: _lgamma_r, _loggamma_r, _lgammaf_r, _loggammaf_r

# lgamma_test block
for (T, lgamma, labsgamma) in ((Float64, _loggamma_r, _lgamma_r),
                               (Float32, _loggammaf_r, _lgammaf_r))
    @testset "lgamma_test, $T" begin
        @test lgamma(T(Inf)) === T(Inf)
        @test lgamma(T(0)) === T(Inf)
        @test lgamma(T(NaN)) === T(NaN)

        @test lgamma(T(-3)) === T(Inf)
        @test lgamma(T(-Inf)) === T(Inf)

        # lgamma(1) == 0, lgamma (1) sets signgam to 1
        y, signgam = labsgamma(T(1))
        @test y === T(0.0)
        @test signgam == 1

        # lgamma(3) == log(2), lgamma (3) sets signgam to 1
        y, signgam = labsgamma(T(3))
        @test y === log(T(2.0))
        @test signgam == 1

        # lgamma(0.5) == log(sqrt(pi)), lgamma(0.5) sets signgam to 1
        y, signgam = labsgamma(T(0.5))
        @test y === T(0.5log(π))
        @test signgam == 1

        # lgamma(-0.5) == log(2sqrt(pi)), lgamma(-0.5) sets signgam to -1
        y, signgam = labsgamma(T(-0.5))
        @test y === T(0.5log(4π))
        @test signgam == -1
        @test_throws DomainError lgamma(T(-0.5))

        # In the two "broken" tests, an exact match not possible, even
        # in Float64, thus, we check for as close a tolerance as
        # possible.

        # lgamma(0.7) == 0.26086724653166651439, lgamma(0.7) sets signgam to 1
        y, signgam = labsgamma(T(0.7))
        # @test_broken y === 0.26086724653166651439
        if T === Float64
            @test y ≈ 0.26086724653166651439 atol=6e-17
        else
            @test y ≈ 0.26086724653166651439 atol=3e-8
        end
        @test signgam == 1

        # lgamma(1.2) == -0.853740900033158497197e-1, lgamma(1.2) sets signgam to 1
        y, signgam = labsgamma(T(1.2))
        # @test_broken y === -0.853740900033158497197e-1
        if T === Float64
            @test y ≈ -0.853740900033158497197e-1 atol=2e-17
        else
            @test y ≈ -0.853740900033158497197e-1 atol=2e-8
        end
        @test signgam == 1
    end
end
