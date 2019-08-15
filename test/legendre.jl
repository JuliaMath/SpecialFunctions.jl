using Polynomials: Poly

@testset "legendre and related functions" begin
    @testset "legendre" begin
        n_poly  = 6
        c       = zeros(n_poly, n_poly)
        c[1,1]  = 1
        c[2,2]  = 1
        c[3, 1:3] .= [-1,  0,   3               ] / 2
        c[4, 1:4] .= [ 0, -3,   0,   5          ] / 2
        c[5, 1:5] .= [ 3,  0, -30,   0, 35      ] / 8
        c[6, 1:6] .= [ 0, 15,   0, -70,  0, 63  ] / 8

        n_x = 20
        x_arr = range(-1, 1, length=n_x)
        for n = 0:n_poly-1
            P = Poly(c[n+1,:])
            for x in x_arr
                @test legendreP(n, x) ≈ P(x)        rtol=1e-14
            end
        end

        @test_throws DomainError legendreP(-1, 0)
        @test_throws DomainError legendreP( 1, 2)
        @test_throws DomainError legendreP(-1, 2)

        @test_throws MethodError legendreP(0, Complex(1))
    end
end
