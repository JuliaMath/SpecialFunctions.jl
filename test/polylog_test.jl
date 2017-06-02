using SpecialFunctions
using Base.Test

const SF = SpecialFunctions

# test functions, see runtests.jl
relerr(z, x) = z == x ? 0.0 : abs(z - x) / abs(x)
relerrc(z, x) = max(relerr(real(z),real(x)), relerr(imag(z),imag(x)))
≅(a,b) = relerrc(a,b) ≤ 1e-13

# test functions, similar to runtests.jl, but can't do relative error when some values are zero
abserr(z, x) = z == x ? 0.0 : abs(z - x)
abserrc(z, x) = max(abserr(real(z),real(x)), abserr(imag(z),imag(x)))
≒(a,b) = abserrc(a,b) ≤ 1e-13

# golden ratio is used in some tests
ϕ = (sqrt(5)+1)/2.0

@testset "Polylogarithm polylog" begin
    @testset "Check different input formats" begin
        @test SF.polylog(-1.0, 0.0) ≈ 0.0
        @test SF.polylog(-1, 0.0) ≈ 0.0
        @test SF.polylog(-1.0, Complex(0.0)) ≈ 0.0 + 0.0im # don't seem to need imaginary bit here
        @test SF.polylog(Complex(-1.0), 0.0) ≈ 0.0 + 0.0im
        @test SF.polylog(Complex(-1.0), 0.0) ≈ 0.0
        @test SF.polylog(1, 0) ≈ 0.0
        x = collect(0.0:0.1:0.9)
        @test all([SF.polylog.(1, x)[i] ≈ SF.polylog(1, x[i]) for i=1:length(x)])
        @test SF.polylog(Complex(-1.0), Complex(0.3)) ≈ SF.polylog(-1.0, 0.3)
        @test SF.polylog(Complex(-1.0), Complex(0.3)) ≈ SF.polylog(-1.0, 0.3)
    end

    @testset "Test s = n (a real integer)" begin
        # simple cases
        @test SF.polylog(1, 0.5) ≈ log(2)
        @test !isfinite( SF.polylog(1, 1.0) )
        @test SF.polylog(1, 2) ≒ -pi*im # comes from -log(complex(1-2))
    
        @testset "   dilogarithm for real z" begin
            @test SF.polylog(2,-1.0) ≈ -pi^2/12.0
            @test SF.polylog(2, 0.0) ≈ 0.0
            @test SF.polylog(2, 0.5) ≈ pi^2/12 - 0.5*log(2)^2
            @test SF.polylog(2, 1.0) ≈ pi^2/6.0
            @test SF.polylog(2, 2.0) ≈ pi^2/4 - im*pi*log(2)
            @test SF.polylog(2, -ϕ)     ≈ -pi^2/10 - log(ϕ)^2
            @test SF.polylog(2, -1/ϕ)   ≈ -pi^2/15 + log(ϕ)^2/2
            @test SF.polylog(2, 1/ϕ^2)  ≈  pi^2/15 - log(ϕ)^2
            @test SF.polylog(2, 1/ϕ)    ≈  pi^2/10 - log(ϕ)^2
            @test SF.polylog(2, ϕ)      ≈  11*pi^2/15 + log(Complex(-1/ϕ))^2/2 # wiki has this one, but no ref
            @test SF.polylog(2, ϕ^2)    ≈ -11*pi^2/15 - log(Complex(-ϕ))^2

            # identities
            Z = [3.0 + 0.4im, -3.0 + 0.4im, 3.0 - 0.4im, -3.0 + -0.4im]
            for i=1:length(Z)
                z = Z[i]
                @test SF.polylog(2, z) + SF.polylog(2, 1/z) ≈ -pi^2/6.0 - log(Complex(-z))^2/2.0
            end
        end

        @testset "   trilogarithm for real z" begin
            @test SF.polylog(3,-1.0) ≈ -3*SF.zeta(3)/4
            @test SF.polylog(3, 0.0) ≈ 0.0
            @test SF.polylog(3, 0.5) ≈ log(2)^3/6.0 - pi^2*log(2)/12.0 + (7.0/8.0)*SF.zeta(3)
            @test SF.polylog(3, 1.0) ≈ SF.zeta(3)
            @test SF.polylog(3, ϕ^(-2)) ≈ 4*SF.zeta(3)/5 + 2*log(ϕ)^3/3 - 2*pi^2*log(ϕ)/15
        end
        
        @testset "   general case for real z" begin
            # X = collect(-0.95:0.1:0.95)
            X = collect(-3.0:0.1:3.0)
            for i=1:length(X)
                x = X[i]
                # println(x)
                @test SF.polylog(1, x) ≈ -log(Complex(1-x))
                @test SF.polylog(0, x) ≈ x ./ (1-x)
                @test SF.polylog(-1, x) ≈ x ./ (1-x).^2
                @test SF.polylog(-2, x) ≈ x .* (1+x) ./ (1-x).^3
                @test SF.polylog(-3, x) ≈ x .* (1+4*x+x.^2) ./ (1-x).^4
                @test SF.polylog(-4, x) ≈ x .* (1+x) .* (1+10*x+x.^2) ./ (1-x).^5
            end
        end
        
        @testset "   general case for complex z" begin
            X = collect(-3.0:0.5:3.0)
            Y = [-1.3, -0.4, 0.4, 1.5]
            for i=1:length(X)
                for j=1:length(Y)
                    z = Complex(X[i], Y[j])
                    # println(z)
                    @test SF.polylog(1, z) ≈ -log(Complex(1-z))
                    @test SF.polylog(0, z) ≈ z ./ (1-z)
                    @test SF.polylog(-1, z) ≈ z ./ (1-z).^2
                    @test SF.polylog(-2, z) ≈ z .* (1+z) ./ (1-z).^3
                    @test SF.polylog(-3, z) ≈ z .* (1+4*z+z.^2) ./ (1-z).^4
                    @test SF.polylog(-4, z) ≈ z .* (1+z) .* (1+10*z+z.^2) ./ (1-z).^5
                end
            end
        end    
    end

    @testset "Particular values |z| == 1" begin
        S_r = [2.1 2.5 3.0]
        S_i = [-1.3, -1.0, -0.5, 0.0, 0.5, 1.0, 1.3]
        for i=1:length(S_r)
            for j=1:length(S_i)
                s = Complex(S_r[i], S_i[j])
                @test SF.polylog(s, -1.0) ≈ -SF.eta(s)
                @test SF.polylog(s,  1.0) ≈  SF.zeta(s)
                @test SF.polylog(s,  im) ≈  - 2.0.^(-s).*SF.eta.(s) + im*SF.Dbeta.(s) 
            end
        end

        # the singularity
        @test !isfinite( SF.polylog(0.5, 1.0) )
    end

    @testset "Additional Identities" begin
        z = 0.5
        for n=1:5
            @test SF.polylog(-n,z) + (-1)^n * SF.polylog(-n, 1/z) ≈ 0.0 
        end

        # for real s, and real z<1, SF.polylog should be real
        S = [-1, 0.1, 2]
        Z = [-2, -1.0, 0.1, 0.95]
        for i=1:length(S)
            for j=1:length(Z)
                s = S[i]
                z = Z[j]
                # println("s = $s; z = $z")
                @test abs( imag( SF.polylog(s,  z) ) ) < 1.0e-14
            end
        end

        # for real s, and real z>=1, the imaginary part is give
        S = [-1.5, 0.1, 2]
        Z = [1.05, 3.0]
        for i=1:length(S)
            for j=1:length(Z)
                s = S[i]
                z = Z[j]
                # println("s = $s; z = $z")
                μ = log(z)
                @test imag( SF.polylog(s,  z) ) ≈ -pi*μ^(s-1)/gamma(s)
            end
        end

        S_r = [2.1 2.5 3.0]
        S_i = [-1.3, -1.0, -0.5, 0.0, 0.5, 1.0, 1.3]
        z = 3.0 - 0.1im
        for i=1:length(S_r)
            for j=1:length(S_i)
                s = Complex(S_r[i], S_i[j])
                @test SF.polylog(s,z) + SF.polylog(s,-z) ≈ complex(2)^(1-s) * SF.polylog(s, z^2)
            end
        end
    end
end
