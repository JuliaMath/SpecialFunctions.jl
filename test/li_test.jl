using Base.Test

include("li.jl")
ϕ = (sqrt(5)+1)/2.0

@testset "Check different input formats" begin
    @test Li(-1.0, 0.0) ≈ 0.0
    @test Li(-1, 0.0) ≈ 0.0
    @test Li(-1.0, Complex(0.0)) ≈ 0.0 + 0.0im # don't seem to need imaginary bit here
    @test Li(Complex(-1.0), 0.0) ≈ 0.0 + 0.0im
    @test Li(Complex(-1.0), 0.0) ≈ 0.0
    @test Li(1, 0) ≈ 0.0
    x = collect(0.0:0.1:0.9)
    @test all([Li.(1, x)[i] ≈ Li(1, x[i]) for i=1:length(x)])
    @test Li(Complex(-1.0), Complex(0.3)) ≈ Li(-1.0, 0.3)
end

@testset "Test s = n (a real integer)" begin
    # simple case
    @test Li(1, 0.5) ≈ log(2)

    @testset "   dilogarithm for real z" begin
        @test Li(2,-1.0) ≈ -pi^2/12.0
        @test Li(2, 0.0) ≈ 0.0
        @test Li(2, 0.5) ≈ pi^2/12 - 0.5*log(2)^2
        @test Li(2, 1.0) ≈ pi^2/6.0
        @test Li(2, 2.0) ≈ pi^2/4 - im*pi*log(2)
        @test Li(2, -ϕ)     ≈ -pi^2/10 - log(ϕ)^2
        @test Li(2, -1/ϕ)   ≈ -pi^2/15 + log(ϕ)^2/2
        @test Li(2, 1/ϕ^2)  ≈  pi^2/15 - log(ϕ)^2
        @test Li(2, 1/ϕ)    ≈  pi^2/10 - log(ϕ)^2
        # wiki has this one, but no ref:
        @test Li(2, ϕ)      ≈  11*pi^2/15 + log(Complex(-1/ϕ))^2/2
        @test Li(2, ϕ^2)    ≈ -11*pi^2/15 - log(Complex(-ϕ))^2

        # identities
        Z = [3.0 + 0.4im, -3.0 + 0.4im, 3.0 - 0.4im, -3.0 + -0.4im]
        for i=1:length(Z)
            z = Z[i]
            @test Li(2, z) + Li(2, 1/z) ≈ -pi^2/6.0 - log(Complex(-z))^2/2.0
        end
    end

    @testset "   trilogarithm for real z" begin
        @test Li(3,-1.0) ≈ -3*zeta(3)/4
        @test Li(3, 0.0) ≈ 0.0
        @test Li(3, 0.5) ≈ log(2)^3/6.0 - pi^2*log(2)/12.0 + (7.0/8.0)*zeta(3)
        @test Li(3, 1.0) ≈ zeta(3)
        @test Li(3, ϕ^(-2)) ≈ 4*zeta(3)/5 + 2*log(ϕ)^3/3 - 2*pi^2*log(ϕ)/15
        
        # identities
        z = 3.0 + 0.4im # test for random complex z
        z = 1.5 # test for random complex z
        # @test Li(3, z) + Li(3, 1-z) + Li(3, 1 - 1/z) ≈ zeta(3) + log(z)^3/6  + pi^2*log(z)/6 - 0.5*log(z)^2*log(complex(1-z))
        # @test Li(3, z) - Li(3, -z) ≈ 0.25 * Li(3, z^2)
    end
        
    @testset "   general case for real z" begin
        # X = collect(-0.95:0.1:0.95)
        X = collect(-3.0:0.1:3.0)
        for i=1:length(X)
            x = X[i]
            # println(x)
            @test Li(1, x) ≈ -log(Complex(1-x))
            @test Li(0, x) ≈ x ./ (1-x)
            @test Li(-1, x) ≈ x ./ (1-x).^2
            @test Li(-2, x) ≈ x .* (1+x) ./ (1-x).^3
            @test Li(-3, x) ≈ x .* (1+4*x+x.^2) ./ (1-x).^4
            @test Li(-4, x) ≈ x .* (1+x) .* (1+10*x+x.^2) ./ (1-x).^5
        end
    end
   
    @testset "   general case for complex z" begin
        X = collect(-3.0:0.5:3.0)
        Y = [-1.3, -0.4, 0.4, 1.5]
        for i=1:length(X)
            for j=1:length(Y)
                z = Complex(X[i], Y[j])
                # println(z)
                @test Li(1, z) ≈ -log(Complex(1-z))
                @test Li(0, z) ≈ z ./ (1-z)
                @test Li(-1, z) ≈ z ./ (1-z).^2
                @test Li(-2, z) ≈ z .* (1+z) ./ (1-z).^3
                @test Li(-3, z) ≈ z .* (1+4*z+z.^2) ./ (1-z).^4
                @test Li(-4, z) ≈ z .* (1+z) .* (1+10*z+z.^2) ./ (1-z).^5
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
            @test Li(s, -1.0) ≈ -eta(s)
            @test Li(s,  1.0) ≈ zeta(s)
            @test Li(s,  im) ≈  - 2.0.^(-s).*eta.(s) + im*Dbeta.(s) 
        end
    end
end

@testset "Additional Identities" begin
    z = 0.5
    for n=1:5
        @test Li(-n,z) + (-1)^n * Li(-n, 1/z) ≈ 0.0 
    end

    # for real s, and real z<1, Li should be real
    S = [-1, 0.1, 2]
    Z = [-2, -1.0, 0.1, 0.95]
    for i=1:length(S)
        for j=1:length(Z)
            s = S[i]
            z = Z[j]
            # println("s = $s; z = $z")
            @test abs( imag( Li(s,  z) ) ) < 1.0e-14
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
            @test imag( Li(s,  z) ) ≈ -pi*μ^(s-1)/gamma(s)
        end
    end

    S_r = [2.1 2.5 3.0]
    S_i = [-1.3, -1.0, -0.5, 0.0, 0.5, 1.0, 1.3]
    z = 3.0 - 0.1im
    for i=1:length(S_r)
        for j=1:length(S_i)
            s = Complex(S_r[i], S_i[j])
            @test Li(s,z) + Li(s,-z) ≈ complex(2)^(1-s) * Li(s, z^2)
        end
    end
end
