using SpecialFunctions
using Base.Test

const SF = SpecialFunctions

@testset "Harmonic numbers" begin
    @testset "    basics" begin
        @test SF.harmonic(1) ≈ 1.0
        @test SF.harmonic(2) ≈ 1.5
        @test SF.harmonic(3) ≈ 11.0/6.0
        @test SF.harmonic(4) ≈ 25.0/12.0
        @test SF.harmonic(5) ≈ 137.0/60.0
    end

    @testset "    identities" begin
        for n=4:10
            @test SF.harmonic(n) ≈ 1.0/n + SF.harmonic(n-1)
        end
        
        for n=20:30
            @test SF.harmonic(n) ≈ sum( 1.0./(1:n) )
        end
    end


    @testset "Generalized harmonic numbers" begin
        for n=1:10
            @test SF.harmonic(n,1.0) ≈ SF.harmonic(n)
            @test SF.harmonic(n,0.0) ≈ n
            
            for r = 3:2:9
                @test SF.harmonic(n,r) ≈ Float64(n)^(-r) + SF.polygamma(r-1,n)/gamma(r) + SF.zeta(r)
            end
        end
    end
    
end
