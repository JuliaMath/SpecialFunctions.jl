@testset "elliptic" begin

setprecision(BigFloat,256)

function Base.read(s::IO, ::Type{BigFloat})
    setprecision(BigFloat, 256) do
        x = BigFloat()
        x.sign = Cint(read(s,Int8))
        unsafe_read(s, reinterpret(Ptr{UInt8},x.d), 32)
        x.exp = Clong(ntoh(read(s,Int64)))

        if x.sign == 0
            return big(0)
        elseif x.sign == 2
            return big(-Inf)
        elseif x.sign == 3
            return big(Inf)
        elseif x.sign == 4
            return big(NaN)
        else
            @assert unsafe_load(reinterpret(Ptr{UInt8},x.d), 32) >= 128
            # ^ MPFR crashes if the first bit in the mantissa is not set
            return x
        end
    end
end
Base.read(s::IO, ::Type{Complex{BigFloat}}) =
    complex(read(s,BigFloat),read(s,BigFloat))

@testset "io" begin
    open("elliptic/io.bin","r") do f
        @test ntoh(read(f,Int64)) ==  0
        @test ntoh(read(f,Int64)) ==  1
        @test ntoh(read(f,Int64)) ==  -1
        @test read(f,BigFloat) ==  0
        @test read(f,BigFloat) ==  1
        @test read(f,BigFloat) ==  2
        @test read(f,BigFloat) == -1
        @test read(f,Complex{BigFloat}) == im
        @test read(f,BigFloat) == sqrt(big(2))
        @test read(f,BigFloat) == Inf
        @test read(f,BigFloat) == -Inf
        @test isnan(read(f,BigFloat))
    end
end

@testset "ellipj" begin

    @testset "type stability" begin
        realtypes = (Int,Float32,Float64,BigFloat)
        types = (realtypes..., complex.(realtypes)...)
        @testset "typeof(u) = $U, typeof(m) = $M" for U in types, M in types
            @inferred ellipj(zero(U),zero(M))
        end
    end

    @testset "special values" begin
        vals = (0.0,Inf,NaN)
        @testset "u = $u, m = $m" for u in vals, m in vals
            if !isfinite(u) || !isfinite(m)
                @test all(isnan.(ellipj(u,m)))
            end
        end
    end

    @testset "mpmath" begin
        open("elliptic/ellipj.bin","r") do f
            npts = ntoh(read(f,Int64))
            for i = 1:npts
                u = read(f, Complex{BigFloat})
                m = read(f, Complex{BigFloat})
                snref = read(f, Complex{BigFloat})
                cnref = read(f, Complex{BigFloat})
                dnref = read(f, Complex{BigFloat})
                sn,cn,dn = ellipj(u,m)
                atol = sqrt(eps(BigFloat))*maximum(abs.((snref,cnref,dnref)))
                @test sn ≈ snref atol=atol
                @test cn ≈ cnref atol=atol
                @test dn ≈ dnref atol=atol
            end
        end
    end

    @testset "precision ($T)" for T in (Float16,Float32,Float64)
        n = 16
        s = @. exp(2T(π)*im*(0:n-1)/n)
        r = T[eps(T), π*eps(T), sqrt(eps(T)), π*sqrt(eps(T)), 1/π, 0.5, 3/π, 1]
        x = Complex{T}[0; vec(s.*r')]
        for u = x
            for m0 = (0,1)
                for m = m0 .+ x
                    sn,cn,dn = ellipj(u,m)
                    snref,cnref,dnref = ellipj(big(u),big(m))

                    rtol = 0
                    atol = 10*eps(T)*maximum(abs.((snref,cnref,dnref)))
                    @test sn ≈ snref   rtol=rtol atol = atol
                    @test cn ≈ cnref   rtol=rtol atol = atol
                    @test dn ≈ dnref   rtol=rtol atol = atol
                end
            end
        end
    end

end

@testset "ellipK" begin
    @testset "type stability" begin
        realtypes = (Int,Float32,Float64,BigFloat)
        types = (realtypes..., complex.(realtypes)...)
        @testset "typeof(m) = $M" for M in types
            @inferred ellipK(zero(M))
        end
    end

    @testset "special values" begin
        @test isnan(ellipK(NaN))
        @test ellipK(Inf+0im) == 0.0
        @test ellipK(-Inf) == 0.0
    end

    @testset "mpmath" begin
        open("elliptic/ellipK.bin","r") do f
            npts = ntoh(read(f,Int64))
            for i = 1:npts
                m = read(f, Complex{BigFloat})
                ellipKref = read(f, Complex{BigFloat})
                @test ellipK(m) ≈ ellipKref
            end
        end
    end

    @testset "precision ($T)" for T in (Float16,Float32,Float64)
        n = 16
        s = @. exp(2T(π)*im*(0:n-1)/n)
        r = T[eps(T), π*eps(T), sqrt(eps(T)), π*sqrt(eps(T)), 1/T(π), 0.5, 1, π, 1/sqrt(eps(T)), 1/eps(T)]
        x = Complex{T}[0; vec(s.*r')]
        for m0 = (0,1)
            for m = m0 .+ x
                @test ellipK(m) ≈ ellipK(big(m)) rtol=3*eps(T)
            end
        end
    end
end

end
