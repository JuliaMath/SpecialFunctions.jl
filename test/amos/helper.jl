
"Special Float32"
const SPECIAL_FLOAT32 = Float64[
    0.0,        # +0.0
    0x1p-149,   # Min Subnormal Float
    0x1p-127,   # mid Subnormal Float
    0x1.fffffcp-127,  # Max Subnormal Float
    0x1p-126,   # Min Normalized Float
    0x1p+1,     # 2
    0x1.fffffep+127,  # Max Normalized Float
    Inf,        # Inf
    reinterpret(Float32, 0x7f800001),  # NaN (+min)
    NaN,        # NaN (+min)
    reinterpret(Float32, 0x7fffffff),  # NaN (+max)

    # eps( 10^0 ~ 10^16 )
    [ eps(10.0^i) for i in 0:16 ]...,
] # SPECIAL_FLOAT32
# TOOD: add Float64 const


"""
input:
    zarr[ (+) ]

output:
    zarr[
        (+), (-),
        inv.(+), inv.(-),
    ]
"""
function gen_neg_inv(zarr::Vector{Float64})
    vcat(
        zarr, -1.0*zarr,
        inv.(zarr), -1.0*inv.(zarr),
    )
end

"""
input:
    zarr[ (+,+) ]

output:
    zarr[
        (-,+),  (+,+),
        (-,-),  (+,-),
    ]
"""
function gen_phase4(zarr::Vector{ComplexF64})
    vcat(
        1.0im * zarr,
        zarr
        -1.0 * zarr,
        -1.0im * zarr,
    )
end
function gen_phase4(zvv::Vector{Vector{ComplexF64}})
    gen_phase4.(zvv)
end

"""test if ComplexF64 contains NaN or Inf in real or imag parts.
"""
function contains_inf_nan(z::Vector{ComplexF64})
    z |> real .|> isinf |> any ||
        z |> imag .|> isinf |> any ||
        z |> real .|> isnan |> any ||
        z |> imag .|> isnan |> any
end


# ==== Test Set ====#

@testset "AMOS.uchk" begin
    function tuchk(y::ComplexF64, ascle::Float64, tol::Float64)
        @test AMOS.uchk(y,ascle,tol) == AMOS._uchk(y,ascle,tol)
    end

    test_tol = [
        SPECIAL_FLOAT32...,
        # 1e-0 ~ 1e-16
        [ 10.0^i for i in 0:-1:-16 ]...,
    ]
    test_ascle = [
        test_tol...,
    ]
    test_y = [
        complex(SPECIAL_FLOAT32)...,
        # [0, 1)+i[0, 1)
        rand(ComplexF64, 10)...,
        # (1+im) * 1e0 ~ 1e16
        [ (1.0+1im)*10^i for i in 0:16 ]...,
    ]

    for tol in gen_neg_inv(test_tol),
        ascle in gen_neg_inv(test_ascle),
        y in gen_phase4(test_y)
        tuchk(y, ascle, tol)
    end
end

"Pure julia log(gamma(z))"
function gammalog(z::Float64)
    if z > 0.0
        gam = gamma(z)
        if gam > 0.0
            return log(gam)
        end
    end

    return NaN
end

@testset "AMOS.gammaln" begin
    test_y = [
        SPECIAL_FLOAT32...,
        # [0, 1)
        rand(Float64, 10)...,
        # 1e-0 ~ 1e-16
        [ 10.0^i for i in 0:-1:-16 ]...,
        -1.0, -10., -Inf,
        NaN,
    ]

    for y in (test_y)
        ref_jl = gammalog(y)    # pure julia impl, ground truth
        ref = AMOS._gammaln(y)  # call AMOS, baseline
        res = AMOS.gammaln(y)   # fortran translation

        if isnan(ref)
            @test isnan(res)
            if isinf(ref_jl)
                # broken on 0x1.fffffep+127
                @test_broken isnan(ref_jl)
            else
                @test isnan(ref_jl)
            end
        else
            @test ref ≈ res
            @test ref == res

            if isinf(ref_jl)
                # broken on 0x1.fffffcp-127
                @test_broken ref_jl ≈ res
            else
                @test ref_jl ≈ res
            end
        end
    end

    @test_broken AMOS.gammaln(0x1.fffffep+127) == Inf
    @test_broken AMOS.gammaln(0x1.fffffcp-127) == gammalog(0x1.fffffcp-127)
end

@testset "AMOS.kscl!" begin
    function tkscl!(
        y::Vector{ComplexF64},
        zr::ComplexF64,
        fnu::Float64=2.0,
        n::Int=2,
        rz::ComplexF64=0.5+0.5*im,
        ascle::Float64=1e-30,
        tol::Float64=1e-15,
        elim::Float64=20.,
    )
        # @info "input" y zr fnu n rz ascle tol elim
        y_ref = copy(y)
        y_res = copy(y)
        ref = AMOS._kscl!(y_ref,zr,fnu,n,rz,ascle,tol,elim)
        res = AMOS.kscl!(y_res,zr,fnu,n,rz,ascle,tol,elim)

        @test ref == res
        if ref != res
            @info "ref != res" ref res
            @info "params" zr fnu n rz ascle elim
            @show y y_ref y_res
        end

        if contains_inf_nan(y_ref)
            # use === to test NaN/Inf
            @test all(y_ref .=== y_res)
        else
            @test isapprox(y_ref, y_res)
            # @info "isapprox"  ref zr y y_ref y_res
            # @show y_ref
            # @show y_res
            if !isapprox(y_ref, y_res)
                @info "!isapprox"  ref zr y y_ref y_res
                @info "params" zr fnu n rz ascle elim
                @show y y_ref y_res
            end
        end
    end

    test_y = [
        # ComplexF64[],  # TODO: n==0
        # n=1
        ComplexF64[0.0],
        ComplexF64[1.0],
        ComplexF64[pi],
        [ rand(ComplexF64, 1) for _ in 1:10 ]...,
        # n=2
        ComplexF64[0.0, 0.0],
        ComplexF64[1., 1.],
        ComplexF64[1., 0.],
        ComplexF64[0., 1.],
        ComplexF64[1., 2.],
        [ rand(ComplexF64, 2) for _ in 1:10 ]...,
        # n=5
        ComplexF64[ complex(i,i) for i in 1:5 ],
    ]
    test_zr = [
        SPECIAL_FLOAT32...,
        # [0, 1)
        rand(Float64, 10)...,
        # 1e-0 ~ 1e-16
        [ 10.0^i for i in 0:-1:-16 ]...,
        1.0,
        pi, ℯ,
    ]

    for y in gen_phase4(test_y),
        zr in gen_neg_inv(test_zr)
        # TODO: @test_broken
        isnan(zr) && continue
        isinf(zr) && continue
        0x1.fffffep+127==zr && continue

        tkscl!(y, complex(zr), 2.0, length(y))
    end
end
