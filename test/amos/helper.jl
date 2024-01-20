
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
    zarr[ (+), (-) ]
"""
function gen_neg_inv(zarr::Vector{Float64})
    vcat(
        zarr,
        -1.0 * zarr,
        inv.(zarr)
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
