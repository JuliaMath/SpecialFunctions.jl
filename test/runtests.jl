using SpecialFunctions
using Base.Test

if VERSION < v"0.5.0-dev+2955"
    rtoldefault{T<:Number,S<:Number}(x::Union{T,Type{T}}, y::Union{S,Type{S}}) =
        max(Base.rtoldefault(real(T)), Base.rtoldefault(real(S)))

    _isapprox(x, y; rtol::Real=rtoldefault(x, y), atol::Real=0) =
        Base.isapprox(x, y; rtol=rtol, atol=atol)

    # This will cause a warning
    const ≈ = _isapprox
end

include("bessel.jl")
