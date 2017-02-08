# This file contains code that was formerly a part of Julia. License is MIT: http://julialang.org/license

using Base: @deprecate, depwarn
using Compat: @dep_vectorize_1arg, @dep_vectorize_2arg

for f in (:digamma, :trigamma, :zeta, :eta, :erfcx, :erfi, :dawson, :airyai, :airyaiprime,
          :airybi, :airybiprime, :besselj0, :besselj1, :bessely0, :bessely1, :erf, :erfc)
    @eval @dep_vectorize_1arg Number $f
end

for f in (:invdigamma, :erfinc, :erfcinv)
    @eval @dep_vectorize_1arg Real $f
end

for f in (:polygamma, :zeta, :besseli, :besselix, :besselj, :besseljx, :besselk, :besselkx,
          :bessely, :besselyx, :besselh, :besselhx, :hankelh1, :hankelh2, :hankelh1x, :hankelh2x)
    @eval @dep_vectorize_2arg Number $f
end

@deprecate airy(z::Number) airyai(z)
@deprecate airyx(z::Number) airyaix(z)
@deprecate airyprime(z::Number) airyaiprime(z)
@deprecate airy{T<:Number}(x::AbstractArray{T}) airyai.(x)
@deprecate airyx{T<:Number}(x::AbstractArray{T}) airyaix.(x)
@deprecate airyprime{T<:Number}(x::AbstractArray{T}) airyprime.(x)

function _airy(k::Integer, z::Complex128)
    depwarn("`airy(k,x)` is deprecated, use `airyai(x)`, `airyaiprime(x)`, `airybi(x)` or `airybiprime(x)` instead.",:airy)
    id = Int32(k==1 || k==3)
    if k == 0 || k == 1
        return _airy(z, id, Int32(1))
    elseif k == 2 || k == 3
        return _biry(z, id, Int32(1))
    else
        throw(ArgumentError("k must be between 0 and 3"))
    end
end
function _airyx(k::Integer, z::Complex128)
    depwarn("`airyx(k,x)` is deprecated, use `airyaix(x)`, `airyaiprimex(x)`, `airybix(x)` or `airybiprimex(x)` instead.",:airyx)
    id = Int32(k==1 || k==3)
    if k == 0 || k == 1
        return _airy(z, id, Int32(2))
    elseif k == 2 || k == 3
        return _biry(z, id, Int32(2))
    else
        throw(ArgumentError("k must be between 0 and 3"))
    end
end

for afn in (:airy,:airyx)
    _afn = Symbol("_"*string(afn))
    suf  = string(afn)[5:end]
    @eval begin
        function $afn(k::Integer, z::Complex128)
            afn = $(QuoteNode(afn))
            suf = $(QuoteNode(suf))
            depwarn("`$afn(k,x)` is deprecated, use `airyai$suf(x)`, `airyaiprime$suf(x)`, `airybi$suf(x)` or `airybiprime$suf(x)` instead.",$(QuoteNode(afn)))
            $_afn(k,z)
        end

        $afn(k::Integer, z::Complex) = $afn(k, float(z))
        $afn{T<:AbstractFloat}(k::Integer, z::Complex{T}) = throw(MethodError($afn,(k,z)))
        $afn(k::Integer, z::Complex64) = Complex64($afn(k, Complex128(z)))
        $afn(k::Integer, x::Real) = $afn(k, float(x))
        $afn(k::Integer, x::AbstractFloat) = real($afn(k, complex(x)))

        function $afn{T<:Number}(k::Number, x::AbstractArray{T})
            $afn.(k,x)
        end
        function $afn{S<:Number}(k::AbstractArray{S}, x::Number)
            $afn.(k,x)
        end
        function $afn{S<:Number,T<:Number}(k::AbstractArray{S}, x::AbstractArray{T})
            $afn.(k,x)
        end
    end
end
