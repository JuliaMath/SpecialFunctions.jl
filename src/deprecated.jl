# This file contains code that was formerly a part of Julia. License is MIT: http://julialang.org/license

@deprecate airy(z::Number) airyai(z)
@deprecate airyx(z::Number) airyaix(z)
@deprecate airyprime(z::Number) airyaiprime(z)

function _airy(k::Integer, z::Complex128)
    Base.depwarn("`airy(k,x)` is deprecated, use `airyai(x)`, `airyaiprime(x)`, `airybi(x)` or `airybiprime(x)` instead.",:airy)
    id = Int32(k==1 || k==3)
    if k == 0 || k == 1
        return _airy(z, id, Int32(1))
    elseif k == 2 || k == 3
        return _biry(z, id, Int32(1))
    else
        throw(DomainError(k, "`k` must be between 0 and 3."))
    end
end
function _airyx(k::Integer, z::Complex128)
    Base.depwarn("`airyx(k,x)` is deprecated, use `airyaix(x)`, `airyaiprimex(x)`, `airybix(x)` or `airybiprimex(x)` instead.",:airyx)
    id = Int32(k==1 || k==3)
    if k == 0 || k == 1
        return _airy(z, id, Int32(2))
    elseif k == 2 || k == 3
        return _biry(z, id, Int32(2))
    else
        throw(DomainError(k, "`k` must be between 0 and 3."))
    end
end

for afn in (:airy,:airyx)
    _afn = Symbol("_"*string(afn))
    suf  = string(afn)[5:end]
    @eval begin
        function $afn(k::Integer, z::Complex128)
            afn = $(QuoteNode(afn))
            suf = $(QuoteNode(suf))
            Base.depwarn("`$afn(k,x)` is deprecated, use `airyai$suf(x)`, `airyaiprime$suf(x)`, `airybi$suf(x)` or `airybiprime$suf(x)` instead.",$(QuoteNode(afn)))
            $_afn(k,z)
        end

        $afn(k::Integer, z::Complex) = $afn(k, float(z))
        $afn(k::Integer, z::Complex{<:AbstractFloat}) = throw(MethodError($afn,(k,z)))
        $afn(k::Integer, z::Complex64) = Complex64($afn(k, Complex128(z)))
        $afn(k::Integer, x::Real) = $afn(k, float(x))
        $afn(k::Integer, x::AbstractFloat) = real($afn(k, complex(x)))
    end
end
