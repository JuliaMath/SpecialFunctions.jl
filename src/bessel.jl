# This file contains code that was formerly a part of Julia. License is MIT: http://julialang.org/license

using Base.Math: nan_dom_err

struct AmosException <: Exception
    id::Int32
end

function Base.showerror(io::IO, ex::AmosException)
    print(io, "AmosException with id $(ex.id): ")
    if ex.id == 0
        print(io, "normal return, computation complete.")
    elseif ex.id == 1
        print(io, "input error.")
    elseif ex.id == 2
        print(io, "overflow.")
    elseif ex.id == 3
        print(io, "input argument magnitude large, less than half machine accuracy loss by argument reduction.")
    elseif ex.id == 4
        print(io, "input argument magnitude too large, complete loss of accuracy by argument reduction.")
    elseif ex.id == 5
        print(io, "algorithm termination condition not met.")
    else
        print(io, "invalid error flag.")
    end
end

## Airy functions
function _airy(z::Complex{Float64}, id::Int32, kode::Int32)
    ai1, ai2 = Ref{Float64}(), Ref{Float64}()
    ae1, ae2 = Ref{Int32}(), Ref{Int32}()

    ccall((:zairy_,libopenspecfun), Cvoid,
          (Ref{Float64}, Ref{Float64}, Ref{Int32}, Ref{Int32},
           Ref{Float64}, Ref{Float64}, Ref{Int32}, Ref{Int32}),
           real(z), imag(z), id, kode,
           ai1, ai2, ae1, ae2)

    if ae2[] == 0 || ae2[] == 3 # ignore underflow and less than half machine accuracy loss
        return complex(ai1[], ai2[])
    else
        throw(AmosException(ae2[]))
    end
end

function _biry(z::Complex{Float64}, id::Int32, kode::Int32)
    ai1, ai2 = Ref{Float64}(), Ref{Float64}()
    ae1 = Ref{Int32}()

    ccall((:zbiry_,libopenspecfun), Cvoid,
          (Ref{Float64}, Ref{Float64}, Ref{Int32}, Ref{Int32},
           Ref{Float64}, Ref{Float64}, Ref{Int32}),
           real(z), imag(z), id, kode,
           ai1, ai2, ae1)

    if ae1[] == 0 || ae1[] == 3 # ignore less than half machine accuracy loss
        return complex(ai1[], ai2[])
    else
        throw(AmosException(ae1[]))
    end
end


"""
    airyai(x)

Airy function of the first kind ``\\operatorname{Ai}(x)``.

External links: [DLMF](https://dlmf.nist.gov/9.2), [Wikipedia](https://en.wikipedia.org/wiki/Airy_function)

See also: [`airyaix`](@ref), [`airyaiprime`](@ref), [`airybi`](@ref)
"""
function airyai end
airyai(z::Complex{Float64}) = _airy(z, Int32(0), Int32(1))

"""
    airyaiprime(x)

Derivative of the Airy function of the first kind ``\\operatorname{Ai}'(x)``.

External links: [DLMF](https://dlmf.nist.gov/9.2), [Wikipedia](https://en.wikipedia.org/wiki/Airy_function)

See also: [`airyaiprimex`](@ref), [`airyai`](@ref), [`airybi`](@ref)
"""
function airyaiprime end
airyaiprime(z::Complex{Float64}) =  _airy(z, Int32(1), Int32(1))

"""
    airybi(x)

Airy function of the second kind ``\\operatorname{Bi}(x)``.

External links: [DLMF](https://dlmf.nist.gov/9.2), [Wikipedia](https://en.wikipedia.org/wiki/Airy_function)

See also: [`airybix`](@ref), [`airybiprime`](@ref),  [`airyai`](@ref)
"""
function airybi end
airybi(z::Complex{Float64}) = _biry(z, Int32(0), Int32(1))

"""
    airybiprime(x)

Derivative of the Airy function of the second kind ``\\operatorname{Bi}'(x)``.

External links: [DLMF](https://dlmf.nist.gov/9.2), [Wikipedia](https://en.wikipedia.org/wiki/Airy_function)

See also: [`airybiprimex`](@ref), [`airybi`](@ref), [`airyai`](@ref)
"""
function airybiprime end
airybiprime(z::Complex{Float64}) = _biry(z, Int32(1), Int32(1))

"""
    airyaix(x)

Scaled Airy function of the first kind ``\\operatorname{Ai}(x) e^{\\frac{2}{3} x
\\sqrt{x}}``.  Throws `DomainError` for negative `Real` arguments.

External links: [DLMF](https://dlmf.nist.gov/9.2), [Wikipedia](https://en.wikipedia.org/wiki/Airy_function)

See also: [`airyai`](@ref), [`airyaiprime`](@ref), [`airybi`](@ref)
"""
function airyaix end
airyaix(z::Complex{Float64}) = _airy(z, Int32(0), Int32(2))

"""
    airyaiprimex(x)

Scaled derivative of the Airy function of the first kind ``\\operatorname{Ai}'(x)
e^{\\frac{2}{3} x \\sqrt{x}}``.  Throws `DomainError` for negative `Real` arguments.

External links: [DLMF](https://dlmf.nist.gov/9.2), [Wikipedia](https://en.wikipedia.org/wiki/Airy_function)

See also: [`airyaiprime`](@ref), [`airyai`](@ref), [`airybi`](@ref)
"""
function airyaiprimex end
airyaiprimex(z::Complex{Float64}) =  _airy(z, Int32(1), Int32(2))

"""
    airybix(x)

Scaled Airy function of the second kind ``\\operatorname{Bi}(x) e^{- \\left| \\operatorname{Re} \\left( \\frac{2}{3} x \\sqrt{x} \\right) \\right|}``.

External links: [DLMF](https://dlmf.nist.gov/9.2), [Wikipedia](https://en.wikipedia.org/wiki/Airy_function)

See also: [`airybi`](@ref), [`airybiprime`](@ref), [`airyai`](@ref)
"""
function airybix end
airybix(z::Complex{Float64}) = _biry(z, Int32(0), Int32(2))

"""
    airybiprimex(x)

Scaled derivative of the Airy function of the second kind ``\\operatorname{Bi}'(x) e^{- \\left| \\operatorname{Re} \\left( \\frac{2}{3} x \\sqrt{x} \\right) \\right|}``.

External links: [DLMF](https://dlmf.nist.gov/9.2), [Wikipedia](https://en.wikipedia.org/wiki/Airy_function)

See also: [`airybiprime`](@ref), [`airybi`](@ref), [`airyai`](@ref)
"""
function airybiprimex end
airybiprimex(z::Complex{Float64}) = _biry(z, Int32(1), Int32(2))

for afn in (:airyai, :airyaiprime, :airybi, :airybiprime,
            :airyaix, :airyaiprimex, :airybix, :airybiprimex)
    @eval begin
        $afn(z::Complex) = $afn(float(z))
        $afn(z::Complex{Float16}) = Complex{Float16}($afn(Complex{Float32}(z)))
        $afn(z::Complex{Float32}) = Complex{Float32}($afn(Complex{Float64}(z)))
        $afn(z::Complex{<:AbstractFloat}) = throw(MethodError($afn,(z,)))
    end
    if afn in (:airyaix, :airyaiprimex)
        @eval $afn(x::Real) = x < 0 ? throw(DomainError(x, "`x` must be nonnegative.")) : real($afn(complex(float(x))))
    else
        @eval $afn(x::Real) = real($afn(complex(float(x))))
    end
end

function airyai(x::BigFloat)
    z = BigFloat()
    ccall((:mpfr_ai, :libmpfr), Int32, (Ref{BigFloat}, Ref{BigFloat}, Int32), z, x, ROUNDING_MODE[])
    return z
end

## Bessel functions

# besselj0, besselj1, bessely0, bessely1
for jy in ("j","y"), nu in (0,1)
    jynu = Expr(:quote, Symbol(jy,nu))
    jynuf = Expr(:quote, Symbol(jy,nu,"f"))
    bjynu = Symbol("bessel",jy,nu)
    if jy == "y"
        @eval begin
            $bjynu(x::Float64) = nan_dom_err(ccall(($jynu,libopenlibm),  Float64, (Float64,), x), x)
            $bjynu(x::Float32) = nan_dom_err(ccall(($jynuf,libopenlibm), Float32, (Float32,), x), x)
            $bjynu(x::Float16) = Float16($bjynu(Float32(x)))
        end
    else
        @eval begin
            $bjynu(x::Float64) = ccall(($jynu,libopenlibm),  Float64, (Float64,), x)
            $bjynu(x::Float32) = ccall(($jynuf,libopenlibm), Float32, (Float32,), x)
            $bjynu(x::Float16) = Float16($bjynu(Float32(x)))
        end
    end
    @eval begin
        $bjynu(x::Real) = $bjynu(float(x))
        $bjynu(x::Complex) = $(Symbol("bessel",jy))($nu,x)
    end
end


function _besselh(nu::Float64, k::Int32, z::Complex{Float64}, kode::Int32)
    ai1, ai2 = Ref{Float64}(), Ref{Float64}()
    ae1, ae2 = Ref{Int32}(), Ref{Int32}()

    ccall((:zbesh_,libopenspecfun), Cvoid,
           (Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Int32}, Ref{Int32}, Ref{Int},
            Ref{Float64}, Ref{Float64}, Ref{Int32}, Ref{Int32}),
            real(z), imag(z), nu, kode, k, 1,
            ai1, ai2, ae1, ae2)

    if ae2[] == 0 || ae2[] == 3
        return complex(ai1[], ai2[])
    else
        throw(AmosException(ae2[]))
    end
end

function _besseli(nu::Float64, z::Complex{Float64}, kode::Int32)
    ai1, ai2 = Ref{Float64}(), Ref{Float64}()
    ae1, ae2 = Ref{Int32}(), Ref{Int32}()

    ccall((:zbesi_,libopenspecfun), Cvoid,
          (Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Int32}, Ref{Int32},
           Ref{Float64}, Ref{Float64}, Ref{Int32}, Ref{Int32}),
           real(z), imag(z), nu, kode, 1,
           ai1, ai2, ae1, ae2)

    if ae2[] == 0 || ae2[] == 3
        return complex(ai1[], ai2[])
    else
        throw(AmosException(ae2[]))
    end
end

function _besselj(nu::Float64, z::Complex{Float64}, kode::Int32)
    ai1, ai2 = Ref{Float64}(), Ref{Float64}()
    ae1, ae2 = Ref{Int32}(), Ref{Int32}()

    ccall((:zbesj_,libopenspecfun), Cvoid,
          (Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Int32}, Ref{Int32},
           Ref{Float64}, Ref{Float64}, Ref{Int32}, Ref{Int32}),
           real(z), imag(z), nu, kode, 1,
           ai1, ai2, ae1, ae2)

    if ae2[] == 0 || ae2[] == 3
        return complex(ai1[], ai2[])
    else
        throw(AmosException(ae2[]))
    end
end

function _besselk(nu::Float64, z::Complex{Float64}, kode::Int32)
    ai1, ai2 = Ref{Float64}(), Ref{Float64}()
    ae1, ae2 = Ref{Int32}(), Ref{Int32}()

    ccall((:zbesk_,libopenspecfun), Cvoid,
          (Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Int32}, Ref{Int32},
           Ref{Float64}, Ref{Float64}, Ref{Int32}, Ref{Int32}),
           real(z), imag(z), nu, kode, 1,
           ai1, ai2, ae1, ae2)

    if ae2[] == 0 || ae2[] == 3
        return complex(ai1[], ai2[])
    else
        throw(AmosException(ae2[]))
    end
end

function _bessely(nu::Float64, z::Complex{Float64}, kode::Int32)
    ai1, ai2 = Ref{Float64}(), Ref{Float64}()
    ae1, ae2 = Ref{Int32}(), Ref{Int32}()
    wrk1, wrk2 = Ref{Float64}(), Ref{Float64}()

    ccall((:zbesy_,libopenspecfun), Cvoid,
          (Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Int32}, Ref{Int32},
           Ref{Float64}, Ref{Float64}, Ref{Int32}, Ref{Float64}, Ref{Float64}, Ref{Int32}),
           real(z), imag(z), nu, kode, 1,
           ai1, ai2, ae1, wrk1, wrk2, ae2)

    if ae2[] == 0 || ae2[] == 3
        return complex(ai1[], ai2[])
    else
        throw(AmosException(ae2[]))
    end
end

"""
    besselh(nu, [k=1,] x)

Bessel function of the third kind of order `nu` (the Hankel function). `k` is either 1 or 2,
selecting [`hankelh1`](@ref) or [`hankelh2`](@ref), respectively.
`k` defaults to 1 if it is omitted.

External links: [DLMF](https://dlmf.nist.gov/10.2.5) and [DLMF](https://dlmf.nist.gov/10.2.6), [Wikipedia](https://en.wikipedia.org/wiki/Bessel_function#Hankel_functions:_H(1)%CE%B1,_H(2)%CE%B1)

See also: [`besselhx`](@ref) for an exponentially scaled variant.
"""
function besselh end

function besselh(nu::Float64, k::Integer, z::Complex{Float64})
    if nu < 0
        s = (k == 1) ? 1 : -1
        return _besselh(-nu,Int32(k),z,Int32(1)) * complex(cospi(nu),-s*sinpi(nu))
    end
    return _besselh(nu,Int32(k),z,Int32(1))
end

function besselh(nu::Float64, k::Integer, x::Float64)
    # Given that x is real, Jnu(x) and Ynu(x) are also real.
    if k == 1
        return complex(besselj(nu, x), bessely(nu, x))
    elseif k == 2
        return complex(besselj(nu, x), -bessely(nu, x))
    else
        # We emulate ZBESH's behaviour
        throw(AmosException(1))
    end
end

"""
    besselhx(nu, [k=1,] z)

Compute the scaled Hankel function ``\\exp(∓iz) H_ν^{(k)}(z)``, where
``k`` is 1 or 2, ``H_ν^{(k)}(z)`` is `besselh(nu, k, z)`, and ``∓`` is
``-`` for ``k=1`` and ``+`` for ``k=2``.  `k` defaults to 1 if it is omitted.

The reason for this function is that ``H_ν^{(k)}(z)`` is asymptotically
proportional to ``\\exp(∓iz)/\\sqrt{z}`` for large ``|z|``, and so the
[`besselh`](@ref) function is susceptible to overflow or underflow
when `z` has a large imaginary part.  The `besselhx` function cancels this
exponential factor (analytically), so it avoids these problems.

External links: [DLMF](https://dlmf.nist.gov/10.2), [Wikipedia](https://en.wikipedia.org/wiki/Bessel_function#Hankel_functions:_H(1)%CE%B1,_H(2)%CE%B1)

See also: [`besselh`](@ref)
"""
function besselhx end

function besselhx(nu::Float64, k::Integer, z::Complex{Float64})
    if nu < 0
        s = (k == 1) ? 1 : -1
        return _besselh(-nu,Int32(k),z,Int32(2)) * complex(cospi(nu),-s*sinpi(nu))
    end
    return _besselh(nu,Int32(k),z,Int32(2))
end

function besseli(nu::Float64, z::Complex{Float64})
    if nu < 0
        if isinteger(nu)
            return _besseli(-nu,z,Int32(1))
        else
            return _besseli(-nu,z,Int32(1)) - 2_besselk(-nu,z,Int32(1))*sinpi(nu)/pi
        end
    else
        return _besseli(nu,z,Int32(1))
    end
end

function besselix(nu::Float64, z::Complex{Float64})
    if nu < 0
        if isinteger(nu)
            return _besseli(-nu,z,Int32(2))
        else
            return _besseli(-nu,z,Int32(2)) - 2_besselk(-nu,z,Int32(2))*exp(-abs(real(z))-z)*sinpi(nu)/pi
        end
    else
        return _besseli(nu,z,Int32(2))
    end
end

function besselj(nu::Float64, z::Complex{Float64})
    if nu < 0
        if isinteger(nu)
            return _besselj(-nu,z,Int32(1))*cospi(nu)
        else
            return _besselj(-nu,z,Int32(1))*cospi(nu) + _bessely(-nu,z,Int32(1))*sinpi(nu)
        end
    else
        return _besselj(nu,z,Int32(1))
    end
end

besselj(nu::Cint, x::Float64) = ccall((:jn, libopenlibm), Float64, (Cint, Float64), nu, x)
besselj(nu::Cint, x::Float32) = ccall((:jnf, libopenlibm), Float32, (Cint, Float32), nu, x)


function besseljx(nu::Float64, z::Complex{Float64})
    if nu < 0
        if isinteger(nu)
            return _besselj(-nu,z,Int32(2))*cospi(nu)
        else
            return _besselj(-nu,z,Int32(2))*cospi(nu) + _bessely(-nu,z,Int32(2))*sinpi(nu)
        end
    else
        return _besselj(nu,z,Int32(2))
    end
end

besselk(nu::Float64, z::Complex{Float64}) = _besselk(abs(nu), z, Int32(1))

besselkx(nu::Float64, z::Complex{Float64}) = _besselk(abs(nu), z, Int32(2))

function bessely(nu::Cint, x::Float64)
    if x < 0
        throw(DomainError(x, "`x` must be nonnegative."))
    end
    ccall((:yn, libopenlibm), Float64, (Cint, Float64), nu, x)
end
function bessely(nu::Cint, x::Float32)
    if x < 0
        throw(DomainError(x, "`x` must be nonnegative."))
    end
    ccall((:ynf, libopenlibm), Float32, (Cint, Float32), nu, x)
end

function bessely(nu::Float64, z::Complex{Float64})
    if nu < 0
        return _bessely(-nu,z,Int32(1))*cospi(nu) - _besselj(-nu,z,Int32(1))*sinpi(nu)
    else
        return _bessely(nu,z,Int32(1))
    end
end

function besselyx(nu::Float64, z::Complex{Float64})
    if nu < 0
        return _bessely(-nu,z,Int32(2))*cospi(nu) - _besselj(-nu,z,Int32(2))*sinpi(nu)
    else
        return _bessely(nu,z,Int32(2))
    end
end

"""
    besseli(nu, x)

Modified Bessel function of the first kind of order `nu`, ``I_\\nu(x)``.

External links: [DLMF](https://dlmf.nist.gov/10.25.2), [Wikipedia](https://en.wikipedia.org/wiki/Bessel_function#Modified_Bessel_functions:_I%CE%B1,_K%CE%B1)

See also: [`besselix(nu,x)`](@ref SpecialFunctions.besselix), [`besselj(nu,x)`](@ref SpecialFunctions.besselj), [`besselk(nu,x)`](@ref SpecialFunctions.besselk)
"""
function besseli(nu::Real, x::AbstractFloat)
    if x < 0 && !isinteger(nu)
        throw(DomainError(x, "`x` must be nonnegative and `nu` must be an integer."))
    end
    real(besseli(float(nu), complex(x)))
end

"""
    besselix(nu, x)

Scaled modified Bessel function of the first kind of order `nu`, ``I_\\nu(x) e^{- | \\operatorname{Re}(x) |}``.

External links: [DLMF](https://dlmf.nist.gov/10.25.2), [Wikipedia](https://en.wikipedia.org/wiki/Bessel_function#Modified_Bessel_functions:_I%CE%B1,_K%CE%B1)

See also: [`besseli(nu,x)`](@ref SpecialFunctions.besseli), [`besselj(nu,x)`](@ref SpecialFunctions.besselj), [`besselk(nu,x)`](@ref SpecialFunctions.besselk)
"""
function besselix(nu::Real, x::AbstractFloat)
    if x < 0 && !isinteger(nu)
        throw(DomainError(x, "`x` must be nonnegative and `nu` must be an integer."))
    end
    real(besselix(float(nu), complex(x)))
end

"""
    besselj(nu, x)

Bessel function of the first kind of order `nu`, ``J_\\nu(x)``.

External links: [DLMF](https://dlmf.nist.gov/10.2.2), [Wikipedia](https://en.wikipedia.org/wiki/Bessel_function#Bessel_functions_of_the_first_kind:_J%CE%B1)

See also: [`besseljx(nu,x)`](@ref SpecialFunctions.besseljx), [`besseli(nu,x)`](@ref SpecialFunctions.besseli), [`besselk(nu,x)`](@ref SpecialFunctions.besselk)
"""
function besselj(nu::Real, x::AbstractFloat)
    if isinteger(nu)
        if typemin(Cint) <= nu <= typemax(Cint)
            return besselj(Cint(nu), x)
        end
    elseif x < 0
        throw(DomainError(x, "`x` must be nonnegative."))
    end
    real(besselj(float(nu), complex(x)))
end

"""
    besseljx(nu, x)

Scaled Bessel function of the first kind of order `nu`, ``J_\\nu(x) e^{- | \\operatorname{Im}(x) |}``.

External links: [DLMF](https://dlmf.nist.gov/10.2.2), [Wikipedia](https://en.wikipedia.org/wiki/Bessel_function#Bessel_functions_of_the_first_kind:_J%CE%B1)

See also: [`besselj(nu,x)`](@ref SpecialFunctions.besselj), [`besseli(nu,x)`](@ref SpecialFunctions.besseli), [`besselk(nu,x)`](@ref SpecialFunctions.besselk)
"""
function besseljx(nu::Real, x::AbstractFloat)
    if x < 0 && !isinteger(nu)
        throw(DomainError(x, "`x` must be nonnegative and `nu` must be an integer."))
    end
    real(besseljx(float(nu), complex(x)))
end

"""
    besselk(nu, x)

Modified Bessel function of the second kind of order `nu`, ``K_\\nu(x)``.

External links: [DLMF](https://dlmf.nist.gov/10.25.3), [Wikipedia](https://en.wikipedia.org/wiki/Bessel_function#Modified_Bessel_functions:_I%CE%B1,_K%CE%B1)

See also: See also: [`besselkx(nu,x)`](@ref SpecialFunctions.besselkx), [`besseli(nu,x)`](@ref SpecialFunctions.besseli), [`besselj(nu,x)`](@ref SpecialFunctions.besselj)
"""
function besselk(nu::Real, x::AbstractFloat)
    if x < 0
        throw(DomainError(x, "`x` must be nonnegative."))
    elseif x == 0
        return oftype(x, Inf)
    end
    real(besselk(float(nu), complex(x)))
end

"""
    besselkx(nu, x)

Scaled modified Bessel function of the second kind of order `nu`, ``K_\\nu(x) e^x``.

External links: [DLMF](https://dlmf.nist.gov/10.25.3), [Wikipedia](https://en.wikipedia.org/wiki/Bessel_function#Modified_Bessel_functions:_I%CE%B1,_K%CE%B1)

See also: [`besselk(nu,x)`](@ref SpecialFunctions.besselk), [`besseli(nu,x)`](@ref SpecialFunctions.besseli), [`besselj(nu,x)`](@ref SpecialFunctions.besselj)
"""
function besselkx(nu::Real, x::AbstractFloat)
    if x < 0
        throw(DomainError(x, "`x` must be nonnegative."))
    elseif x == 0
        return oftype(x, Inf)
    end
    real(besselkx(float(nu), complex(x)))
end

"""
    bessely(nu, x)

Bessel function of the second kind of order `nu`, ``Y_\\nu(x)``.

External links: [DLMF](https://dlmf.nist.gov/10.2.3), [Wikipedia](https://en.wikipedia.org/wiki/Bessel_function#Bessel_functions_of_the_second_kind:_Y%CE%B1)

See also [`besselyx(nu,x)`](@ref SpecialFunctions.besselyx) for a scaled variant.
"""
function bessely(nu::Real, x::AbstractFloat)
    if x < 0
        throw(DomainError(x, "`x` must be nonnegative."))
    elseif isinteger(nu) && typemin(Cint) <= nu <= typemax(Cint)
        return bessely(Cint(nu), x)
    end
    real(bessely(float(nu), complex(x)))
end

"""
    besselyx(nu, x)

Scaled Bessel function of the second kind of order `nu`,
``Y_\\nu(x) e^{- | \\operatorname{Im}(x) |}``.

External links: [DLMF](https://dlmf.nist.gov/10.2.3), [Wikipedia](https://en.wikipedia.org/wiki/Bessel_function#Bessel_functions_of_the_second_kind:_Y%CE%B1)

See also [`bessely(nu,x)`](@ref SpecialFunctions.bessely)
"""
function besselyx(nu::Real, x::AbstractFloat)
    if x < 0
        throw(DomainError(x, "`x` must be nonnegative."))
    end
    real(besselyx(float(nu), complex(x)))
end

for f in ("i", "ix", "j", "jx", "k", "kx", "y", "yx")
    bfn = Symbol("bessel", f)
    @eval begin
        $bfn(nu::Real, x::Real) = $bfn(nu, float(x))
        function $bfn(nu::Real, z::Complex)
            Tf = promote_type(float(typeof(nu)),float(typeof(real(z))))
            $bfn(Tf(nu), Complex{Tf}(z))
        end
        $bfn(nu::Float16, x::Complex{Float16}) = Complex{Float16}($bfn(Float32(nu), Complex{Float32}(x)))
        $bfn(nu::Float32, x::Complex{Float32}) = Complex{Float32}($bfn(Float64(nu), Complex{Float64}(x)))
        $bfn(k::T, z::Complex{T}) where {T<:AbstractFloat} = throw(MethodError($bfn,(k,z)))
    end
end


for bfn in (:besselh, :besselhx)
    @eval begin
        $bfn(nu, z) = $bfn(nu, 1, z)
        $bfn(nu::Real, k::Integer, x::Real) = $bfn(float(nu), k, float(x))
        $bfn(nu::AbstractFloat, k::Integer, x::AbstractFloat) = $bfn(float(nu), k, complex(x))
        function $bfn(nu::Real, k::Integer, z::Complex)
            Tf = promote_type(float(typeof(nu)),float(typeof(real(z))))
            $bfn(Tf(nu), k, Complex{Tf}(z))
        end
        $bfn(nu::Float16, k::Integer, x::Complex{Float16}) = Complex{Float16}($bfn(Float32(nu), k, Complex{Float32}(x)))
        $bfn(nu::Float32, k::Integer, x::Complex{Float32}) = Complex{Float32}($bfn(Float64(nu), k, Complex{Float64}(x)))
        $bfn(nu::T, k::Integer, z::Complex{T}) where {T<:AbstractFloat} = throw(MethodError($bfn,(nu,k,z)))
    end
end

besselh(nu::Float16, k::Integer, x::Float16) = Complex{Float16}(besselh(Float32(nu), k, Float32(x)))
besselh(nu::Float32, k::Integer, x::Float32) = Complex{Float32}(besselh(Float64(nu), k, Float64(x)))

"""
    besselj0(x)

Bessel function of the first kind of order 0, ``J_0(x)``.

External links: [DLMF](https://dlmf.nist.gov/10.2.2), [Wikipedia](https://en.wikipedia.org/wiki/Bessel_function#Bessel_functions_of_the_first_kind:_J%CE%B1)

See also: [`besselj(nu,x)`](@ref SpecialFunctions.besselj)
"""
function besselj0(x::BigFloat)
    z = BigFloat()
    ccall((:mpfr_j0, :libmpfr), Int32, (Ref{BigFloat}, Ref{BigFloat}, Int32), z, x, ROUNDING_MODE[])
    return z
end

"""
    besselj1(x)

Bessel function of the first kind of order 1, ``J_1(x)``.

External links: [DLMF](https://dlmf.nist.gov/10.2.2), [Wikipedia](https://en.wikipedia.org/wiki/Bessel_function#Bessel_functions_of_the_first_kind:_J%CE%B1)

See also: [`besselj(nu,x)`](@ref SpecialFunctions.besselj)
"""
function besselj1(x::BigFloat)
    z = BigFloat()
    ccall((:mpfr_j1, :libmpfr), Int32, (Ref{BigFloat}, Ref{BigFloat}, Int32), z, x, ROUNDING_MODE[])
    return z
end

function besselj(n::Integer, x::BigFloat)
    z = BigFloat()
    ccall((:mpfr_jn, :libmpfr), Int32, (Ref{BigFloat}, Clong, Ref{BigFloat}, Int32), z, n, x, ROUNDING_MODE[])
    return z
end

"""
    bessely0(x)

Bessel function of the second kind of order 0, ``Y_0(x)``.

External links: [DLMF](https://dlmf.nist.gov/10.2.3), [Wikipedia](https://en.wikipedia.org/wiki/Bessel_function#Bessel_functions_of_the_second_kind:_Y%CE%B1)

See also: [`bessely(nu,x)`](@ref SpecialFunctions.bessely)
"""
function bessely0(x::BigFloat)
    if x < 0
        throw(DomainError(x, "`x` must be nonnegative."))
    end
    z = BigFloat()
    ccall((:mpfr_y0, :libmpfr), Int32, (Ref{BigFloat}, Ref{BigFloat}, Int32), z, x, ROUNDING_MODE[])
    return z
end

"""
    bessely1(x)

Bessel function of the second kind of order 1, ``Y_1(x)``.

External links: [DLMF](https://dlmf.nist.gov/10.2.3), [Wikipedia](https://en.wikipedia.org/wiki/Bessel_function#Bessel_functions_of_the_second_kind:_Y%CE%B1)

See also: [`bessely(nu,x)`](@ref SpecialFunctions.bessely)
"""
function bessely1(x::BigFloat)
    if x < 0
        throw(DomainError(x, "`x` must be nonnegative."))
    end
    z = BigFloat()
    ccall((:mpfr_y1, :libmpfr), Int32, (Ref{BigFloat}, Ref{BigFloat}, Int32), z, x, ROUNDING_MODE[])
    return z
end

function bessely(n::Integer, x::BigFloat)
    if x < 0
        throw(DomainError(x, "`x` must be nonnegative."))
    end
    z = BigFloat()
    ccall((:mpfr_yn, :libmpfr), Int32, (Ref{BigFloat}, Clong, Ref{BigFloat}, Int32), z, n, x, ROUNDING_MODE[])
    return z
end

"""
    sphericalbesselj(nu, x)

Spherical bessel function of the first kind at order `nu`, ``j_ν(x)``. This is the non-singular
solution to the radial part of the Helmholz equation in spherical coordinates.
"""
function sphericalbesselj(nu, x::T) where {T}
    besselj_nuhalf_x = besselj(nu + one(nu)/2, x)
    if abs(x) ≤ sqrt(eps(real(zero(besselj_nuhalf_x))))
        nu == 0 ? one(besselj_nuhalf_x) : zero(besselj_nuhalf_x)
    else
        √((float(T))(π)/2x) * besselj_nuhalf_x
    end
end

"""
    sphericalbessely(nu, x)

Spherical bessel function of the second kind at order `nu`, ``y_ν(x)``. This is the singular
solution to the radial part of the Helmholz equation in spherical coordinates. Sometimes
known as a spherical Neumann function.
"""
sphericalbessely(nu, x::T) where {T} = √((float(T))(π)/2x) * bessely(nu + one(nu)/2, x)

"""
    hankelh1(nu, x)

Bessel function of the third kind of order `nu`, ``H^{(1)}_\\nu(x)``.

External links: [DLMF](https://dlmf.nist.gov/10.2.5), [Wikipedia](https://en.wikipedia.org/wiki/Bessel_function#Hankel_functions:_H(1)%CE%B1,_H(2)%CE%B1)

See also: [`hankelh1x`](@ref SpecialFunctions.hankelh1x)
"""
hankelh1(nu, z) = besselh(nu, 1, z)

"""
    hankelh2(nu, x)

Bessel function of the third kind of order `nu`, ``H^{(2)}_\\nu(x)``.

External links: [DLMF](https://dlmf.nist.gov/10.2.6), [Wikipedia](https://en.wikipedia.org/wiki/Bessel_function#Hankel_functions:_H(1)%CE%B1,_H(2)%CE%B1)

See also: [`hankelh2x(nu,x)`](@ref SpecialFunctions.hankelh2x)
"""
hankelh2(nu, z) = besselh(nu, 2, z)

"""
    hankelh1x(nu, x)

Scaled Bessel function of the third kind of order `nu`, ``H^{(1)}_\\nu(x) e^{-x i}``.

External links: [DLMF](https://dlmf.nist.gov/10.2.5), [Wikipedia](https://en.wikipedia.org/wiki/Bessel_function#Hankel_functions:_H(1)%CE%B1,_H(2)%CE%B1)

See also: [`hankelh1`](@ref SpecialFunctions.hankelh1)
"""
hankelh1x(nu, z) = besselhx(nu, 1, z)

"""
    hankelh2x(nu, x)

Scaled Bessel function of the third kind of order `nu`, ``H^{(2)}_\\nu(x) e^{x i}``.

External links: [DLMF](https://dlmf.nist.gov/10.2.6), [Wikipedia](https://en.wikipedia.org/wiki/Bessel_function#Hankel_functions:_H(1)%CE%B1,_H(2)%CE%B1)

See also: [`hankelh2(nu,x)`](@ref SpecialFunctions.hankelh2)
"""
hankelh2x(nu, z) = besselhx(nu, 2, z)

"""
    jinc(x)

Bessel function of the first kind divided by `x`.
Following convention: ``\\operatorname{jinc}{x} = \\frac{2 \\cdot J_1{\\pi x}}{\\pi x}``.
Sometimes known as sombrero or besinc function.

External links: [Wikipedia](https://en.wikipedia.org/wiki/Sombrero_function)
"""
jinc(x::Number) = _jinc(float(x))
_jinc(x::Number) = iszero(x) ? one(x) : _jinc_core(x)
_jinc(x::Float16) = Float16(_jinc(Float32(x)))
_jinc(x::ComplexF16) = ComplexF16(_jinc(ComplexF32(x)))

# below these thresholds we evaluate a fourth order Taylor polynomial
_jinc_threshold(::Type{Float64}) = 0.002
_jinc_threshold(::Type{Float32}) = 0.05f0


 # for small arguments, a Taylor series is faster
@inline function _jinc(x::Union{T,Complex{T}}) where {T<:Union{Float32,Float64}}
    if fastabs(x) < _jinc_threshold(T)
        return @evalpoly(x^2, T(1), -T(π)^2/8, T(π)^4/192)
    else
        return _jinc_core(x)
    end
end
 # the actual definition of jinc
_jinc_core(x::Number) =  2 * besselj1(π*x) / (π*x)
