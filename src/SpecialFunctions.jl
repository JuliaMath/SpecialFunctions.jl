module SpecialFunctions

# export
#     airy,
#     airyai,
#     airyaiprime,
#     airybi,
#     airybiprime,
#     airyprime,
#     airyx,
#     besselh,
#     besseli,
#     besselix,
#     besselj,
#     besselj0,
#     besselj1,
#     besseljx,
#     besselk,
#     besselkx,
#     bessely,
#     bessely0,
#     bessely1,
#     besselyx,
#     hankelh1,
#     hankelh1x,
#     hankelh2,
#     hankelh2x

include("amos/Amos.jl")

using Base.Math.libm

const int32 = Int32
const float32 = Float32
const float64 = Float64
const complex64 = Complex64
const complex128 = Complex128

type AmosException <: Exception
    info::Int32
end

## Airy functions

function _airy(z::Complex128, id::Int32, kode::Int32)
    nz = ierr = Int32(0)
    (air,aii,nz,ierr) = Amos.ZAIRY(real(z), imag(z), id, kode, 0., 0., nz, ierr)
    if ierr == 0 || ierr == 3
        return complex(air,aii)
    else
        throw(AmosException(ierr))
    end
end

function _biry(z::Complex128, id::Int32, kode::Int32)
    (air,aii,ierr) = Amos.ZBIRY(real(z), imag(z), id, kode, 0., 0., Int32(0))
    if ierr == 0 || ierr == 3
        return complex(air,aii)
    else
        throw(AmosException(ierr))
    end
end

function airy(k::Int, z::Complex128)
    id = Int32(k==1 || k==3)
    if k == 0 || k == 1
        return _airy(z, id, Int32(1))
    elseif k == 2 || k == 3
        return _biry(z, id, Int32(1))
    else
        throw(ArgumentError("k must be between 0 and 3, got $k"))
    end
end

airy(z) = airy(0,z)
@vectorize_1arg Number airy
airyprime(z) = airy(1,z)
@vectorize_1arg Number airyprime
airyai(z) = airy(0,z)
@vectorize_1arg Number airyai
airyaiprime(z) = airy(1,z)
@vectorize_1arg Number airyaiprime
airybi(z) = airy(2,z)
@vectorize_1arg Number airybi
airybiprime(z) = airy(3,z)
@vectorize_1arg Number airybiprime

function airyx(k::Int, z::Complex128)
    id = Int32(k==1 || k==3)
    if k == 0 || k == 1
        return _airy(z, id, Int32(2))
    elseif k == 2 || k == 3
        return _biry(z, id, Int32(2))
    else
        throw(ArgumentError("k must be between 0 and 3, got $k"))
    end
end

for afn in (:airy, :airyx)
    @eval begin
        $afn(k::Integer, z::Complex) = $afn(k, float(z))
        $afn{T<:AbstractFloat}(k::Integer, z::Complex{T}) = throw(MethodError($afn,(k,z)))
        $afn(k::Integer, z::Complex64) = Complex64($afn(k, Complex128(z)))

        $afn(k::Integer, x::Real) = $afn(k, float(x))
        $afn(k::Integer, x::AbstractFloat) = real($afn(k, complex(x)))

        $afn(z) = $afn(0, z)
        @vectorize_1arg Number $afn
        @vectorize_2arg Number $afn
    end
end

## Bessel functions

# besselj0, besselj1, bessely0, bessely1
for jy in ("j", "y"), nu in (0, 1)
    jynu = Expr(:quote, Symbol(string(jy, nu)))
    jynuf = Expr(:quote, Symbol(string(jy, nu, "f")))
    bjynu = Symbol(string("bessel", jy, nu))
    if jy == "y"
        @eval begin
            $bjynu(x::Float64) = nan_dom_err(ccall(($jynu, libm),  Float64, (Float64,), x), x)
            $bjynu(x::Float32) = nan_dom_err(ccall(($jynuf, libm), Float32, (Float32,), x), x)
        end
    else
        @eval begin
            $bjynu(x::Float64) = ccall(($jynu, libm),  Float64, (Float64,), x)
            $bjynu(x::Float32) = ccall(($jynuf, libm), Float32, (Float32,), x)
        end
    end
    @eval begin
        $bjynu(x::Real) = $bjynu(float(x))
        $bjynu(x::Complex) = $(Symbol(string("bessel", jy)))($nu, x)
        @vectorize_1arg Number $bjynu
    end
end

const cy1 = [0.]
const cy2 = [0.]
const wrk1 = [0.]
const wrk2 = [0.]

function _besselh(nu::Float64, k::Int32, z::Complex128, kode::Int32)
    nz = ierr = Int32(0)
    (nz,ierr) = Amos.ZBESH(real(z), imag(z), nu, kode, k, Int32(1), cy1, cy2, nz, ierr)
    if ierr == 0 || ierr == 3
        return complex(cy1[1],cy2[1])
    else
        throw(AmosException(ierr))
    end
end

function _besseli(nu::Float64, z::Complex128, kode::Int32)
    nz = ierr = Int32(0)
    (nz,ierr) = Amos.ZBESI(real(z), imag(z), nu, kode, Int32(1), cy1, cy2, nz, ierr)
    if ierr == 0 || ierr == 3
        return complex(cy1[1],cy2[1])
    else
        throw(AmosException(ierr))
    end
end

function _besselj(nu::Float64, z::Complex128, kode::Int32)
    nz = ierr = Int32(0)
    (nz,ierr) = Amos.ZBESJ(real(z), imag(z), nu, kode, Int32(1), cy1, cy2, nz, ierr)
    if ierr == 0 || ierr == 3
        return complex(cy1[1],cy2[1])
    else
        throw(AmosException(ierr))
    end

end

function _besselk(nu::Float64, z::Complex128, kode::Int32)
    nz = ierr = Int32(0)
    (nz,ierr) = Amos.ZBESK(real(z), imag(z), nu, kode, Int32(1), cy1, cy2, nz, ierr)
    if ierr == 0 || ierr == 3
        return complex(cy1[1],cy2[1])
    else
        throw(AmosException(ierr))
    end
end

function _bessely(nu::Float64, z::Complex128, kode::Int32)
    nz = ierr = Int32(0)
    (nz,ierr) = Amos.ZBESY(real(z), imag(z), nu, kode, Int32(1), cy1, cy2, nz, wrk1, wrk2, ierr)
    if ierr == 0 || ierr == 3
        return complex(cy1[1],cy2[1])
    else
        throw(AmosException(ierr))
    end
end

function besselh(nu::Float64, k::Integer, z::Complex128)
    if nu < 0
        s = (k == 1) ? 1 : -1
        return _besselh(-nu,Int32(k),z,Int32(1)) * complex(cospi(nu),-s*sinpi(nu))
    end
    return _besselh(nu,Int32(k),z,Int32(1))
end

function besselhx(nu::Float64, k::Integer, z::Complex128)
    if nu < 0
        s = (k == 1) ? 1 : -1
        return _besselh(-nu,Int32(k),z,Int32(2)) * complex(cospi(nu),-s*sinpi(nu))
    end
    return _besselh(nu,Int32(k),z,Int32(2))
end

function besseli(nu::Float64, z::Complex128)
    if nu < 0
        return _besseli(-nu,z,Int32(1)) - 2_besselk(-nu,z,Int32(1))*sinpi(nu)/pi
    else
        return _besseli(nu,z,Int32(1))
    end
end

function besselix(nu::Float64, z::Complex128)
    if nu < 0
        return _besseli(-nu,z,Int32(2)) - 2_besselk(-nu,z,Int32(2))*exp(-abs(real(z))-z)*sinpi(nu)/pi
    else
        return _besseli(nu,z,Int32(2))
    end
end

function besselj(nu::Float64, z::Complex128)
    if nu < 0
        return _besselj(-nu,z,Int32(1))*cospi(nu) + _bessely(-nu,z,Int32(1))*sinpi(nu)
    else
        return _besselj(nu,z,Int32(1))
    end
end

besselj(nu::Integer, x::AbstractFloat) = typemin(Int32) <= nu <= typemax(Int32) ?
    oftype(x, ccall((:jn, libm), Float64, (Cint, Float64), nu, x)) :
    besselj(Float64(nu), x)

besselj(nu::Integer, x::Float32) = typemin(Int32) <= nu <= typemax(Int32) ?
    ccall((:jnf, libm), Float32, (Cint, Float32), nu, x) :
    besselj(Float64(nu), x)

function besseljx(nu::Float64, z::Complex128)
    if nu < 0
        return _besselj(-nu,z,Int32(2))*cospi(nu) + _bessely(-nu,z,Int32(2))*sinpi(nu)
    else
        return _besselj(nu,z,Int32(2))
    end
end

besselk(nu::Float64, z::Complex128) = _besselk(abs(nu), z, Int32(1))

besselkx(nu::Float64, z::Complex128) = _besselk(abs(nu), z, Int32(2))

function bessely(nu::Float64, z::Complex128)
    if nu < 0
        return _bessely(-nu,z,Int32(1))*cospi(nu) - _besselj(-nu,z,Int32(1))*sinpi(nu)
    else
        return _bessely(nu,z,Int32(1))
    end
end

function besselyx(nu::Float64, z::Complex128)
    if nu < 0
        return _bessely(-nu,z,Int32(2))*cospi(nu) - _besselj(-nu,z,Int32(2))*sinpi(nu)
    else
        return _bessely(nu,z,Int32(2))
    end
end


function besseli(nu::Real, x::AbstractFloat)
    if x < 0 && !isinteger(nu)
        throw(DomainError())
    end
    oftype(x, real(besseli(Float64(nu), Complex128(x))))
end

function besselix(nu::Real, x::AbstractFloat)
    if x < 0 && !isinteger(nu)
        throw(DomainError())
    end
    oftype(x, real(besselix(Float64(nu), Complex128(x))))
end

function besselj(nu::AbstractFloat, x::AbstractFloat)
    if isinteger(nu)
        if typemin(Int32) <= nu <= typemax(Int32)
            return besselj(int(nu), x)
        end
    elseif x < 0
        throw(DomainError())
    end
    oftype(x, real(besselj(Float64(nu), Complex128(x))))
end

function besseljx(nu::Real, x::AbstractFloat)
    if x < 0 && !isinteger(nu)
        throw(DomainError())
    end
    oftype(x, real(besseljx(Float64(nu), Complex128(x))))
end

function besselk(nu::Real, x::AbstractFloat)
    if x < 0
        throw(DomainError())
    end
    if x == 0
        return oftype(x, Inf)
    end
    oftype(x, real(besselk(Float64(nu), Complex128(x))))
end

function besselkx(nu::Real, x::AbstractFloat)
    if x < 0
        throw(DomainError())
    end
    if x == 0
        return oftype(x, Inf)
    end
    oftype(x, real(besselkx(Float64(nu), Complex128(x))))
end

function bessely(nu::Real, x::AbstractFloat)
    if x < 0
        throw(DomainError())
    end
    if isinteger(nu) && typemin(Int32) <= nu <= typemax(Int32)
        return bessely(int(nu), x)
    end
    oftype(x, real(bessely(Float64(nu), Complex128(x))))
end
function bessely(nu::Integer, x::AbstractFloat)
    if x < 0
        throw(DomainError())
    end
    return oftype(x, ccall((:yn, libm), Float64, (Cint, Float64), nu, x))
end
function bessely(nu::Integer, x::Float32)
    if x < 0
        throw(DomainError())
    end
    return ccall((:ynf, libm), Float32, (Cint, Float32), nu, x)
end

function besselyx(nu::Real, x::AbstractFloat)
    if x < 0
        throw(DomainError())
    end
    oftype(x, real(besselyx(Float64(nu), Complex128(x))))
end

for f in ("i", "ix", "j", "jx", "k", "kx", "y", "yx")
    bfn = Symbol(string("bessel", f))
    @eval begin
        $bfn(nu::Real, z::Complex64) = Complex64($bfn(Float64(nu), Complex128(z)))
        $bfn(nu::Real, z::Complex) = $bfn(Float64(nu), Complex128(z))
        $bfn(nu::Real, x::Integer) = $bfn(nu, Float64(x))
        @vectorize_2arg Number $bfn
    end
end

for bfn in (:besselh, :besselhx)
    @eval begin
        $bfn(nu, z) = $bfn(nu, 1, z)
        $bfn(nu::Real, k::Integer, x::Real) = $bfn(nu, k, float(x))
        $bfn(nu::Real, k::Integer, x::AbstractFloat) = $bfn(float(nu), k, complex(x))

        function $bfn(nu::Real, k::Integer, z::Complex)
            Tf = promote_type(float(typeof(nu)),float(typeof(real(z))))
            $bfn(Tf(nu), k, Complex{Tf}(z))
        end

        $bfn{T<:AbstractFloat}(nu::T, k::Integer, z::Complex{T}) = throw(MethodError($bfn,(nu,k,z)))
        $bfn(nu::Float32, k::Integer, x::Complex64) = Complex64($bfn(Float64(nu), k, Complex128(x)))
        @vectorize_2arg Number $bfn
    end
end

## Hankel functions

hankelh1(nu, z) = besselh(nu, 1, z)
@vectorize_2arg Number hankelh1

hankelh2(nu, z) = besselh(nu, 2, z)
@vectorize_2arg Number hankelh2

hankelh1x(nu, z) = besselhx(nu, 1, z)
@vectorize_2arg Number hankelh1x

hankelh2x(nu, z) = besselhx(nu, 2, z)
@vectorize_2arg Number hankelh2x

end # module
