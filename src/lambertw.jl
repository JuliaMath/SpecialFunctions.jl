import Base: convert

using Compat
import Compat.MathConstants  # For clarity, we use MathConstants.e for Euler's number

#### Lambert W function ####

# Use Halley's root-finding method to find
# x = lambertw(z) with initial point x.
function _lambertw(z::T, x::T, maxits) where T <: Number
    two_t = convert(T,2)
    lastx = x
    lastdiff = zero(T)
    converged::Bool = false
    for i in 1:maxits
        ex = exp(x)
        xexz = x * ex - z
        x1 = x + 1
        x -= xexz / (ex * x1 - (x + two_t) * xexz / (two_t * x1 ) )
        xdiff = abs(lastx - x)
        if xdiff <= 3*eps(abs(lastx)) || lastdiff == xdiff  # second condition catches two-value cycle
            converged = true
            break
        end
        lastx = x
        lastdiff = xdiff
    end
    converged || warn("lambertw with z=", z, " did not converge in ", maxits, " iterations.")
    return x
end

### Real z ###

# Real x, k = 0
# This appears to be inferrable with T=Float64 and T=BigFloat, including if x=Inf.
# The fancy initial condition selection does not seem to help speed, but we leave it for now.
function lambertwk0(x::T, maxits)::T where T<:AbstractFloat
    isnan(x) && return(NaN)
    x == Inf && return Inf # appears to return convert(BigFloat,Inf) for x == BigFloat(Inf)
    one_t = one(T)
    oneoe = -one_t/convert(T,MathConstants.e)  # The branch point
    x == oneoe && return -one_t
    oneoe <= x || throw(DomainError(x))
    itwo_t = 1/convert(T,2)
    if x > one_t
        lx = log(x)
        llx = log(lx)
        x1 = lx - llx - log(one_t - llx/lx) * itwo_t
    else
        x1 = (567//1000) * x
    end
    return _lambertw(x, x1, maxits)
end

# Real x, k = -1
function lambertwkm1(x::T, maxits) where T<:Real
    oneoe = -one(T)/convert(T,MathConstants.e)
    x == oneoe && return -one(T) # W approaches -1 as x -> -1/e from above
    oneoe <= x || throw(DomainError(x))  # branch domain exludes x < -1/e
    x == zero(T) && return -convert(T,Inf) # W decreases w/o bound as x -> 0 from below
    x < zero(T) || throw(DomainError(x))
    return _lambertw(x, log(-x), maxits)
end

"""
    lambertw(z::Complex{T}, k::V=0, maxits=1000) where {T<:Real, V<:Integer}
    lambertw(z::T, k::V=0, maxits=1000) where {T<:Real, V<:Integer}

Compute the `k`th branch of the Lambert W function of `z`. If `z` is real, `k` must be
either `0` or `-1`. For `Real` `z`, the domain of the branch `k = -1` is `[-1/e,0]` and the
domain of the branch `k = 0` is `[-1/e,Inf]`. For `Complex` `z`, and all `k`, the domain is
the complex plane. When using root finding to compute `W`, a value for `W` is returned
with a warning if it has not converged after `maxits` iterations.

```jldoctest
julia> lambertw(-1/e,-1)
-1.0

julia> lambertw(-1/e,0)
-1.0

julia> lambertw(0,0)
0.0

julia> lambertw(0,-1)
-Inf

julia> lambertw(Complex(-10.0,3.0), 4)
-0.9274337508660128 + 26.37693445371142im
```

"""
lambertw(z, k::Integer=0, maxits::Integer=1000) = lambertw_(z, k, maxits)

function lambertw_(x::Real, k, maxits)
    k == 0 && return lambertwk0(x, maxits)
    k == -1 && return lambertwkm1(x, maxits)
    throw(DomainError(k, "lambertw: real x must have branch k == 0 or k == -1"))
end

function lambertw_(x::Union{Integer,Rational}, k, maxits)
    if k == 0
        x == 0 && return float(zero(x))
        x == 1 && return convert(typeof(float(x)), omega) # must be a more efficient way
    end
    return lambertw_(float(x), k, maxits)
end

### Complex z ###

# choose initial value inside correct branch for root finding
function lambertw_(z::Complex{T}, k, maxits) where T<:Real
    one_t = one(T)
    local w::Complex{T}
    pointseven = 7//10
    if abs(z) <= one_t/convert(T,MathConstants.e)
        if z == 0
            k == 0 && return z
            return complex(-convert(T,Inf),zero(T))
        end
        if k == 0
            w = z
        elseif k == -1 && imag(z) == 0 && real(z) < 0
            w = complex(log(-real(z)),1//10^7) # need offset for z ≈ -1/e.
        else
            w = log(z)
            k != 0 ? w += complex(0,k * 2 * pi) : nothing
        end
    elseif k == 0 && imag(z) <= pointseven && abs(z) <= pointseven
        w = abs(z+ 1//2) < 1//10 ? imag(z) > 0 ? complex(pointseven,pointseven) : complex(pointseven,-pointseven) : z
    else
        if real(z) == convert(T,Inf)
            k == 0 && return z
            return z + complex(0,2*k*pi)
        end
        real(z) == -convert(T,Inf) && return -z + complex(0,(2*k+1)*pi)
        w = log(z)
        k != 0 ? w += complex(0, 2*k*pi) : nothing
    end
    return _lambertw(z, w, maxits)
end

lambertw_(z::Complex{T}, k, maxits) where T<:Integer = lambertw_(float(z), k, maxits)
lambertw_(n::Irrational, k, maxits) = lambertw_(float(n), k, maxits)

# lambertw(e + 0im,k) is ok for all k
# Maybe this should return a float. But, this should cause no type instability in any case
function lambertw_(::typeof(MathConstants.e), k, maxits)
    k == 0 && return 1
    throw(DomainError(k))
end

### omega constant ###

const omega_const_ = 0.567143290409783872999968662210355
# The BigFloat `omega_const_bf_` is set via a literal in the function __init__ to prevent a segfault

# maybe compute higher precision. converges very quickly
function omega_const(::Type{BigFloat})
    precision(BigFloat) <= 256 && return omega_const_bf_[]
    myeps = eps(BigFloat)
    oc = omega_const_bf_[]
    for i in 1:100
        nextoc = (1 + oc) / (1 + exp(oc))
        abs(oc - nextoc) <= myeps && break
        oc = nextoc
    end
    return oc
end

"""
    omega
    ω

The constant defined by `ω exp(ω) = 1`.

```jldoctest
julia> ω
ω = 0.5671432904097...

julia> omega
ω = 0.5671432904097...

julia> ω * exp(ω)
1.0

julia> big(omega)
5.67143290409783872999968662210355549753815787186512508135131079223045793086683e-01
```
"""
const ω = Irrational{:ω}()
@doc (@doc ω) omega = ω

# The following two lines may be removed when support for v0.6 is dropped
Base.convert(::Type{AbstractFloat}, o::Irrational{:ω}) = Float64(o)
Base.convert(::Type{Float16}, o::Irrational{:ω}) = Float16(o)
Base.convert(::Type{T}, o::Irrational{:ω}) where T <:Number = T(o)

Base.Float64(::Irrational{:ω}) = omega_const_  # FIXME: This is very slow. Why ?
Base.Float32(::Irrational{:ω}) = Float32(omega_const_)
Base.Float16(::Irrational{:ω}) = Float16(omega_const_)
Base.BigFloat(o::Irrational{:ω}) = omega_const(BigFloat)

### Expansion about branch point x = -1/e  ###

# Refer to the paper "On the Lambert W function". In (4.22)
# coefficients μ₀ through μ₃ are given explicitly. Recursion relations
# (4.23) and (4.24) for all μ are also given. This code implements the
# recursion relations.

# (4.23) and (4.24) give zero based coefficients.
cset(a,i,v) = a[i+1] = v
cget(a,i) = a[i+1]

# (4.24)
function compa(k,m,a)
    sum0 = zero(eltype(m))
    for j in 2:k-1
        sum0 += cget(m,j) * cget(m,k+1-j)
    end
    cset(a,k,sum0)
    return sum0
end

# (4.23)
function compm(k,m,a)
    kt = convert(eltype(m),k)
    mk = (kt-1)/(kt+1) *(cget(m,k-2)/2 + cget(a,k-2)/4) -
        cget(a,k)/2 - cget(m,k-1)/(kt+1)
    cset(m,k,mk)
    return mk
end

# We plug the known value μ₂ == -1//3 for (4.22) into (4.23) and
# solve for α₂. We get α₂ = 0.
# compute array of coefficients μ in (4.22).
# m[1] is μ₀
function lamwcoeff(T::DataType, n::Int)
    # a = @compat Array{T}(undef,n)
    # m = @compat Array{T}(undef,n)
    a = zeros(T,n) # We don't need initialization, but Compat is a huge PITA.
    m = zeros(T,n)
    cset(a,0,2)  # α₀ literal in paper
    cset(a,1,-1) # α₁ literal in paper
    cset(a,2,0)  # α₂ get this by solving (4.23) for alpha_2 with values printed in paper
    cset(m,0,-1) # μ₀ literal in paper
    cset(m,1,1)  # μ₁ literal in paper
    cset(m,2,-1//3) # μ₂ literal in paper, but only in (4.22)
    for i in 3:n-1  # coeffs are zero indexed
        compa(i,m,a)
        compm(i,m,a)
    end
    return m
end

const LAMWMU_FLOAT64 = lamwcoeff(Float64,500)

# Base.Math.@horner requires literal coefficients
# But, we have an array `p` of computed coefficients
function horner(x, p::AbstractArray, n)
    n += 1
    ex = p[n]
    for i = n-1:-1:2
        ex = :(muladd(t, $ex, $(p[i])))
    end
    ex = :( t * $ex)
    return Expr(:block, :(t = $x), ex)
end

function mkwser(name, n)
    iex = horner(:x,LAMWMU_FLOAT64,n)
    return :(function ($name)(x) $iex  end)
end

eval(mkwser(:wser3, 3))
eval(mkwser(:wser5, 5))
eval(mkwser(:wser7, 7))
eval(mkwser(:wser12, 12))
eval(mkwser(:wser19, 19))
eval(mkwser(:wser26, 26))
eval(mkwser(:wser32, 32))
eval(mkwser(:wser50, 50))
eval(mkwser(:wser100, 100))
eval(mkwser(:wser290, 290))

# Converges to Float64 precision
# We could get finer tuning by separating k=0,-1 branches.
function wser(p,x)
    x < 4e-11 && return wser3(p)
    x < 1e-5 && return wser7(p)
    x < 1e-3 && return wser12(p)
    x < 1e-2 && return wser19(p)
    x < 3e-2 && return wser26(p)
    x < 5e-2 && return wser32(p)
    x < 1e-1 && return wser50(p)
    x < 1.9e-1 && return wser100(p)
    x > 1/MathConstants.e && throw(DomainError(x))  # radius of convergence
    return wser290(p)  # good for x approx .32
end

# These may need tuning.
function wser(p::Complex{T},z) where T<:Real
    x = abs(z)
    x < 4e-11 && return wser3(p)
    x < 1e-5 && return wser7(p)
    x < 1e-3 && return wser12(p)
    x < 1e-2 && return wser19(p)
    x < 3e-2 && return wser26(p)
    x < 5e-2 && return wser32(p)
    x < 1e-1 && return wser50(p)
    x < 1.9e-1 && return wser100(p)
    x > 1/MathConstants.e && throw(DomainError(x))  # radius of convergence
    return wser290(p)
end

@inline function _lambertw0(x) # 1 + W(-1/e + x)  , k = 0
    ps = 2*MathConstants.e*x;
    p = sqrt(ps)
    return  wser(p,x)
end

@inline function _lambertwm1(x) # 1 + W(-1/e + x)  , k = -1
    ps = 2*MathConstants.e*x;
    p = -sqrt(ps)
    return wser(p,x)
end

"""
    lambertwbp(z,k=0)

Compute accurate value of `1 + W(-1/e + z)`, for `abs(z)` in `[0,1/e]` for `k` either `0` or `-1`.
The result is accurate to Float64 precision for abs(z) < 0.32.
If `k=-1` and `imag(z) < 0`, the value on the branch `k=1` is returned.

```jldoctest
julia> lambertw(-1/e + 1e-18, -1)
-1.0

julia> lambertwbp(1e-18, -1)
-2.331643983409312e-9

# Same result, but 1000 times slower
julia> convert(Float64,(lambertw(-BigFloat(1)/e + BigFloat(10)^(-18),-1) + 1))
-2.331643983409312e-9
```

!!! note
    `lambertwbp` uses a series expansion about the branch point `z=-1/e` to avoid loss of precision.
    The loss of precision in `lambertw` is analogous to the loss of precision
    in computing the `sqrt(1-x)` for `x` close to `1`.
"""
function lambertwbp(x::Number,k::Integer)
    k == 0 && return _lambertw0(x)
    k == -1 && return _lambertwm1(x)
    throw(ArgumentError("expansion about branch point only implemented for k = 0 and -1."))
end

lambertwbp(x::Number) = _lambertw0(x)
