#### Lambert W function ####

"""
    lambertw(z::Complex{T}, k::V=0, maxits=1000) where {T<:Real, V<:Integer}
    lambertw(z::T, k::V=0, maxits=1000) where {T<:Real, V<:Integer}

Compute the `k`th branch of the Lambert W function of `z`. If `z` is real, `k` must be
either `0` or `-1`. For `Real` `z`, the domain of the branch `k = -1` is `[-1/e, 0]` and the
domain of the branch `k = 0` is `[-1/e, Inf]`. For `Complex` `z`, and all `k`, the domain is
the complex plane.

```jldoctest
julia> lambertw(-1/e, -1)
-1.0

julia> lambertw(-1/e, 0)
-1.0

julia> lambertw(0, 0)
0.0

julia> lambertw(0, -1)
-Inf

julia> lambertw(Complex(-10.0, 3.0), 4)
-0.9274337508660128 + 26.37693445371142im
```
"""
lambertw(z, k::Integer=0, maxits::Integer=1000) = _lambertw(float(z), k, maxits)

# lambertw(e + 0im, k) is ok for all k
# Maybe this should return a float. But, this should cause no type instability in any case
function _lambertw(::typeof(MathConstants.e), k, maxits)
    k == 0 && return 1
    throw(DomainError(k))
end
_lambertw(x::Irrational, k, maxits) = _lambertw(float(x), k, maxits)
function _lambertw(x::Union{Integer, Rational}, k, maxits)
    if k == 0
        x == 0 && return float(zero(x))
        x == 1 && return convert(typeof(float(x)), omega) # must be a more efficient way
    end
    return _lambertw(float(x), k, maxits)
end

### Real x

function _lambertw(x::Real, k, maxits)
    k == 0 && return lambertw_branch_zero(x, maxits)
    k == -1 && return lambertw_branch_one(x, maxits)
    throw(DomainError(k, "lambertw: real x must have branch k == 0 or k == -1"))
end

# Real x, k = 0
# There is a magic number here. It could be noted, or possibly removed.
# In particular, the fancy initial condition selection does not seem to help speed.
function lambertw_branch_zero(x::T, maxits) where T<:Real
    isfinite(x) || return x
    one_t = one(T)
    oneoe = -one_t / convert(T, MathConstants.e)  # The branch point
    x == oneoe && return -one_t
    oneoe < x || throw(DomainError(x))
    itwo_t = 1 / convert(T, 2)
    if x > one_t
        lx = log(x)
        llx = log(lx)
        x0 = lx - llx - log(one_t - llx / lx) * itwo_t
    else
        x0 = (567//1000) * x
    end
    return lambertw_root_finding(x, x0, maxits)
end

# Real x, k = -1
function lambertw_branch_one(x::T, maxits) where T<:Real
    oneoe = -one(T) / convert(T, MathConstants.e)
    x == oneoe && return -one(T) # W approaches -1 as x -> -1/e from above
    oneoe < x || throw(DomainError(x))  # branch domain exludes x < -1/e
    x == zero(T) && return -convert(T, Inf) # W decreases w/o bound as x -> 0 from below
    x < zero(T) || throw(DomainError(x))
    return lambertw_root_finding(x, log(-x), maxits)
end

### Complex z

_lambertw(z::Complex{<:Integer}, k, maxits) = _lambertw(float(z), k, maxits)
# choose initial value inside correct branch for root finding
function _lambertw(z::Complex{T}, k, maxits) where T<:Real
    one_t = one(T)
    local w::Complex{T}
    pointseven = 7//10
    if abs(z) <= one_t/convert(T, MathConstants.e)
        if z == 0
            k == 0 && return z
            return complex(-convert(T, Inf), zero(T))
        end
        if k == 0
            w = z
        elseif k == -1 && imag(z) == 0 && real(z) < 0
            w = complex(log(-real(z)), 1//10^7) # need offset for z ≈ -1/e.
        else
            w = log(z)
            k != 0 ? w += complex(0, k * 2 * pi) : nothing
        end
    elseif k == 0 && imag(z) <= pointseven && abs(z) <= pointseven
        w = abs(z+ 1//2) < 1//10 ? imag(z) > 0 ? complex(pointseven, pointseven) : complex(pointseven, -pointseven) : z
    else
        if real(z) == convert(T, Inf)
            k == 0 && return z
            return z + complex(0, 2*k*pi)
        end
        real(z) == -convert(T, Inf) && return -z + complex(0, (2*k+1)*pi)
        w = log(z)
        k != 0 ? w += complex(0, 2*k*pi) : nothing
    end
    return lambertw_root_finding(z, w, maxits)
end

### root finding, iterative solution

# Use Halley's root-finding method to find
# x = lambertw(z) with initial point x0.
function lambertw_root_finding(z::T, x0::T, maxits) where T <: Number
    two_t = convert(T, 2)
    x = x0
    lastx = x
    lastdiff = zero(real(T))
    converged = false
    for i in 1:maxits
        ex = exp(x)
        xexz = x * ex - z
        x1 = x + 1
        x -= xexz / (ex * x1 - (x + two_t) * xexz / (two_t * x1 ))
        xdiff = abs(lastx - x)
        if xdiff <= 3 * eps(lastdiff) || lastdiff == xdiff  # second condition catches two-value cycle
            converged = true
            break
        end
        lastx = x
        lastdiff = xdiff
    end
    converged || @warn("lambertw with z=", z, " did not converge in ", maxits, " iterations.")
    return x
end

### omega constant

const _omega_const = 0.567143290409783872999968662210355

# The BigFloat `omega_const_bf_` is set via a literal in the function __init__ to prevent a segfault

# compute omega constant via root finding
# We could compute higher precision. This converges very quickly.
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

# Example
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

Base.Float64(::Irrational{:ω}) = _omega_const
Base.Float32(::Irrational{:ω}) = Float32(_omega_const)
Base.Float16(::Irrational{:ω}) = Float16(_omega_const)
Base.BigFloat(::Irrational{:ω}) = omega_const(BigFloat)

### Expansion about branch point x = -1/e

# Refer to the paper "On the Lambert W function". In (4.22)
# coefficients μ₀ through μ₃ are given explicitly. Recursion relations
# (4.23) and (4.24) for all μ are also given. This code implements the
# recursion relations.

# We plug the known value μ₂ == -1//3 for (4.22) into (4.23) and
# solve for α₂. We get α₂ = 0.
# Compute array of coefficients μ in (4.22).
# m[1] is μ₀
function compute_branch_point_coeffs(T::DataType, n::Int)
    a = Vector{T}(undef, n)
    m = Vector{T}(undef, n)

    a[1] = 2  # α₀ literal in paper
    a[2] = -1 # α₁ literal in paper
    a[3] = 0  # α₂ get this by solving (4.23) for alpha_2 with values printed in paper
    m[1] = -1 # μ₀ literal in paper
    m[2] = 1  # μ₁ literal in paper
    m[3] = -1//3 # μ₂ literal in paper, but only in (4.22)

    for i in 4:n
        # (4.24)
        msum = zero(T)
        @inbounds for j in 2:(i - 2)
            msum += m[j + 1] * m[i + 1 - j]
        end
        a[i] = msum

        # (4.23)
        it = convert(T, i)
        m[i] = (it - 2) / it *(m[i - 2] / 2 + a[i - 2] / 4) -
               a[i] / 2 - m[i - 1] / it
    end
    return m
end

const BRANCH_POINT_COEFFS_FLOAT64 = compute_branch_point_coeffs(Float64, 500)

# Base.Math.@horner requires literal coefficients
# It cannot be used here because we have an array of computed coefficients
function horner(x, coeffs::AbstractArray, n)
    n += 1
    ex = coeffs[n]
    for i = (n - 1):-1:2
        ex = :(muladd(t, $ex, $(coeffs[i])))
    end
    ex = :( t * $ex)
    return Expr(:block, :(t = $x), ex)
end

# write functions that evaluate the branch point series
# with `num_terms` number of terms.
for (func_name, num_terms) in (
    (:wser3, 3), (:wser5, 5), (:wser7, 7), (:wser12, 12),
    (:wser19, 19), (:wser26, 26), (:wser32, 32),
    (:wser50, 50), (:wser100, 100), (:wser290, 290))
    iex = horner(:x, BRANCH_POINT_COEFFS_FLOAT64, num_terms)
    @eval function ($func_name)(x) $iex end
end

# Converges to Float64 precision
# We could get finer tuning by separating k=0, -1 branches.
# Why is wser5 omitted ?
# p is the argument to the series which is computed
# from x before calling `branch_point_series`.
function branch_point_series(p, x)
    x < 4e-11 && return wser3(p)
    x < 1e-5 && return wser7(p)
    x < 1e-3 && return wser12(p)
    x < 1e-2 && return wser19(p)
    x < 3e-2 && return wser26(p)
    x < 5e-2 && return wser32(p)
    x < 1e-1 && return wser50(p)
    x < 1.9e-1 && return wser100(p)
    x > 1 / MathConstants.e && throw(DomainError(x))  # radius of convergence
    return wser290(p)  # good for x approx .32
end

# These may need tuning.
function branch_point_series(p::Complex{T}, z) where T<:Real
    x = abs(z)
    x < 4e-11 && return wser3(p)
    x < 1e-5 && return wser7(p)
    x < 1e-3 && return wser12(p)
    x < 1e-2 && return wser19(p)
    x < 3e-2 && return wser26(p)
    x < 5e-2 && return wser32(p)
    x < 1e-1 && return wser50(p)
    x < 1.9e-1 && return wser100(p)
    x > 1 / MathConstants.e && throw(DomainError(x))  # radius of convergence
    return wser290(p)
end

function _lambertw0(x) # 1 + W(-1/e + x)  , k = 0
    ps = 2 * MathConstants.e * x
    series_arg = sqrt(ps)
    branch_point_series(series_arg, x)
end

function _lambertwm1(x) # 1 + W(-1/e + x)  , k = -1
    ps = 2 * MathConstants.e * x
    series_arg = -sqrt(ps)
    branch_point_series(series_arg, x)
end

"""
    lambertwbp(z, k=0)

Compute accurate value of `1 + W(-1/e + z)`, for `abs(z)` in `[0, 1/e]` for `k` either `0` or `-1`.
This function is faster and more accurate near the branch point `-1/e` between `k=0` and `k=1`.
The result is accurate to Float64 precision for abs(z) < 0.32.
If `k=-1` and `imag(z) < 0`, the value on the branch `k=1` is returned.

# Example
```jldoctest
julia> lambertw(-1/e + 1e-18, -1)
-1.0

julia> lambertwbp(1e-18, -1)
-2.331643983409312e-9

# Same result, but 1000 times slower
julia> convert(Float64, (lambertw(-BigFloat(1)/e + BigFloat(10)^(-18), -1) + 1))
-2.331643983409312e-9
```

!!! note
    `lambertwbp` uses a series expansion about the branch point `z=-1/e` to avoid loss of precision.
    The loss of precision in `lambertw` is analogous to the loss of precision
    in computing the `sqrt(1-x)` for `x` close to `1`.
"""
function lambertwbp(x::Number, k::Integer)
    k == 0 && return _lambertw0(x)
    k == -1 && return _lambertwm1(x)
    throw(ArgumentError("expansion about branch point only implemented for k = 0 and -1."))
end

lambertwbp(x::Number) = _lambertw0(x)
