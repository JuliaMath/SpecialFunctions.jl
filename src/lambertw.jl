#### Lambert W function ####

"""
    lambertw(z::Number, k::Integer=0; [maxiter=1000])

Compute the `k`th branch of the [Lambert W function](https://en.wikipedia.org/wiki/Lambert_W_function) of `z`. If `z` is real, `k` must be
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
lambertw(z::Number, k::Integer=0; maxiter::Integer=1000) = _lambertw(z, k, maxiter)

# lambertw(e + 0im, k) is ok for all k
# Maybe this should return a float. But, this should cause no type instability in any case
function _lambertw(::typeof(MathConstants.e), k::Integer, maxits::Integer)
    k == 0 && return 1
    throw(DomainError(k))
end
_lambertw(x::Irrational, k::Integer, maxits::Integer) = _lambertw(float(x), k, maxits)
function _lambertw(x::Union{Integer, Rational}, k::Integer, maxits::Integer)
    if k == 0
        x == 0 && return float(zero(x))
        x == 1 && return convert(typeof(float(x)), LambertW.Omega) # must be a more efficient way
    end
    return _lambertw(float(x), k, maxits)
end

### Real x

_lambertw(x::Real, k::Integer, maxits::Integer) = _lambertw(x, Val(Int(k)), maxits)

# Real x, k = 0
# There is a magic number here. It could be noted, or possibly removed.
# In particular, the fancy initial condition selection does not seem to help speed.
function _lambertw(x::T, ::Val{0}, maxits::Integer) where T<:Real
    isfinite(x) || return x
    one_t = one(T)
    oneoe = -T(inve)  # The branch point
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
function _lambertw(x::T, ::Val{-1}, maxits::Integer) where T<:Real
    oneoe = -T(inve)
    x == oneoe && return -one(T) # W approaches -1 as x -> -1/e from above
    oneoe < x || throw(DomainError(x))  # branch domain exludes x < -1/e
    x == zero(T) && return -convert(T, Inf) # W decreases w/o bound as x -> 0 from below
    x < zero(T) || throw(DomainError(x))
    return lambertw_root_finding(x, log(-x), maxits)
end

_lambertw(x::Real, k::Val, maxits::Integer) =
    throw(DomainError(x, "lambertw: for branch k=$k not defined, real x must have branch k == 0 or k == -1"))

### Complex z

_lambertw(z::Complex{<:Integer}, k::Integer, maxits::Integer) = _lambertw(float(z), k, maxits)
# choose initial value inside correct branch for root finding
function _lambertw(z::Complex{T}, k::Integer, maxits::Integer) where T<:Real
    local w::Complex{T}
    pointseven = 7//10
    if abs(z) <= T(inve)
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
        w = abs(z+ 1//2) < 1//10 ? imag(z) > 0 ?
                complex(pointseven, pointseven) :
                complex(pointseven, -pointseven) : z
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
function lambertw_root_finding(z::T, x0::T, maxits::Integer) where T <: Number
    x = x0
    lastx = x
    lastdiff = zero(real(T))
    converged = false
    for _ in 1:maxits
        ex = exp(x)
        xexz = x * ex - z
        x1 = x + 1
        x -= xexz / (ex * x1 - (x + 2) * xexz / (2 * x1))
        xdiff = abs(lastx - x)
        if xdiff <= 3 * eps(lastdiff) || lastdiff == xdiff  # second condition catches two-value cycle
            converged = true
            break
        end
        lastx = x
        lastdiff = xdiff
    end
    converged || @warn "lambertw($z) did not converge in $maxits iterations."
    return x
end

### Lambert's Omega constant

# compute BigFloat Omega constant at arbitrary precision
function compute_lambertw_Omega()
    oc = BigFloat("0.5671432904097838729999686622103555497538157871865125081351310792230457930866845666932194")
    precision(oc) <= 256 && return oc
    # iteratively improve the precision
    # see https://en.wikipedia.org/wiki/Omega_constant#Computation
    myeps = eps(BigFloat)
    for _ in 1:1000
        nextoc = (1 + oc) / (1 + exp(oc))
        abs(oc - nextoc) <= myeps && return oc
        oc = nextoc
    end
    @warn "Omega precision is less than current BigFloat precision ($(precision(BigFloat)))"
    return oc
end

# "private" declaration of Omega constant
Base.@irrational lambertw_Omega 0.567143290409783872999968662210355 compute_lambertw_Omega()

module LambertW

"""
Lambert's Omega (*Ω*) constant.

Lambert's *Ω* is the solution to *W(Ω) = 1* equation,
where *W(t) = t exp(t)* is the
[Lambert's *W* function](https://en.wikipedia.org/wiki/Lambert_W_function).

# See also
  * https://en.wikipedia.org/wiki/Omega_constant
  * [`lambertw()`][@ref SpecialFunctions.lambertw]
"""
const Ω = Irrational{:lambertw_Omega}()
const Omega = Ω # ASCII alias

end

### Expansion about branch point x = -1/e

# Refer to the paper "On the Lambert W function". In (4.22)
# coefficients μ₀ through μ₃ are given explicitly. Recursion relations
# (4.23) and (4.24) for all μ are also given. This code implements the
# recursion relations.

# We plug the known value μ₂ == -1//3 for (4.22) into (4.23) and
# solve for α₂. We get α₂ = 0.
# Compute array of coefficients μ in (4.22).
# m[1] is μ₀
function lambertw_coeffs(T::Type{<:Number}, n::Integer)
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

const LAMBERTW_COEFFS_FLOAT64 = lambertw_coeffs(Float64, 500)

(lambertwbp_evalpoly(x::T, ::Val{N})::T) where {T<:Number, N} =
    # assume that Julia compiler is smart to decide for which N to unroll at compile time
    # note that we skip μ₀=-1
    evalpoly(x, ntuple(i -> LAMBERTW_COEFFS_FLOAT64[i+1], N-1))*x

# how many coefficients of the series to use
# to converge to Float64 precision for given x
# We could get finer tuning by separating k=0, -1 branches.
function lambertwbp_series_length(x::Real)
    x < 4e-11 && return 3
    # Why N = 5 is omitted?
    x < 1e-5 && return 7
    x < 1e-3 && return 12
    x < 1e-2 && return 19
    x < 3e-2 && return 26
    x < 5e-2 && return 32
    x < 1e-1 && return 50
    x < 1.9e-1 && return 100
    x > typeof(x)(inve) && throw(DomainError(x))  # radius of convergence
    return 290 # good for x approx .32
end

# These may need tuning.
lambertwbp_series_length(z::Complex) = lambertwbp_series_length(abs(z))

# p is the argument to the series which is computed from x,
# see `_lambertwbp()`.
lambertwbp_series(p::Number, x::Number) =
    lambertwbp_evalpoly(p, Val{lambertwbp_series_length(x)}())

_lambertwbp(x::Number, ::Val{0}) =
    lambertwbp_series(sqrt(2 * MathConstants.e * x), x)

_lambertwbp(x::Number, ::Val{-1}) =
    lambertwbp_series(-sqrt(2 * MathConstants.e * x), x)

_lambertwbp(_::Number, k::Val) =
    throw(ArgumentError("lambertw() expansion about branch point for k=$k not implemented (only implemented for 0 and -1)."))

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
lambertwbp(x::Number, k::Integer=0) = _lambertwbp(x, Val(Int(k)))
