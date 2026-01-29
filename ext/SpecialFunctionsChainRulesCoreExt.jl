module SpecialFunctionsChainRulesCoreExt

using SpecialFunctions
using ChainRulesCore: ChainRulesCore

using SpecialFunctions: sqrtπ, invπ

const BESSEL_ORDER_INFO = """
derivatives of Bessel functions with respect to the order are not implemented currently:
https://github.com/JuliaMath/SpecialFunctions.jl/issues/160
"""

const INCOMPLETE_GAMMA_INFO = """
derivatives of the incomplete Gamma functions with respect to parameter `a` are not
implemented currently:
https://github.com/JuliaMath/SpecialFunctions.jl/issues/317
"""

const INCOMPLETE_EXPINT_INFO = """
derivatives of the exponential integral with respect to parameter `ν` are not implemented
currently:
https://github.com/JuliaMath/SpecialFunctions.jl/issues/321
"""

ChainRulesCore.@scalar_rule(airyai(x), airyaiprime(x))
ChainRulesCore.@scalar_rule(airyaix(x), airyaiprimex(x) + sqrt(x) * Ω)
ChainRulesCore.@scalar_rule(airyaiprime(x), x * airyai(x))
ChainRulesCore.@scalar_rule(airyaiprimex(x), x * airyaix(x) + sqrt(x) * Ω)
ChainRulesCore.@scalar_rule(airybi(x), airybiprime(x))
ChainRulesCore.@scalar_rule(airybiprime(x), x * airybi(x))
ChainRulesCore.@scalar_rule(besselj0(x), -besselj1(x))
ChainRulesCore.@scalar_rule(
    besselj1(x),
    (besselj0(x) - besselj(2, x)) / 2,
)
ChainRulesCore.@scalar_rule(bessely0(x), -bessely1(x))
ChainRulesCore.@scalar_rule(
    bessely1(x),
    (bessely0(x) - bessely(2, x)) / 2,
)
ChainRulesCore.@scalar_rule(dawson(x), 1 - (2 * x * Ω))
ChainRulesCore.@scalar_rule(digamma(x), trigamma(x))

# TODO: use `invsqrtπ` if it is added to IrrationalConstants
ChainRulesCore.@scalar_rule(erf(x), (2 * exp(-x^2)) / sqrtπ)
ChainRulesCore.@scalar_rule(erf(x, y), (- (2 * exp(-x^2)) / sqrtπ, (2 * exp(-y^2)) / sqrtπ))
ChainRulesCore.@scalar_rule(erfc(x), - (2 * exp(-x^2)) / sqrtπ)
ChainRulesCore.@scalar_rule(logerfc(x), - (2 * exp(-x^2 - Ω)) / sqrtπ)
ChainRulesCore.@scalar_rule(erfcinv(x), - (sqrtπ * (exp(Ω^2) / 2)))
ChainRulesCore.@scalar_rule(erfcx(x), 2 * (x * Ω - inv(oftype(Ω, sqrtπ))))
ChainRulesCore.@scalar_rule(logerfcx(x), 2 * (x - exp(-Ω) / sqrtπ))
ChainRulesCore.@scalar_rule(erfi(x), (2 * exp(x^2)) / sqrtπ)
ChainRulesCore.@scalar_rule(erfinv(x), sqrtπ * (exp(Ω^2) / 2))

ChainRulesCore.@scalar_rule(gamma(x), Ω * digamma(x))
ChainRulesCore.@scalar_rule(
    gamma(a, x),
    (
        ChainRulesCore.@not_implemented(INCOMPLETE_GAMMA_INFO),
        - exp(-x) * x^(a - 1),
    ),
)
ChainRulesCore.@scalar_rule(
    gamma_inc(a, x, IND),
    @setup(z = exp(-x) * x^(a - 1) / gamma(a)),
    (
        ChainRulesCore.@not_implemented(INCOMPLETE_GAMMA_INFO),
        z,
        ChainRulesCore.NoTangent(),
    ),
    (
        ChainRulesCore.@not_implemented(INCOMPLETE_GAMMA_INFO),
        -z,
        ChainRulesCore.NoTangent(),
    ),
)
ChainRulesCore.@scalar_rule(
    invdigamma(x),
    inv(trigamma(invdigamma(x))),
)
ChainRulesCore.@scalar_rule(trigamma(x), polygamma(2, x))

# Bessel functions
ChainRulesCore.@scalar_rule(
    besselj(ν, x),
    (
        ChainRulesCore.@not_implemented(BESSEL_ORDER_INFO),
        (besselj(ν - 1, x) - besselj(ν + 1, x)) / 2
    ),
)
ChainRulesCore.@scalar_rule(
    besseli(ν, x),
    (
        ChainRulesCore.@not_implemented(BESSEL_ORDER_INFO),
        (besseli(ν - 1, x) + besseli(ν + 1, x)) / 2,
    ),
)
ChainRulesCore.@scalar_rule(
    bessely(ν, x),
    (
        ChainRulesCore.@not_implemented(BESSEL_ORDER_INFO),
        (bessely(ν - 1, x) - bessely(ν + 1, x)) / 2,
    ),
)
ChainRulesCore.@scalar_rule(
    besselk(ν, x),
    (
        ChainRulesCore.@not_implemented(BESSEL_ORDER_INFO),
        -(besselk(ν - 1, x) + besselk(ν + 1, x)) / 2,
    ),
)
ChainRulesCore.@scalar_rule(
    besselkx(ν, x),
    (
        ChainRulesCore.@not_implemented(BESSEL_ORDER_INFO),
        -(besselkx(ν - 1, x) + besselkx(ν + 1, x)) / 2 + Ω,
    ),
)
ChainRulesCore.@scalar_rule(
    hankelh1(ν, x),
    (
        ChainRulesCore.@not_implemented(BESSEL_ORDER_INFO),
        (hankelh1(ν - 1, x) - hankelh1(ν + 1, x)) / 2,
    ),
)
ChainRulesCore.@scalar_rule(
    hankelh1x(ν, x),
    (
        ChainRulesCore.@not_implemented(BESSEL_ORDER_INFO),
        (hankelh1x(ν - 1, x) - hankelh1x(ν + 1, x)) / 2 - im * Ω,
    ),
)
ChainRulesCore.@scalar_rule(
    hankelh2(ν, x),
    (
        ChainRulesCore.@not_implemented(BESSEL_ORDER_INFO),
        (hankelh2(ν - 1, x) - hankelh2(ν + 1, x)) / 2,
    ),
)
ChainRulesCore.@scalar_rule(
    hankelh2x(ν, x),
    (
        ChainRulesCore.@not_implemented(BESSEL_ORDER_INFO),
        (hankelh2x(ν - 1, x) - hankelh2x(ν + 1, x)) / 2 + im * Ω,
    ),
)

ChainRulesCore.@scalar_rule(
    polygamma(m, x),
    (
        ChainRulesCore.NoTangent(),
        polygamma(m + 1, x),
    ),
)
# todo: setup for common expr
ChainRulesCore.@scalar_rule(
    beta(a, b),
    (Ω*(digamma(a) - digamma(a + b)),
     Ω*(digamma(b) - digamma(a + b)),)
)
ChainRulesCore.@scalar_rule(
    logbeta(a, b),
    (digamma(a) - digamma(a + b),
     digamma(b) - digamma(a + b),)
)

# actually is the absolute value of the logorithm of gamma paired with sign gamma
ChainRulesCore.@scalar_rule(
    logabsgamma(x),
    digamma(x),
    ChainRulesCore.ZeroTangent()
)

ChainRulesCore.@scalar_rule(loggamma(x), digamma(x))
ChainRulesCore.@scalar_rule(
    loggamma(a, x),
    (
        ChainRulesCore.@not_implemented(INCOMPLETE_GAMMA_INFO),
        -exp(- (x + Ω)) * x^(a - 1),
    )
)

# exponential integrals
ChainRulesCore.@scalar_rule(expint(z), - exp(-z) / z)
ChainRulesCore.@scalar_rule(
    expint(ν, z),
    (
        ChainRulesCore.@not_implemented(INCOMPLETE_EXPINT_INFO),
        - expint(ν - 1, z),
    )
)
ChainRulesCore.@scalar_rule(expintx(z), Ω - inv(z))
ChainRulesCore.@scalar_rule(
    expintx(ν, z),
    (
        ChainRulesCore.@not_implemented(INCOMPLETE_EXPINT_INFO),
        Ω - expintx(ν - 1, z),
    )
)
ChainRulesCore.@scalar_rule(expinti(x), exp(x) / x)
ChainRulesCore.@scalar_rule(sinint(x), sinc(invπ * x))
ChainRulesCore.@scalar_rule(cosint(x), cos(x) / x)

# elliptic integrals
ChainRulesCore.@scalar_rule(
    ellipk(m),
    iszero(m) ? oftype(Ω, π) / 8 : (ellipe(m) / (1 - m) - Ω) / (2 * m),
)
ChainRulesCore.@scalar_rule(
    ellipe(m),
    iszero(m) ? -oftype(Ω, π) / 8 : (Ω - ellipk(m)) / (2 * m),
)

# non-holomorphic functions
function ChainRulesCore.frule((_, Δν, Δx), ::typeof(besselix), ν::Number, x::Number)
    # primal
    Ω = besselix(ν, x)

    # derivative
    ∂Ω_∂ν = ChainRulesCore.@not_implemented(BESSEL_ORDER_INFO)
    a = (besselix(ν - 1, x) + besselix(ν + 1, x)) / 2
    ΔΩ = if Δx isa Real
        muladd(muladd(-sign(real(x)), Ω, a), Δx, ∂Ω_∂ν * Δν)
    else
        muladd(a, Δx, muladd(-sign(real(x)) * real(Δx), Ω, ∂Ω_∂ν * Δν))
    end

    return Ω, ΔΩ
end
function ChainRulesCore.rrule(::typeof(besselix), ν::Number, x::Number)
    Ω = besselix(ν, x)
    project_x = ChainRulesCore.ProjectTo(x)
    function besselix_pullback(ΔΩ)
        ν̄ = ChainRulesCore.@not_implemented(BESSEL_ORDER_INFO)
        a = (besselix(ν - 1, x) + besselix(ν + 1, x)) / 2
        x̄ = project_x(muladd(conj(a), ΔΩ, - sign(real(x)) * real(conj(Ω) * ΔΩ)))
        return ChainRulesCore.NoTangent(), ν̄, x̄
    end
    return Ω, besselix_pullback
end

function ChainRulesCore.frule((_, Δν, Δx), ::typeof(besseljx), ν::Number, x::Number)
    # primal
    Ω = besseljx(ν, x)

    # derivative
    ∂Ω_∂ν = ChainRulesCore.@not_implemented(BESSEL_ORDER_INFO)
    a = (besseljx(ν - 1, x) - besseljx(ν + 1, x)) / 2
    ΔΩ = if Δx isa Real
        muladd(a, Δx, ∂Ω_∂ν * Δν)
    else
        muladd(a, Δx, muladd(-sign(imag(x)) * imag(Δx), Ω, ∂Ω_∂ν * Δν))
    end

    return Ω, ΔΩ
end
function ChainRulesCore.rrule(::typeof(besseljx), ν::Number, x::Number)
    Ω = besseljx(ν, x)
    project_x = ChainRulesCore.ProjectTo(x)
    function besseljx_pullback(ΔΩ)
        ν̄ = ChainRulesCore.@not_implemented(BESSEL_ORDER_INFO)
        a = (besseljx(ν - 1, x) - besseljx(ν + 1, x)) / 2
        x̄ = if x isa Real
            project_x(a * ΔΩ)
        else
            project_x(muladd(conj(a), ΔΩ, - sign(imag(x)) * real(conj(Ω) * ΔΩ) * im))
        end
        return ChainRulesCore.NoTangent(), ν̄, x̄
    end
    return Ω, besseljx_pullback
end

function ChainRulesCore.frule((_, Δν, Δx), ::typeof(besselyx), ν::Number, x::Number)
    # primal
    Ω = besselyx(ν, x)

    # derivative
    ∂Ω_∂ν = ChainRulesCore.@not_implemented(BESSEL_ORDER_INFO)
    a = (besselyx(ν - 1, x) - besselyx(ν + 1, x)) / 2
    ΔΩ = if Δx isa Real
        muladd(a, Δx, ∂Ω_∂ν * Δν)
    else
        muladd(a, Δx, muladd(-sign(imag(x)) * imag(Δx), Ω, ∂Ω_∂ν * Δν))
    end

    return Ω, ΔΩ
end
function ChainRulesCore.rrule(::typeof(besselyx), ν::Number, x::Number)
    Ω = besselyx(ν, x)
    project_x = ChainRulesCore.ProjectTo(x)
    function besselyx_pullback(ΔΩ)
        ν̄ = ChainRulesCore.@not_implemented(BESSEL_ORDER_INFO)
        a = (besselyx(ν - 1, x) - besselyx(ν + 1, x)) / 2
        x̄ = if x isa Real
            project_x(a * ΔΩ)
        else
            project_x(muladd(conj(a), ΔΩ, - sign(imag(x)) * real(conj(Ω) * ΔΩ) * im))
        end
        return ChainRulesCore.NoTangent(), ν̄, x̄
    end
    return Ω, besselyx_pullback
end



## Incomplete beta derivatives via Boik & Robinson-Cox
#
# Reference
#   R. J. Boik and J. F. Robinson-Cox (1999).
#   "Derivatives of the incomplete beta function."
#   Journal of Statistical Software, 3(1).
#   URL: https://www.jstatsoft.org/article/view/v003i01
#
# The following implementation computes the regularized incomplete beta
# I_x(a,b) together with its partial derivatives with respect to a, b, and x
# using a continued-fraction representation of ₂F₁ and differentiating through it.
# This is an independent implementation adapted from https://github.com/arzwa/IncBetaDer.jl.

# Generic-typed helpers used by the continued-fraction evaluation of I_x(a,b)
# and its partial derivatives. These implement the scalar prefactor K(x;p,q),
# the auxiliary variable f, the continued-fraction coefficients a_n, b_n, and
# their partial derivatives w.r.t. p (≡ a) and q (≡ b). See Boik & Robinson-Cox (1999).

@inline function _Kfun(x::T, p::T, q::T) where {T<:AbstractFloat}
    # K(x;p,q) = x^p (1-x)^{q-1} / (p * B(p,q)) computed in log-space for stability
    return exp(p * log(x) + (q - 1) * log1p(-x) - log(p) - logbeta(p, q))
end

@inline function _ffun(x::T, p::T, q::T) where {T<:AbstractFloat}
    # f = q x / (p (1-x)) — convenience variable appearing in CF coefficients
    return q * x / (p * (1 - x))
end

@inline function _a1fun(p::T, q::T, f::T) where {T<:AbstractFloat} 
    # a₁ coefficient of the continued fraction for ₂F₁ representation
    return p * f * (q - 1) / (q * (p + 1))
end

@inline function _anfun(p::T, q::T, f::T, n::Int) where {T<:AbstractFloat}
    # a_n coefficient (n ≥ 1) of the continued fraction for ₂F₁ in terms of p=a, q=b, f.
    # For n=1, falls back to a₁; for n≥2 uses the closed-form product from the Gauss CF.
    n == 1 && return _a1fun(p, q, f)
    r = (p * f / q)^2
    pn = p + n
    p2n = pn + n
    return r * (n - 1) * (pn + q - 2) * (pn - 1) * (q - n) / ((p2n - 3) * (p2n - 2)^2 * (p2n - 1))
end

@inline function _bnfun(p::T, q::T, f::T, n::Int) where {T<:AbstractFloat}
    # b_n coefficient (n ≥ 1) of the continued fraction. Derived for the same CF.
    x = 2 * n * (p * f + 2 * q) * (n + p - 1) + p * q * (p - 2 - p * f)
    y = q * (p + 2*n - 2) * (p + 2*n)
    return x / y
end

@inline function _dK_dp(x::T, p::T, q::T, K::T, ψpq::T, ψp::T) where {T<:AbstractFloat} 
    # ∂K/∂p using digamma identities: d/dp log B(p,q) = ψ(p) - ψ(p+q)
    return K * (log(x) - inv(p) + ψpq - ψp)
end

@inline function _dK_dq(x::T, p::T, q::T, K::T, ψpq::T, ψq::T) where {T<:AbstractFloat} 
    # ∂K/∂q using identical pattern
    K * (log1p(-x) + ψpq - ψq)
end

@inline function _dK_dpdq(x::T, p::T, q::T) where {T<:AbstractFloat}
    # Convenience: compute (∂K/∂p, ∂K/∂q) together with shared ψ(p+q)
    ψ = digamma(p + q)
    Kf = _Kfun(x, p, q)
    dKdp = _dK_dp(x, p, q, Kf, ψ, digamma(p))
    dKdq = _dK_dq(x, p, q, Kf, ψ, digamma(q))
    return dKdp, dKdq
end

@inline function _da1_dp(p::T, q::T, f::T) where {T<:AbstractFloat}
    # ∂a₁/∂p from the closed form of a₁
    return - _a1fun(p, q, f) / (p + 1)
end

@inline function _dan_dp(p::T, q::T, f::T, n::Int) where {T<:AbstractFloat}
    # ∂a_n/∂p via log-derivative: d a_n = a_n * d log a_n; for n=1, uses ∂a₁/∂p
    if n == 1
        return _da1_dp(p, q, f)
    end
    an = _anfun(p, q, f, n)
    dlog = inv(p + q + n - 2) + inv(p + n - 1) - inv(p + 2*n - 3) - 2 * inv(p + 2*n - 2) - inv(p + 2*n - 1)
    return an * dlog
end

@inline function _da1_dq(p::T, q::T, f::T) where {T<:AbstractFloat}
    # ∂a₁/∂q
    return _a1fun(p, q, f) / (q - 1)
end


@inline function _dan_dq(p::T, q::T, f::T, n::Int) where {T<:AbstractFloat}
    # ∂a_n/∂q avoiding the removable singularity at q ≈ n for integer q.
    # For n=1, defer to the specific a₁ derivative.
    if n == 1
        return _da1_dq(p, q, f)
    end
    # Use the simplified closed-form of a_n that eliminates explicit q^2 via f:
    #   a_n = (x/(1-x))^2 * (n-1) * (p+n-1) * (p+q+n-2) * (q-n) / D(p,n)
    # where D(p,n) = (p+2n-3)*(p+2n-2)^2*(p+2n-1) and (x/(1-x)) = p*f/q.
    # Differentiate only the q-dependent factor G(q) = (p+q+n-2)*(q-n):
    #   dG/dq = (q-n) + (p+q+n-2) = p + 2q - 2.

    # This is equivalent to  
    #   return _anfun(p,q,f,n) * (inv(p+q+n-2) + inv(q-n))
    # but more precise.

    pfq = (p * f) / q
    C   = (pfq * pfq) * (n - 1) * (p + n - 1) /
          ((p + 2*n - 3) * (p + 2*n - 2)^2 * (p + 2*n - 1))
    return C * (p + 2*q - 2)
end

@inline function _dbn_dp(p::T, q::T, f::T, n::Int) where {T<:AbstractFloat}
    # ∂b_n/∂p via quotient rule on b_n = N/D.
    # Note the internal dependence f(p,q)=q x/(p(1-x)) — terms cancel in N as per derivation.
    g = p * f + 2 * q
    A = 2 * n^2 + 2 * (p - 1) * n
    N1 = g * A
    N2 = p * q * (p - 2 - p * f)
    N = N1 + N2
    D = q * (p + 2*n - 2) * (p + 2*n)
    dN1_dp = 2 * n * g
    dN2_dp = q * (2 * p - 2) - p * q * f
    dN_dp = dN1_dp + dN2_dp
    dD_dp = q * (2 * p + 4 * n - 2)
    return (dN_dp * D - N * dD_dp) / (D^2)
end

@inline function _dbn_dq(p::T, q::T, f::T, n::Int) where {T<:AbstractFloat}
    # ∂b_n/∂q similarly via quotient rule
    g = p * f + 2 * q
    A = 2 * n^2 + 2 * (p - 1) * n
    N1 = g * A
    N2 = p * q * (p - 2 - p * f)
    N = N1 + N2
    D = q * (p + 2*n - 2) * (p + 2*n)
    g_q = p * (f / q) + 2
    dN1_dq = g_q * A
    dN2_dq = p * (p - 2 - p * f) - p^2 * f
    dN_dq = dN1_dq + dN2_dq
    dD_dq = (p + 2*n - 2) * (p + 2*n)
    return (dN_dq * D - N * dD_dq) / (D^2)
end

@inline function _nextapp(f::T, p::T, q::T, n::Int, App::T, Ap::T, Bpp::T, Bp::T) where {T<:AbstractFloat}
    # One step of the continuant recurrences:
    #   A_n = a_n A_{n-2} + b_n A_{n-1}
    #   B_n = a_n B_{n-2} + b_n B_{n-1}
    an = _anfun(p, q, f, n)
    bn = _bnfun(p, q, f, n)
    An = an * App + bn * Ap
    Bn = an * Bpp + bn * Bp
    return An, Bn, an, bn
end

@inline function _dnextapp(an::T, bn::T, dan::T, dbn::T, Xpp::T, Xp::T, dXpp::T, dXp::T) where {T<:AbstractFloat}
    # Derivative propagation for the same recurrences (X∈{A,B})
    return dan * Xpp + an * dXpp + dbn * Xp + bn * dXp
end

function _beta_inc_grad(a::T, b::T, x::T; maxapp::Int=200, minapp::Int=3, err::T=eps(T)*T(1e4)) where {T<:AbstractFloat}
    # Compute I_x(a,b) and partial derivatives (∂I/∂a, ∂I/∂b, ∂I/∂x)
    # using a differentiated continued fraction with convergence control.
    oneT = one(T)
    zeroT = zero(T)

    # 1) Boundary cases for x
    isone(x) && return oneT, zeroT, zeroT, zeroT
    iszero(x) && return zeroT, zeroT, zeroT, zeroT

    # 2) Get tolerence
    ϵ = err

    # 3) Non-boundary path: precompute ∂I/∂x at original (a,b,x) via stable log form
    dx = exp((a - oneT) * log(x) + (b - oneT) * log1p(-x) - logbeta(a,b))

    # 4) Optional tail-swap for symmetry and improved CF convergence:
    #    if x > a/(a+b), evaluate at (p,q,x₀) = (b,a,1-x) and swap back at the end.
    p    = a
    q    = b
    x₀   = x
    swap = false
    if x > a / (a + b)
        x₀   = oneT - x
        p    = b
        q    = a
        swap = true
    end
    
    # 5) Initialize CF state and derivatives
    K                    = _Kfun(x₀, p, q)
    dK_dp_val, dK_dq_val = _dK_dpdq(x₀, p, q)
    f                    = _ffun(x₀, p, q)
    App                  = oneT
    Ap                   = oneT
    Bpp                  = zeroT
    Bp                   = oneT
    dApp_dp              = zeroT
    dBpp_dp              = zeroT
    dAp_dp               = zeroT
    dBp_dp               = zeroT
    dApp_dq              = zeroT
    dBpp_dq              = zeroT
    dAp_dq               = zeroT
    dBp_dq               = zeroT
    dI_dp                = T(NaN)
    dI_dq                = T(NaN)
    Ixpq                 = T(NaN)
    Ixpqn                = T(NaN)
    dI_dp_prev           = T(NaN)
    dI_dq_prev           = T(NaN)
    
    # 6) Main CF loop (n from 1): update continuants, scale, form current approximant Cn=A_n/B_n
    #    and its derivatives to update I and ∂I/∂(p,q). Stop on relative convergence of all.
    for n=1:maxapp
        
        # Update continuants. 
        An, Bn, an, bn = _nextapp(f, p, q, n, App, Ap, Bpp, Bp)
        dan            = _dan_dp(p, q, f, n)
        dbn            = _dbn_dp(p, q, f, n)
        dAn_dp         = _dnextapp(an, bn, dan, dbn, App, Ap, dApp_dp, dAp_dp)
        dBn_dp         = _dnextapp(an, bn, dan, dbn, Bpp, Bp, dBpp_dp, dBp_dp)
        dan            = _dan_dq(p, q, f, n)
        dbn            = _dbn_dq(p, q, f, n)
        dAn_dq         = _dnextapp(an, bn, dan, dbn, App, Ap, dApp_dq, dAp_dq)
        dBn_dq         = _dnextapp(an, bn, dan, dbn, Bpp, Bp, dBpp_dq, dBp_dq)

        # Normalize states to control growth/underflow (scale-invariant transform)
        s = maximum((abs(An), abs(Bn), abs(Ap), abs(Bp), abs(App), abs(Bpp)))
        if isfinite(s) && s > zeroT
            invs     = inv(s)
            An      *= invs
            Bn      *= invs
            Ap      *= invs
            Bp      *= invs
            App     *= invs
            Bpp     *= invs
            dAn_dp  *= invs
            dBn_dp  *= invs
            dAn_dq  *= invs
            dBn_dq  *= invs
            dAp_dp  *= invs
            dBp_dp  *= invs
            dApp_dp *= invs
            dBpp_dp *= invs
            dAp_dq  *= invs
            dBp_dq  *= invs
            dApp_dq *= invs
            dBpp_dq *= invs
        end

        # Form current approximant Cn=A_n/B_n and its derivatives.
        # Guard against tiny/zero Bn to avoid NaNs/Inf in divisions.
        tiny   = sqrt(eps(T))
        absBn  = abs(Bn)
        sgnBn  = ifelse(Bn >= zeroT, oneT, -oneT)
        invBn  = absBn > tiny && isfinite(absBn) ? inv(Bn) : inv(sgnBn * tiny)
        Cn     = An * invBn
        invBn2 = invBn * invBn
        dI_dp  = dK_dp_val * Cn + K * (invBn * dAn_dp - (An * invBn2) * dBn_dp)
        dI_dq  = dK_dq_val * Cn + K * (invBn * dAn_dq - (An * invBn2) * dBn_dq)
        Ixpqn  = K * Cn

        # Decide convergence: 
        if n >= minapp
            # Relative convergence for I, ∂I/∂p, ∂I/∂q (guards against tiny denominators)
            denomI = max(abs(Ixpqn), abs(Ixpq), eps(T))
            denomp = max(abs(dI_dp), abs(dI_dp_prev), eps(T))
            denomq = max(abs(dI_dq), abs(dI_dq_prev), eps(T))
            rI     = abs(Ixpqn - Ixpq) / denomI
            rp     = abs(dI_dp - dI_dp_prev) / denomp
            rq     = abs(dI_dq - dI_dq_prev) / denomq
            if max(rI, rp, rq) < ϵ
                break
            end
        end
        Ixpq       = Ixpqn
        dI_dp_prev = dI_dp
        dI_dq_prev = dI_dq

        # Shift CF state for next iteration
        App        = Ap
        Bpp        = Bp
        Ap         = An
        Bp         = Bn
        dApp_dp    = dAp_dp
        dApp_dq    = dAp_dq
        dBpp_dp    = dBp_dp
        dBpp_dq    = dBp_dq
        dAp_dp     = dAn_dp
        dAp_dq     = dAn_dq
        dBp_dp     = dBn_dp
        dBp_dq     = dBn_dq
    end
    
    # 7) Undo tail-swap if applied; ∂I/∂x is the pdf at original (a,b,x)
    if swap
        return oneT - Ixpqn, -dI_dq, -dI_dp, dx
    else
        return Ixpqn, dI_dp, dI_dq, dx
    end
end



 

# Incomplete beta: beta_inc(a,b,x) -> (p, q) with q=1-p
function ChainRulesCore.frule((_, Δa, Δb, Δx), ::typeof(beta_inc), a::Number, b::Number, x::Number)
    # primal
    p, q = beta_inc(a, b, x)
    # derivatives
    _a, _b, _x = map(float, promote(a, b, x))
    _, dIa, dIb, dIx = _beta_inc_grad(_a, _b, _x)
    Δp = muladd(dIx, Δx, muladd(dIb, Δb, dIa * Δa))
    Δq = -Δp
    Tout = typeof((p, q))
    return (p, q), ChainRulesCore.Tangent{Tout}(Δp, Δq)
end

function ChainRulesCore.rrule(::typeof(beta_inc), a::Number, b::Number, x::Number)
    p, q = beta_inc(a, b, x)
    Ta = ChainRulesCore.ProjectTo(a)
    Tb = ChainRulesCore.ProjectTo(b)
    Tx = ChainRulesCore.ProjectTo(x)
    _a, _b, _x = map(float, promote(a, b, x))
    _, dIa, dIb, dIx = _beta_inc_grad(_a, _b, _x)
    function beta_inc_pullback(Δ)
        Δp, Δq = Δ
        s = Δp - Δq # because q = 1 - p
        ā = Ta(s * dIa)
        b̄ = Tb(s * dIb)
        x̄ = Tx(s * dIx)
        return ChainRulesCore.NoTangent(), ā, b̄, x̄
    end
    return (p, q), beta_inc_pullback
end
function ChainRulesCore.frule((_, Δa, Δb, Δx, Δy), ::typeof(beta_inc), a::Number, b::Number, x::Number, y::Number)
    p, q = beta_inc(a, b, x, y)
    _a, _b, _x, _y = map(float, promote(a, b, x, y))
    _, dIa, dIb, dIx = _beta_inc_grad(_a, _b, _x)
    Δp = muladd(dIx, Δx, muladd(-dIx, Δy, muladd(dIb, Δb, dIa * Δa)))
    Δq = -Δp
    Tout = typeof((p, q))
    return (p, q), ChainRulesCore.Tangent{Tout}(Δp, Δq)
end

function ChainRulesCore.rrule(::typeof(beta_inc), a::Number, b::Number, x::Number, y::Number)
    p, q = beta_inc(a, b, x, y)
    Ta = ChainRulesCore.ProjectTo(a)
    Tb = ChainRulesCore.ProjectTo(b)
    Tx = ChainRulesCore.ProjectTo(x)
    Ty = ChainRulesCore.ProjectTo(y)
    _a, _b, _x, _y = map(float, promote(a, b, x, y))
    _, dIa, dIb, dIx = _beta_inc_grad(_a, _b, _x)
    function beta_inc_pullback(Δ)
        Δp, Δq = Δ
        s = Δp - Δq
        ā = Ta(s * dIa)
        b̄ = Tb(s * dIb)
        x̄ = Tx(s * dIx)
        ȳ = Ty(-s * dIx)
        return ChainRulesCore.NoTangent(), ā, b̄, x̄, ȳ
    end
    return (p, q), beta_inc_pullback
end

# Inverse incomplete beta: beta_inc_inv(a,b,p) -> (x, 1-x)
function ChainRulesCore.frule((_, Δa, Δb, Δp), ::typeof(beta_inc_inv), a::Number, b::Number, p::Number)
    x, y = beta_inc_inv(a, b, p)
    _a, _b, _x, _p = map(float, promote(a, b, x, p))
    # Implicit differentiation at solved x: I_x(a,b) = p
    _, dIa, dIb, _ = _beta_inc_grad(_a, _b, _x)
    # ∂I/∂x at solved x via stable log-space expression
    dIx_acc = exp(muladd(_a - 1, log(_x), muladd(_b - 1, log1p(-_x), -logbeta(_a, _b))))
    inv_dIx = inv(dIx_acc)
    dx_da = -dIa * inv_dIx
    dx_db = -dIb * inv_dIx
    dx_dp = inv_dIx
    Δx = muladd(dx_dp, Δp, muladd(dx_db, Δb, dx_da * Δa))
    Δy = -Δx
    Tout = typeof((x, y))
    return (x, y), ChainRulesCore.Tangent{Tout}(Δx, Δy)
end

function ChainRulesCore.rrule(::typeof(beta_inc_inv), a::Number, b::Number, p::Number)
    x, y = beta_inc_inv(a, b, p)
    Ta = ChainRulesCore.ProjectTo(a)
    Tb = ChainRulesCore.ProjectTo(b)
    Tp = ChainRulesCore.ProjectTo(p)
    _a, _b, _x, _p = map(float, promote(a, b, x, p))
    _, dIa, dIb, _ = _beta_inc_grad(_a, _b, _x)
    # ∂I/∂x at solved x via stable log-space expression
    dIx_acc = exp(muladd(_a - 1, log(_x), muladd(_b - 1, log1p(-_x), -logbeta(_a, _b))))
    inv_dIx = inv(dIx_acc)
    dx_da = -dIa * inv_dIx
    dx_db = -dIb * inv_dIx
    dx_dp = inv_dIx
    function beta_inc_inv_pullback(Δ)
        Δx, Δy = Δ
        s = Δx - Δy
        ā = Ta(s * dx_da)
        b̄ = Tb(s * dx_db)
        p̄ = Tp(s * dx_dp)
        return ChainRulesCore.NoTangent(), ā, b̄, p̄
    end
    return (x, y), beta_inc_inv_pullback
end

end # module
