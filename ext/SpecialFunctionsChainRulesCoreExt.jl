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
# This is an independent implementation adapted from the MIT-licensed
# https://github.com/arzwa/IncBetaDer. Its author also explicitly permitted reuse at
# https://github.com/JuliaMath/SpecialFunctions.jl/pull/506#issuecomment-3369042219.
# Copyright (c) 2025 Arthur Zwaenepoel

# Generic-typed helpers used by the continued-fraction evaluation of I_x(a,b)
# and its partial derivatives. These implement the scalar prefactor K(x;p,q),
# the auxiliary variable f, the continued-fraction coefficients a_n, b_n, and
# their partial derivatives w.r.t. p (≡ a) and q (≡ b). See Boik & Robinson-Cox (1999).

function _Kfun(logx::T, log1mx::T, p::T, q::T, logbetapq::T) where {T}
    # K(x;p,q) = x^p (1-x)^{q-1} / (p * B(p,q)) computed in log-space for stability
    # logx = log(x), log1mx = log(1-x), precomputed
    return exp(p * logx + (q - 1) * log1mx - log(p) - logbetapq)
end

function _ffun(x::T, p::T, q::T) where {T}
    # f = q x / (p (1-x)) — convenience variable appearing in CF coefficients
    return q * x / (p * (1 - x))
end

function _dK_dp(logx::T, p::T, K::T, ψpq::T, ψp::T) where {T}
    # ∂K/∂p using digamma identities: d/dp log B(p,q) = ψ(p) - ψ(p+q)
    # Near x = p/(p+q), cancellation in this bracket grows like p*eps(Float64)
    # and eventually limits the accuracy of the parameter derivative prefactor.
    return K * (logx - inv(p) + ψpq - ψp)
end

function _dK_dq(log1mx::T, K::T, ψpq::T, ψq::T) where {T}
    # ∂K/∂q using identical pattern
    return K * (log1mx + ψpq - ψq)
end

function _dK_dpdq(logx::T, log1mx::T, p::T, q::T, K::T) where {T}
    # Convenience: compute (∂K/∂p, ∂K/∂q) together with shared ψ(p+q)
    ψ = digamma(p + q)
    dKdp = _dK_dp(logx, p, K, ψ, digamma(p))
    dKdq = _dK_dq(log1mx, K, ψ, digamma(q))
    return dKdp, dKdq
end

# n=1 case
function _nextapp1(f::T, p::T, q::T) where {T}
    # One step of the continuant recurrences:
    #   A_n = a_n A_{n-2} + b_n A_{n-1}
    #   B_n = a_n B_{n-2} + b_n B_{n-1}
    an = p * f * (q - 1) / (q * (p + 1))
    bn = (2p*f / q + 2 + p * (1 - f)) / (p + 2)
    An = an + bn
    return An, bn
end

# This recurrence is evaluated four times per continued-fraction iteration. Keeping
# the shared helper inline avoids the slowdown measured in the discussion at
# https://github.com/JuliaMath/SpecialFunctions.jl/pull/506#discussion_r3521698360.
@inline function _dnextapp(an::T, bn::T, dan::T, dbn::T, Xpp::T, Xp::T, dXpp::T, dXp::T) where {T}
    # Derivative propagation for the same recurrences (X∈{A,B})
    return dan * Xpp + an * dXpp + dbn * Xp + bn * dXp
end

function _beta_inc_grad(a::Float64, b::Float64, x::Float64, y::Float64;
                        maxapp::Int=10_000, minapp::Int=10)
    # Compute I_x(a,b) and partial derivatives (∂I/∂a, ∂I/∂b, ∂I/∂x)
    # using a differentiated continued fraction. Boik & Robinson-Cox used a minimum
    # of 3 and a maximum of 200 approximants. The larger limit here is empirical: the
    # loop still exits as soon as all three values converge, while difficult, central,
    # nearly symmetric cases can require thousands of approximants. With this limit,
    # symmetric cases stop converging around a = b = 1.35e8 and strongly asymmetric
    # large shapes around min(a,b) = 6.8e7 (the exact thresholds are platform- and
    # tolerance-dependent). Highly imbalanced shapes can also need more than 10,000
    # approximants, e.g. (a,b,x) = (1e-6,100,1e-7).
    # The stopping tolerance follows the original Float64 implementation.
    # Reaching the iteration ceiling without meeting the convergence test returns NaNs;
    # the warning is limited because this function is often called in loops.
    ϵ = 1e4 * eps(Float64)
    oneT = 1.0
    zeroT = 0.0

    # Match the primal's short-circuit for invalid floating-point inputs. In
    # particular, do not spend all maxapp iterations on a recurrence of NaNs.
    if !isfinite(a) || !isfinite(b) || !isfinite(x) || !isfinite(y)
        return NaN, NaN, NaN
    end

    # At either endpoint I_x is independent of a and b. The x derivative is the
    # endpoint limit of the beta density.
    if iszero(x)
        dx = a < oneT ? Inf : isone(a) ? b : zeroT
        return zeroT, zeroT, dx
    elseif iszero(y)
        dx = b < oneT ? Inf : isone(b) ? a : zeroT
        return zeroT, zeroT, dx
    elseif iszero(a) || iszero(b)
        # Parameter derivatives are undefined at degenerate shapes; the primal is
        # constant in x in the interior, hence its x derivative is zero.
        return NaN, NaN, zeroT
    end

    # Precompute log(x) and log(y) once, preserving the explicit complement supplied
    # to the four-argument API even when 1 - x would round differently.
    logx   = log(x)
    log1mx = log(y)

    # Precompute ∂I/∂x at original (a,b,x). Reuse the DiDonato-Morris
    # machinery from the primal to avoid cancellation for large shape parameters.
    logbetapq = logbeta(a, b)  # Time-consuming step; symmetric in a and b.
    dx = SpecialFunctions.beta_integrand(a, b, x, y, -logx - log1mx)

    # Optional tail-swap for symmetry and improved CF convergence:
    #    if x > a/(a+b), evaluate at (p,q,x₀) = (b,a,y) and swap back at the end.
    swap = x > a / (a + b)
    if swap
        x₀      = y
        p       = b
        q       = a
        logx₀   = log1mx    # log(1-x) = log(x₀)
        log1mx₀ = logx      # log(1-(1-x)) = log(x)
    else
        x₀      = x
        p       = a
        q       = b
        logx₀   = logx
        log1mx₀ = log1mx
    end

    # Initialize CF state and derivatives.
    K                    = _Kfun(logx₀, log1mx₀, p, q, logbetapq)
    dK_dp_val, dK_dq_val = _dK_dpdq(logx₀, log1mx₀, p, q, K)
    f                    = _ffun(x₀, p, q)

    # Precompute loop-invariant expressions (only depend on p, q, f). Although f
    # depends on both parameters, p*f = q*x₀/(1-x₀) is independent of p, and
    # p*f/q = x₀/(1-x₀) is independent of both p and q. Consequently pfq² is
    # constant with respect to the differentiated recurrence parameters; preserving
    # these cancellations keeps the coefficient partials compact.
    pf      = p * f
    pfq     = pf / q                    # p * f / q
    pfq2    = pfq * pfq                 # (p * f / q)^2
    pf_plus_2q = pf + 2 * q
    p_minus_2_minus_pf = p - 2 - pf
    pq_times_p_minus_2_minus_pf = p * q * p_minus_2_minus_pf
    p_plus_2q_minus_2 = p + 2 * q - 2
    a1      = p * f * (q - 1) / (q * (p + 1))           # a₁ coefficient
    da1_dp  = -a1 / (p + 1)             # ∂a₁/∂p
    da1_dq  = pfq / (p + 1)             # ∂a₁/∂q, including the removable q == 1 case

    # Update continuants.
    An, Bn     = _nextapp1(f, p, q)
    dBn_dq     = -pfq / (p + 2)
    dBn_dp     = dBn_dq * (2 - q) / (p + 2)
    dAn_dp     = da1_dp + dBn_dp
    dAn_dq     = da1_dq + dBn_dq

    # Form current approximant Cn=A_n/B_n and its derivatives.
    # Guard against tiny/zero Bn to avoid NaNs/Inf in divisions. This conservative
    # sqrt(eps) clamp sets the documented large-shape ceiling; a smaller threshold
    # would raises it substantially but would require separate accuracy validation.
    tiny = sqrt(eps(Float64))

    invBn  = abs(Bn) > tiny && isfinite(Bn) ? inv(Bn) : inv(copysign(tiny, Bn))
    Cn     = An * invBn
    invBn2 = invBn * invBn
    dI_dp  = dK_dp_val * Cn + K * (invBn * dAn_dp - (An * invBn2) * dBn_dp)
    dI_dq  = dK_dq_val * Cn + K * (invBn * dAn_dq - (An * invBn2) * dBn_dq)
    Ixpqn  = K * Cn
    Ixpq       = Ixpqn
    dI_dp_prev = dI_dp
    dI_dq_prev = dI_dq

    # Shift CF state for next iteration
    App      = oneT
    Bpp      = oneT
    Ap       = An
    Bp       = Bn
    dApp_dp  = zeroT
    dApp_dq  = zeroT
    dBpp_dp  = zeroT
    dBpp_dq  = zeroT
    dAp_dp   = dAn_dp
    dAp_dq   = dAn_dq
    dBp_dp   = dBn_dp
    dBp_dq   = dBn_dq

    # Main CF loop (n from 2): update continuants, scale, form current approximant Cn=A_n/B_n
    #    and its derivatives to update I and ∂I/∂(p,q). Stop on relative convergence of all.
    converged = false
    for n=2:maxapp

        # Continued-fraction coefficients a_n and b_n. These expressions live here
        # because they are only used at this single call site in the hot loop.
        pn = p + n
        p2n = pn + n
        denominator_a = (p2n - 3) * (p2n - 2)^2 * (p2n - 1)
        an = pfq2 * (n - 1) * (pn + q - 2) * (pn - 1) * (q - n) / denominator_a

        A = 2 * n^2 + 2 * (p - 1) * n
        N = pf_plus_2q * A + pq_times_p_minus_2_minus_pf
        D = q * (p2n - 2) * p2n
        bn = N / D

        # Partial derivatives of a_n. The q derivative is written without a division
        # by q - n, avoiding its removable singularity for integer q.
        dlog_an_dp = inv(p + q + n - 2) + inv(pn - 1) - inv(p2n - 3) -
                      2 * inv(p2n - 2) - inv(p2n - 1)
        dan_p = an * dlog_an_dp
        dan_q = pfq2 * (n - 1) * (pn - 1) * p_plus_2q_minus_2 / denominator_a

        # Partial derivatives of b_n, sharing N, D, and A between both directions.
        dN_dp = 2 * n * pf_plus_2q + q * (2 * p - 2) - q * pf
        dD_dp = q * (2 * p + 4 * n - 2)
        dN_dq = (pf / q + 2) * A + p * p_minus_2_minus_pf - p * pf
        dD_dq = (p2n - 2) * p2n
        D2 = D^2
        dbn_p = (dN_dp * D - N * dD_dp) / D2
        dbn_q = (dN_dq * D - N * dD_dq) / D2

        # Update the numerator and denominator continuants.
        An = an * App + bn * Ap
        Bn = an * Bpp + bn * Bp
        dAn_dq         = _dnextapp(an, bn, dan_q, dbn_q, App, Ap, dApp_dq, dAp_dq)
        dBn_dq         = _dnextapp(an, bn, dan_q, dbn_q, Bpp, Bp, dBpp_dq, dBp_dq)
        dAn_dp         = _dnextapp(an, bn, dan_p, dbn_p, App, Ap, dApp_dp, dAp_dp)
        dBn_dp         = _dnextapp(an, bn, dan_p, dbn_p, Bpp, Bp, dBpp_dp, dBp_dp)

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
        invBn  = abs(Bn) > tiny && isfinite(Bn) ? inv(Bn) : inv(copysign(tiny, Bn))
        Cn     = An * invBn
        dI_dp  = dK_dp_val * Cn + K * (dAn_dp - (An * invBn) * dBn_dp) * invBn
        dI_dq  = dK_dq_val * Cn + K * (dAn_dq - (An * invBn) * dBn_dq) * invBn
        Ixpqn  = K * Cn

        # Decide convergence:
        if n >= minapp
            # Mixed relative/absolute convergence for I, ∂I/∂p, and ∂I/∂q;
            # the eps floor guards tiny denominators, where the test becomes absolute.
            denomI = max(abs(Ixpqn), abs(Ixpq), eps(Float64))
            denomp = max(abs(dI_dp), abs(dI_dp_prev), eps(Float64))
            denomq = max(abs(dI_dq), abs(dI_dq_prev), eps(Float64))
            rI     = (Ixpqn - Ixpq) / denomI
            rp     = (dI_dp - dI_dp_prev) / denomp
            rq     = (dI_dq - dI_dq_prev) / denomq
            if -ϵ < rI < ϵ && -ϵ < rp < ϵ && -ϵ < rq < ϵ
                converged = true
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

    if !converged
        @warn "_beta_inc_grad reached maxapp approximants before convergence; returning NaNs" a b x maxapp minapp maxlog=1
        return NaN, NaN, NaN
    end

    # Undo tail-swap if applied; ∂I/∂x is the pdf at original (a,b,x).
    if swap
        return -dI_dq, -dI_dp, dx
    else
        return dI_dp,  dI_dq, dx
    end
end

function _beta_inc_grad(a::Float64, b::Float64, x::Float64; kwargs...)
    return _beta_inc_grad(a, b, x, 1 - x; kwargs...)
end

function _beta_inc_grad(a::T, b::T, x::T; kwargs...) where {T<:Union{Float16, Float32}}
    return map(T, _beta_inc_grad(Float64(a), Float64(b), Float64(x); kwargs...))
end

function _beta_inc_grad(a::T, b::T, x::T, y::T; kwargs...) where {T<:Union{Float16, Float32}}
    return map(T, _beta_inc_grad(Float64(a), Float64(b), Float64(x), Float64(y); kwargs...))
end

# Incomplete beta: beta_inc(a,b,x) -> (p, q) with q=1-p
ChainRulesCore.@scalar_rule(
    beta_inc(a::Number, b::Number, x::Number),
    @setup((dIa, dIb, dIx) = _beta_inc_grad(map(float, promote(a, b, x))...)),
    (dIa, dIb, dIx),
    (-dIa, -dIb, -dIx),
)
# Incomplete beta: beta_inc(a,b,x,y) -> (p, q) with y=1-x, q=1-p
ChainRulesCore.@scalar_rule(
    beta_inc(a::Number, b::Number, x::Number, y::Number),
    @setup(
        (_a, _b, _x, _y) = map(float, promote(a, b, x, y)),
        (dIa, dIb, dIx) = _beta_inc_grad(_a, _b, _x, _y),
    ),
    (dIa, dIb, dIx / 2, -dIx / 2),
    (-dIa, -dIb, -dIx / 2, dIx / 2),
)
# Inverse incomplete beta: beta_inc_inv(a,b,p) -> (x, 1-x)
ChainRulesCore.@scalar_rule(
    beta_inc_inv(a::Number, b::Number, p::Number),
    @setup(
        (_a, _b, _) = map(float, promote(a, b, p)),
        x = first(Ω),
        (dIa, dIb, dIx) = _beta_inc_grad(_a, _b, float(x)),
        inv_dIx = inv(dIx),
        da = iszero(dIa) ? zero(dIa) : -dIa * inv_dIx,
        db = iszero(dIb) ? zero(dIb) : -dIb * inv_dIx,
    ),
    (da, db, inv_dIx),
    (-da, -db, -inv_dIx),
)

# Inverse incomplete beta: beta_inc_inv(a,b,p,q) -> (x, 1-x), with q=1-p.
ChainRulesCore.@scalar_rule(
    beta_inc_inv(a::Number, b::Number, p::Number, q::Number),
    @setup(
        (_a, _b, _, _) = map(float, promote(a, b, p, q)),
        x = first(Ω),
        (dIa, dIb, dIx) = _beta_inc_grad(_a, _b, float(x)),
        inv_dIx = inv(dIx),
        da = iszero(dIa) ? zero(dIa) : -dIa * inv_dIx,
        db = iszero(dIb) ? zero(dIb) : -dIb * inv_dIx,
        dp = inv_dIx / 2,
    ),
    (da, db, dp, -dp),
    (-da, -db, -dp, dp),
)

end # module
