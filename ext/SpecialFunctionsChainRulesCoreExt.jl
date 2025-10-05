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

# Generic-typed version for high-precision evaluation
function _beta_inc_grad_boik(a::T, b::T, x::T,
                             maxapp::Int=200, minapp::Int=3, ϵ::T=convert(T, 1e-12)) where {T<:AbstractFloat}
    oneT = one(T); zeroT = zero(T)
    if x == oneT
        return oneT, zeroT, zeroT, zeroT
    elseif x == zeroT
        return zeroT, zeroT, zeroT, zeroT
    end
    dx = exp((a - oneT) * log(x) + (b - oneT) * log1p(-x) - logbeta(a,b))
    # swap tails if necessary
    p = a; q = b; x₀ = x; swap = false
    if x > a / (a + b)
        x₀ = oneT - x
        p = b
        q = a
        swap = true
    end
    Kfun(x::T, p::T, q::T) = exp(p * log(x) + (q - oneT) * log1p(-x) - log(p) - logbeta(p, q))
    ffun(x::T, p::T, q::T) = q*x/(p*(oneT - x))
    a1fun(p::T, q::T, f::T) = p*f*(q - oneT)/(q*(p + oneT))
    anfun(p::T, q::T, f::T, n::Int) = n == 1 ? a1fun(p, q, f) : 
        p^2 * f^2 * (T(n) - oneT) * (p + q + T(n) - T(2)) * (p + T(n) - oneT) * (q - T(n)) / 
            (q^2 * (p + T(2n) - T(3)) * (p + T(2n) - T(2))^2 * (p + T(2n) - oneT))
    function bnfun(p::T, q::T, f::T, n::Int)
        x = T(2)*(p*f + T(2)*q)*T(n)^2 + T(2)*(p*f + T(2)*q)*(p - oneT)*T(n) + p*q*(p - T(2) - p*f)
        y = (q * (p + T(2n) - T(2)) * (p + T(2n)))
        x/y
    end
    dK_dp(x::T, p::T, q::T, K::T, ψpq::T, ψp::T) = K*(log(x) - inv(p) + ψpq - ψp)
    dK_dq(x::T, p::T, q::T, K::T, ψpq::T, ψq::T) = K*(log1p(-x) + ψpq - ψq)
    function dK_dpdq(x::T, p::T, q::T)
        ψ = digamma(p+q)
        Kf = Kfun(x, p, q)
        dKdp = dK_dp(x, p, q, Kf, ψ, digamma(p))
        dKdq = dK_dq(x, p, q, Kf, ψ, digamma(q))
        dKdp, dKdq
    end
    # a_n derivatives via log-derivative
    da1_dp(p::T, q::T, f::T) = -a1fun(p, q, f) / (p + oneT)
    function dan_dp(p::T, q::T, f::T, n::Int)
        if n == 1
            return da1_dp(p, q, f)
        end
        an = anfun(p, q, f, n)
        dlog = inv(p + q + T(n) - T(2)) + inv(p + T(n) - oneT) - inv(p + T(2n) - T(3)) - T(2) * inv(p + T(2n) - T(2)) - inv(p + T(2n) - oneT)
        return an * dlog
    end
    da1_dq(p::T, q::T, f::T) = a1fun(p, q, f) / (q - oneT)
    function dan_dq(p::T, q::T, f::T, n::Int)
        if n == 1
            return da1_dq(p, q, f)
        end
        an = anfun(p, q, f, n)
        dlog = inv(p + q + T(n) - T(2)) + inv(q - T(n))
        return an * dlog
    end
    # b_n derivatives via quotient rule, accounting for f_p=-f/p, f_q=f/q which cancel in N
    function dbn_dp(p::T, q::T, f::T, n::Int)
        g = p * f + T(2) * q
        A = T(2) * T(n)^2 + T(2) * (p - oneT) * T(n)
        N1 = g * A
        N2 = p * q * (p - T(2) - p * f)
        N = N1 + N2
        D = q * (p + T(2n) - T(2)) * (p + T(2n))
        dN1_dp = T(2) * T(n) * g
        dN2_dp = q * (T(2) * p - T(2)) - p * q * f
        dN_dp = dN1_dp + dN2_dp
        dD_dp = q * (T(2) * p + T(4) * T(n) - T(2))
        return (dN_dp * D - N * dD_dp) / (D^2)
    end
    function dbn_dq(p::T, q::T, f::T, n::Int)
        g = p * f + T(2) * q
        A = T(2) * T(n)^2 + T(2) * (p - oneT) * T(n)
        N1 = g * A
        N2 = p * q * (p - T(2) - p * f)
        N = N1 + N2
        D = q * (p + T(2n) - T(2)) * (p + T(2n))
        g_q = p * (f / q) + T(2)
        dN1_dq = g_q * A
        dN2_dq = p * (p - T(2) - p * f) - p^2 * f
        dN_dq = dN1_dq + dN2_dq
        dD_dq = (p + T(2n) - T(2)) * (p + T(2n))
        return (dN_dq * D - N * dD_dq) / (D^2)
    end
    _nextapp(f::T, p::T, q::T, n::Int, App::T, Ap::T, Bpp::T, Bp::T) = begin
        an = anfun(p, q, f, n)
        bn = bnfun(p, q, f, n)
        An = an*App + bn*Ap
        Bn = an*Bpp + bn*Bp
        An, Bn, an, bn
    end
    _dnextapp(an::T, bn::T, dan::T, dbn::T, Xpp::T, Xp::T, dXpp::T, dXp::T) = dan * Xpp + an * dXpp + dbn * Xp + bn * dXp

    # compute once
    K = Kfun(x₀, p, q)
    dK_dp_val, dK_dq_val = dK_dpdq(x₀, p, q)
    f = ffun(x₀, p, q)
    App = oneT; Ap  = oneT; Bpp = zeroT; Bp  = oneT
    dApp_dp = zeroT; dBpp_dp = zeroT; dAp_dp  = zeroT; dBp_dp  = zeroT
    dApp_dq = zeroT; dBpp_dq = zeroT; dAp_dq  = zeroT; dBp_dq  = zeroT
    dI_dp   = T(NaN); dI_dq   = T(NaN); Ixpq = T(NaN); Ixpqn = T(NaN); dI_dp_prev = T(NaN); dI_dq_prev = T(NaN)
    for n=1:maxapp
        An, Bn, an, bn = _nextapp(f, p, q, n, App, Ap, Bpp, Bp)
        dan = dan_dp(p, q, f, n); dbn = dbn_dp(p, q, f, n)
        dAn_dp = _dnextapp(an, bn, dan, dbn, App, Ap, dApp_dp, dAp_dp)
        dBn_dp = _dnextapp(an, bn, dan, dbn, Bpp, Bp, dBpp_dp, dBp_dp)
        dan = dan_dq(p, q, f, n); dbn = dbn_dq(p, q, f, n)
        dAn_dq = _dnextapp(an, bn, dan, dbn, App, Ap, dApp_dq, dAp_dq)
        dBn_dq = _dnextapp(an, bn, dan, dbn, Bpp, Bp, dBpp_dq, dBp_dq)
        # normalize states to control growth/underflow (scale-invariant)
        s = maximum((abs(An), abs(Bn), abs(Ap), abs(Bp), abs(App), abs(Bpp)))
        if isfinite(s) && s > zeroT
            invs = inv(s)
            An *= invs; Bn *= invs
            Ap *= invs; Bp *= invs
            App *= invs; Bpp *= invs
            dAn_dp *= invs; dBn_dp *= invs
            dAn_dq *= invs; dBn_dq *= invs
            dAp_dp  *= invs; dBp_dp  *= invs
            dApp_dp *= invs; dBpp_dp *= invs
            dAp_dq  *= invs; dBp_dq  *= invs
            dApp_dq *= invs; dBpp_dq *= invs
        end
        Cn = An/Bn
        dI_dp = dK_dp_val * Cn + K * (inv(Bn) * dAn_dp - (An/(Bn^2)) * dBn_dp)
        dI_dq = dK_dq_val * Cn + K * (inv(Bn) * dAn_dq - (An/(Bn^2)) * dBn_dq)
        Ixpqn = K * Cn
        if n >= minapp
            denomI = max(abs(Ixpqn), abs(Ixpq), eps(T))
            denomp = max(abs(dI_dp), abs(dI_dp_prev), eps(T))
            denomq = max(abs(dI_dq), abs(dI_dq_prev), eps(T))
            rI = abs(Ixpqn - Ixpq) / denomI
            rp = abs(dI_dp - dI_dp_prev) / denomp
            rq = abs(dI_dq - dI_dq_prev) / denomq
            if max(rI, rp, rq) < ϵ
                break
            end
        end
        Ixpq = Ixpqn
        dI_dp_prev = dI_dp
        dI_dq_prev = dI_dq
        App  = Ap; Bpp  = Bp; Ap   = An; Bp   = Bn
        dApp_dp = dAp_dp; dApp_dq = dAp_dq; dBpp_dp = dBp_dp; dBpp_dq = dBp_dq
        dAp_dp  = dAn_dp; dAp_dq  = dAn_dq; dBp_dp  = dBn_dp; dBp_dq  = dBn_dq
    end
    if swap
        return oneT - Ixpqn, -dI_dq, -dI_dp, dx
    else
        return Ixpqn, dI_dp, dI_dq, dx
    end
end

# Generic wrapper preserving the previous interface/name
function _ibeta_grad_splus(a::T, b::T, x::T; maxapp::Int=200, minapp::Int=3, err::T=eps(T)*T(1e4)) where {T<:AbstractFloat}
    tol = min(err, T(1e-14))
    maxit = max(1000, maxapp)
    minit = max(5, minapp)
    I, dIa, dIb, dIx = _beta_inc_grad_boik(a, b, x, maxit, minit, tol)
    return I, dIa, dIb, dIx
end

 

# Incomplete beta: beta_inc(a,b,x) -> (p, q) with q=1-p
function ChainRulesCore.frule((_, Δa, Δb, Δx), ::typeof(beta_inc), a::Number, b::Number, x::Number)
    # primal
    p, q = beta_inc(a, b, x)
    # derivatives
    T = promote_type(float(typeof(a)), float(typeof(b)), float(typeof(x)))
    _, dIa_, dIb_, dIx_ = _ibeta_grad_splus(T(a), T(b), T(x))
    dIa::T = dIa_; dIb::T = dIb_; dIx::T = dIx_
    ΔaT::T = Δa isa Real ? T(Δa) : zero(T)
    ΔbT::T = Δb isa Real ? T(Δb) : zero(T)
    ΔxT::T = Δx isa Real ? T(Δx) : zero(T)
    Δp = dIa * ΔaT + dIb * ΔbT + dIx * ΔxT
    Δq = -Δp
    Tout = typeof((p, q))
    return (p, q), ChainRulesCore.Tangent{Tout}(Δp, Δq)
end

function ChainRulesCore.rrule(::typeof(beta_inc), a::Number, b::Number, x::Number)
    p, q = beta_inc(a, b, x)
    Ta = ChainRulesCore.ProjectTo(a)
    Tb = ChainRulesCore.ProjectTo(b)
    Tx = ChainRulesCore.ProjectTo(x)
    T = promote_type(float(typeof(a)), float(typeof(b)), float(typeof(x)))
    _, dIa_, dIb_, dIx_ = _ibeta_grad_splus(T(a), T(b), T(x))
    dIa::T = dIa_; dIb::T = dIb_; dIx::T = dIx_
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
    T = promote_type(float(typeof(a)), float(typeof(b)), float(typeof(x)), float(typeof(y)))
    _, dIa_, dIb_, dIx_ = _ibeta_grad_splus(T(a), T(b), T(x))
    dIa::T = dIa_; dIb::T = dIb_; dIx::T = dIx_
    ΔaT::T = Δa isa Real ? T(Δa) : zero(T)
    ΔbT::T = Δb isa Real ? T(Δb) : zero(T)
    ΔxT::T = Δx isa Real ? T(Δx) : zero(T)
    ΔyT::T = Δy isa Real ? T(Δy) : zero(T)
    Δp = dIa * ΔaT + dIb * ΔbT + dIx * (ΔxT - ΔyT)
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
    T = promote_type(float(typeof(a)), float(typeof(b)), float(typeof(x)), float(typeof(y)))
    _, dIa_, dIb_, dIx_ = _ibeta_grad_splus(T(a), T(b), T(x))
    dIa::T = dIa_; dIb::T = dIb_; dIx::T = dIx_
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
    T = promote_type(float(typeof(a)), float(typeof(b)), float(typeof(p)))
    # Implicit differentiation at solved x: I_x(a,b) = p
    _, dIa_, dIb_, _ = _ibeta_grad_splus(T(a), T(b), T(x))
    dIa::T = dIa_; dIb::T = dIb_
    # ∂I/∂x at solved x via stable log-space expression
    dIx_acc::T = exp(muladd(T(a) - one(T), log(T(x)), muladd(T(b) - one(T), log1p(-T(x)), -logbeta(T(a), T(b)))))
    inv_dIx::T = inv(dIx_acc)
    dx_da::T = -dIa * inv_dIx
    dx_db::T = -dIb * inv_dIx
    dx_dp::T = inv_dIx
    ΔaT::T = Δa isa Real ? T(Δa) : zero(T)
    ΔbT::T = Δb isa Real ? T(Δb) : zero(T)
    ΔpT::T = Δp isa Real ? T(Δp) : zero(T)
    Δx = dx_da * ΔaT + dx_db * ΔbT + dx_dp * ΔpT
    Δy = -Δx
    Tout = typeof((x, y))
    return (x, y), ChainRulesCore.Tangent{Tout}(Δx, Δy)
end

function ChainRulesCore.rrule(::typeof(beta_inc_inv), a::Number, b::Number, p::Number)
    x, y = beta_inc_inv(a, b, p)
    Ta = ChainRulesCore.ProjectTo(a)
    Tb = ChainRulesCore.ProjectTo(b)
    Tp = ChainRulesCore.ProjectTo(p)
    T = promote_type(float(typeof(a)), float(typeof(b)), float(typeof(p)))
    _, dIa_, dIb_, _ = _ibeta_grad_splus(T(a), T(b), T(x))
    dIa::T = dIa_; dIb::T = dIb_
    # ∂I/∂x at solved x via stable log-space expression
    dIx_acc::T = exp(muladd(T(a) - one(T), log(T(x)), muladd(T(b) - one(T), log1p(-T(x)), -logbeta(T(a), T(b)))))
    inv_dIx::T = inv(dIx_acc)
    dx_da::T = -dIa * inv_dIx
    dx_db::T = -dIb * inv_dIx
    dx_dp::T = inv_dIx
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
