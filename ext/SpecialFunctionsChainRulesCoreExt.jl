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


# Note on incomplete beta derivatives implementation
# --------------------------------------------------
# The rules for the regularized incomplete beta I_x(a,b) and its inverse are
# implemented using a direct translation of the original S-PLUS/MATLAB code by
# Boik & Robinson-Cox. See:
#   Boik, R. J., & Robinson-Cox, J. F. (1998).
#   Derivatives of the incomplete beta function with respect to its parameters.
#   Computational Statistics & Data Analysis, 27(1), 85–106.
# The coefficient recurrences and derivative accumulation are ported verbatim
# (scalar form) from inc.beta.deriv.S/inbeder.m.

@inline function _derconf_coeffs(n::Int, p::T, q::T, w::T) where {T<:AbstractFloat}
    F = w * q / p
    if n == 1
        t1 = 1 - inv(p + 1)
        t2 = 1 - inv(q)
        t3 = 1 - 2 / (p + 2)
        t4 = 1 - 2 / q
        an1 = t1 * t2 * F
        an2 = -an1 / (p + 1)
        an4 = t1 * F / q
        bn1 = 1 - t3 * t4 * F
        bn2 = t3 * t4 * F / (p + 2)
        bn4 = -t3 * F / q
        return an1, an2, an4, bn1, bn2, bn4
    end
    t2  = F^2
    t3  = 2n - 2
    t5  = p * q
    t7  = inv(t3 * q + t5)
    t8  = t2 * t7
    t9  = n^2
    t10 = t9^2
    t11 = t2 * t10
    t12 = 4n - 2
    t13 = q^2
    t14 = t12 * t13
    t15 = p * t13
    t17 = inv(t14 + 2t15)
    t19 = t9 * n
    t20 = t19 * t2
    t22 = inv(p + (2n - 1))
    t23 = t20 * t22
    t24 = 2n - 1
    t27 = inv(t24 * q + t5)
    t28 = t20 * t27
    t30 = t10 * n * t2
    t32 = n * t2
    t33 = 2n - 3
    t36 = inv(t33 * t13 + t15)
    t37 = t32 * t36
    t38 = t9 * t2
    t39 = inv(t13)
    t41 = t32 * t39
    t43 = (-8 + 4n) * n
    t47 = inv(4 + t43 + (4n - 4 + p) * p)
    t49 = t38 * t17
    t50 = t38 * t47
    t51 = t20 * t47
    t52 = inv(q)
    t54 = t2 * t47
    t55 = t32 * t47
    t57 = inv(2p + (4n - 6))
    t59 = 4t8 - 3t11 * t17 - 4t23 - t28 - 4t30 * t27 + 9t37 - t38 * t39 + t41 + 4t11 * t47 - t49 + 24t50 - 16t51 - t2 * t52 + 4t54 - 16t55 - 53t38 * t57
    t62 = inv(p + (2n - 2))
    t63 = t32 * t62
    t65 = inv(2p + (4n - 2))
    t69 = t2 * inv(p + (2n - 3))
    t70 = t69 * t19
    t73 = inv(t3 * t13 + t15)
    t74 = t11 * t73
    t76 = t10 * t9 * t2
    t79 = inv(t24 * t13 + t15)
    t81 = t2 * t62
    t82 = 4 + t43
    t84 = 4n - 4
    t89 = inv(t82 * t13 + (t84 * t13 + t15) * p)
    t91 = t20 * t36
    t92 = t11 * t27
    t96 = t20 * t89
    t97 = t20 * t7
    t98 = t12 * q
    t100 = inv(t98 + 2t5)
    t102 = 51t32 * t57 - 24t63 + 5t38 * t65 + 12t70 + 40t74 + 2t76 * t79 + 8t81 + 4t76 * t89 + 52t91 + 6t92 - 2t69 * t10 - 8t20 * t62 + 2t11 * t22 - 16t96 - 64t97 + t32 * t100
    t104 = t38 * t62
    t105 = t30 * t36
    t107 = 4n - 6
    t108 = t107 * q
    t110 = inv(t108 + 2t5)
    t113 = t38 * t73
    t116 = inv(t33 * q + t5)
    t117 = t11 * t116
    t118 = t20 * t116
    t119 = t30 * t79
    t120 = t32 * t73
    t122 = t20 * t73
    t123 = t20 * t79
    t126 = 24t104 + 14t105 + t32 * t52 + 87t32 * t110 - 9t69 - 12t30 * t73 + 24t113 - 26t117 + 65t118 - 2t119 - 4t120 + 4t30 * t116 - 48t122 + 2t123 - 2t76 * t36 - 3t38 * t100
    t132 = inv(t82 * q + (t84 * q + t5) * p)
    t133 = t20 * t132
    t135 = t38 * t89
    t136 = t11 * t89
    t137 = t30 * t89
    t138 = t11 * t132
    t139 = t107 * t13
    t141 = inv(t139 + 2t15)
    t142 = t38 * t141
    t143 = t32 * t132
    t144 = t32 * t7
    t145 = t38 * t7
    t149 = t38 * t132
    t151 = t2 * t116
    t152 = -48t133 - 8t30 * t132 + 4t135 + 24t136 - 16t137 + 32t138 - 69t142 - 8t143 - 32t144 + 72t145 - t32 * t65 + 20t11 * t7 - 77t11 * t141 + 32t149 - 155t38 * t110 - 9t151
    an1 = t59 + t102 + t126 + t152
    # an2 (∂/∂p)
    t155 = (4n - 4) * n
    t156 = 1 + t155
    t161 = inv(t156 * t13 + (t14 + t15) * p)
    t162 = t30 * t161
    t163 = -8 + 8n
    t164 = t163 * n
    t165 = 2 + t164
    t167 = -4 + 8n
    t172 = inv(t165 * t13 + (t167 * t13 + 2t15) * p)
    t175 = (-24 + 8n) * n
    t179 = inv(18 + t175 + (-12 + 8n + 2p) * p)
    t181 = t20 * t161
    t182 = t38 * t22
    t184 = (24 + t175) * n
    t186 = (-24 + 12n) * n
    t192 = inv(-8 + t184 + (12 + t186 + (-6 + 6n + p) * p) * p)
    t198 = inv(t156 * q + (t98 + t5) * p)
    t199 = t11 * t198
    t200 = t20 * t192
    t201 = -4t8 + 2t162 + 3t11 * t172 - 51t32 * t179 + 2t23 + 4t28 - 2t181 - 3t182 - 8t11 * t192 - 6t199 + 32t200 - 6t37
    t207 = inv(t165 * q + (t167 * q + 2t5) * p)
    t210 = (-12 + 4n) * n
    t211 = 9 + t210
    t216 = inv(t211 * t13 + (t139 + t15) * p)
    t217 = t32 * t216
    t218 = -8 + t184
    t220 = 12 + t186
    t222 = -6 + 6n
    t229 = inv(t218 * t13 + (t220 * t13 + (t222 * t13 + t15) * p) * p)
    t230 = t11 * t229
    t231 = t20 * t216
    t232 = t69 * n
    t233 = t30 * t216
    t234 = 18 + t175
    t236 = -12 + 8n
    t241 = inv(t234 * t13 + (t236 * t13 + 2t15) * p)
    t242 = t38 * t241
    t243 = 3t38 * t207 - 36t50 + 12t51 - 12t54 - 9t217 + 36t55 + 12t63 - 48t230 - 52t231 - 13t232 - 14t233 + 69t242
    t245 = t32 * t192
    t251 = inv(t234 * q + (t236 * q + 2t5) * p)
    t256 = inv(1 + t155 + (4n - 2 + p) * p)
    t257 = t20 * t256
    t258 = 32t245 - 2t70 - 10t74 - 6t81 - 22t91 - 4t92 + 60t96 + 16t97 - 6t104 - 87t32 * t251 - 2t105 + 4t257
    t267 = inv(t218 * q + (t220 * q + (t222 * q + t5) * p) * p)
    t268 = t11 * t267
    t269 = t11 * t79
    t270 = t30 * t229
    t271 = t32 * t267
    t272 = 6t69 - 64t268 - 18t113 + 4t117 - 20t118 - t269 + 32t270 + 2t119 + 4t120 + 24t122 - 2t123 + 16t271
    t276 = t32 * t27
    t277 = t69 * t9
    t278 = t38 * t116
    t279 = t38 * t192
    t281 = 77t11 * t241 - t276 + 88t133 - 28t135 - 52t136 + 16t137 + 9t277 + 35t278 - 28t138 - 48t279 + 40t143 + 155t38 * t251
    t286 = inv(t211 * q + (t108 + t5) * p)
    t287 = t20 * t286
    t288 = t2 * t192
    t292 = inv(9 + t210 + (4n - 6 + p) * p)
    t293 = t2 * t292
    t294 = t2 * t286
    t295 = t20 * t267
    t296 = t2 * t132
    t297 = t32 * t89
    t299 = 24t144 - 36t145 - 96t149 - 65t287 + 6t151 - 8t288 + 9t293 + 9t294 + 96t295 - 4t296 + 4t297 - 4t30 * t286
    t304 = t11 * t286
    t305 = t32 * t116
    t308 = t38 * t267
    t309 = t11 * t36
    t311 = t38 * t79
    t315 = inv(2 + t164 + (-4 + 8n + 2p) * p)
    t317 = 2t11 * t292 - t32 * t207 - 2t11 * t256 + 26t304 - 25t305 + 4t30 * t198 + 16t30 * t267 - 64t308 + 11t309 - 8t76 * t229 + t311 - 5t38 * t315
    t319 = t32 * t22
    t320 = t20 * t198
    t321 = t20 * t292
    t322 = t38 * t229
    t323 = t38 * t27
    t324 = t20 * t229
    t328 = t38 * t36
    t329 = t38 * t172
    t330 = t32 * t315 + t319 + t320 - 12t321 - 8t322 + t323 + 32t324 - 2t76 * t161 + 2t76 * t216 + 53t38 * t179 + 19t328 + t329
    an2 = t201 + t243 + t258 + t272 + t281 + t299 + t317 + t330
    # an4
    t521  = 16t8 - 8t28 + t41 - 3t49 + 20t74 + 65t91 + 4t92 - 48t96 - 16t97 + 4t105 + 72t113 - 4t117 + 24t118 +
            6t269 - 4t119 - 32t120 - 64t122 - t123 - t276 - 32t133
    t526  = t2 * t73
    t527  = t2 * t36
    t528  = 48t149 - 18t151 + 8t296 - 8t297 + 51t305 - 26t309 + 5t323 + t32 * t17 + 87t32 * t141 + 4t526 - 9t527
    an4 = t521 + 32t135 + 32t136 - 8t137 - t2 * t39 - 53t278 + 8t138 - 155t142 - 32t143 - 48t144 + 48t145 + t528

    # bn1, bn2, bn4
    t544 = t9 * F
    t546 = inv(p + 2n)
    t548 = q * n
    t550 = inv(t5 + 2t548)
    t551 = t544 * t550
    t552 = t544 * t7
    t553 = n * F
    t554 = t553 * t7
    t555 = t19 * F
    t557 = F * t62
    t559 = t557 * n
    bn1 = 1 - F + 2t544 * t546 - 2t551 - 4t552 + 2t554 + 2t555 * t7 - 2t557 - 2t557 * t9 + 4t559 - 2t555 * t550 + 2t553 * t52
    t563  = t553 * t550
    t564  = t553 * t132
    t567  = t544 * t132
    t568  = F * t47
    t572  = inv(4 * t9 + (4n + p) * p)
    t574  = q * t9
    t578  = inv(4 * t574 + (4 * t548 + t5) * p)
    t580  = t544 * t578
    t582  = t553 * t47
    bn2 = -t563 - 2t564 + 2t544 * t47 - 2t555 * t132 + 4t567 + 2t568 - 2t544 * t572 + 2t555 * t578 - t551 + 2t580 + t552 - t554 + t557 - t559 + t553 * t546 - 4t582
    bn4 = -F * t52 - 2t552 + 4t554 - 2(F * t7) + 2t551
    return an1, an2, an4, bn1, bn2, bn4
end

function _ibeta_grad_splus(a::T, b::T, x::T; maxapp::Int=200, minapp::Int=3, err::T=eps(T)*T(1e4)) where {T<:AbstractFloat}
    if x <= zero(T)
        return zero(T), zero(T), zero(T), zero(T)
    elseif x >= one(T)
        return one(T), zero(T), zero(T), zero(T)
    end
    # ∂I/∂x at original params
    dI_dx = exp(muladd(a - one(T), log(x), muladd(b - one(T), log1p(-x), -logbeta(a, b))))
    # psi
    lbet = logbeta(a, b)
    pa = digamma(a); pa1 = trigamma(a)
    pb = digamma(b); pb1 = trigamma(b)
    pab = digamma(a + b); pab1 = trigamma(a + b)
    # possibly swap
    x1 = x; omx = one(T) - x; pp = a; qq = b
    swapped = false
    if x > a / (a + b)
        swapped = true
        x1 = one(T) - x
        omx = x
        pp, qq = b, a
        pa, pb = pb, pa
        pa1, pb1 = pb1, pa1
    end
    w = x1 / omx
    logx1 = log(x1); logomx = log(omx)
    cc1 = muladd(pp, logx1, muladd(qq - one(T), logomx, -lbet - log(pp)))
    c0 = exp(cc1)
    cc2 = logx1 - inv(pp) - pa + pab
    cc4 = logomx - pb + pab
    # init recurrences
    an1_1 = one(T); an1_p = zero(T); an1_q = zero(T)
    an2_1 = one(T); an2_p = zero(T); an2_q = zero(T)
    bn1_1 = one(T); bn1_p = zero(T); bn1_q = zero(T)
    bn2_1 = zero(T); bn2_p = zero(T); bn2_q = zero(T)
    I = zero(T); Ip = zero(T); Iq = zero(T)
    prevI = T(NaN); prevIp = T(NaN); prevIq = T(NaN)
    d = one(T); n = 0
    while (n < minapp) || ((d >= err) && (n < maxapp))
        n += 1
        a1, ap, aq, b1, bp, bq = _derconf_coeffs(n, pp, qq, w)
        # forward recurrences
        dan1 = a1 * an2_1 + b1 * an1_1
        dbn1 = a1 * bn2_1 + b1 * bn1_1
        danp = ap * an2_1 + a1 * an2_p + bp * an1_1 + b1 * an1_p
        dbnp = ap * bn2_1 + a1 * bn2_p + bp * bn1_1 + b1 * bn1_p
        danq = aq * an2_1 + a1 * an2_q + bq * an1_1 + b1 * an1_q
        dbnq = aq * bn2_1 + a1 * bn2_q + bq * bn1_1 + b1 * bn1_q
        # scale
        Rn = dan1
        if abs(dbn1) > abs(dan1)
            Rn = dbn1
        end
        if Rn != 0
            invRn = inv(Rn)
            an1_1 *= invRn; an1_p *= invRn; an1_q *= invRn
            bn1_1 *= invRn; bn1_p *= invRn; bn1_q *= invRn
            danp *= invRn; dbnp *= invRn; danq *= invRn; dbnq *= invRn
            if abs(dbn1) > abs(dan1)
                dan1 *= invRn; dbn1 = one(T)
            else
                dbn1 *= invRn; dan1 = one(T)
            end
        else
            dbn1 = one(T); dan1 = one(T)
        end
        # approximant components
        dr1 = dan1 / dbn1
        drp = (danp - dr1 * dbnp) / dbn1
        drq = (danq - dr1 * dbnq) / dbn1
        # shift n-1/n-2
        an2_1, an2_p, an2_q = an1_1, an1_p, an1_q
        an1_1, an1_p, an1_q = dan1, danp, danq
        bn2_1, bn2_p, bn2_q = bn1_1, bn1_p, bn1_q
        bn1_1, bn1_p, bn1_q = dbn1, dbnp, dbnq
        # nth approximant
        pr = dr1 > 0 ? exp(cc1 + log(dr1)) : zero(T)
        I  = pr
        Ip = muladd(pr, cc2, c0 * drp)
        Iq = muladd(pr, cc4, c0 * drq)
        # convergence
        d1 = max(err, abs(I)); d2 = max(err, abs(Ip)); d4 = max(err, abs(Iq))
        dI = isfinite(prevI) ? abs(prevI - I) / d1 : one(T)
        dP = isfinite(prevIp) ? abs(prevIp - Ip) / d2 : one(T)
        dQ = isfinite(prevIq) ? abs(prevIq - Iq) / d4 : one(T)
        d = max(dI, max(dP, dQ))
        prevI, prevIp, prevIq = I, Ip, Iq
    end
    if swapped
        I = one(T) - I
        Ip, Iq = -Iq, -Ip
    end
    return I, Ip, Iq, dI_dx
end

 

# Incomplete beta: beta_inc(a,b,x) -> (p, q) with q=1-p
function ChainRulesCore.frule((_, Δa, Δb, Δx), ::typeof(beta_inc), a::Number, b::Number, x::Number)
    # primal
    p, q = beta_inc(a, b, x)
    # derivatives
    T = promote_type(float(typeof(a)), float(typeof(b)), float(typeof(x)))
    _, dIa, dIb, dIx = _ibeta_grad_splus(T(a), T(b), T(x))
    Δp = dIa * convert(T, Δa) + dIb * convert(T, Δb) + dIx * convert(T, Δx)
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
    _, dIa, dIb, dIx = _ibeta_grad_splus(T(a), T(b), T(x))
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
    _, dIa, dIb, dIx = _ibeta_grad_splus(T(a), T(b), T(x))
    Δp = dIa * convert(T, Δa) + dIb * convert(T, Δb) + dIx * (convert(T, Δx) - convert(T, Δy))
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
    _, dIa, dIb, dIx = _ibeta_grad_splus(T(a), T(b), T(x))
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
    _, dIa, dIb, _ = _ibeta_grad_splus(T(a), T(b), T(x))
    # ∂I/∂x at solved x via stable log-space expression
    dIx_acc = exp(muladd(T(a) - one(T), log(T(x)), muladd(T(b) - one(T), log1p(-T(x)), -logbeta(T(a), T(b)))))
    inv_dIx = inv(dIx_acc)
    dx_da = -dIa * inv_dIx
    dx_db = -dIb * inv_dIx
    dx_dp = inv_dIx
    Δx = dx_da * convert(T, Δa) + dx_db * convert(T, Δb) + dx_dp * convert(T, Δp)
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
    _, dIa, dIb, _ = _ibeta_grad_splus(T(a), T(b), T(x))
    # ∂I/∂x at solved x via stable log-space expression
    dIx_acc = exp(muladd(T(a) - one(T), log(T(x)), muladd(T(b) - one(T), log1p(-T(x)), -logbeta(T(a), T(b)))))
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
