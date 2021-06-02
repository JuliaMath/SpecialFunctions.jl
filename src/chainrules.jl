const BESSEL_ORDER_INFO = """
derivatives of Bessel functions with respect to the order are not implemented currently:
https://github.com/JuliaMath/SpecialFunctions.jl/issues/160
"""

const INCOMPLETE_GAMMA_INFO = """
derivatives of the incomplete Gamma functions with respect to parameter `a` are not
implemented currently:
https://github.com/JuliaMath/SpecialFunctions.jl/issues/317
"""

ChainRulesCore.@scalar_rule(airyai(x), airyaiprime(x))
ChainRulesCore.@scalar_rule(airyaiprime(x), x * airyai(x))
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
ChainRulesCore.@scalar_rule(erf(x), (2 / sqrt(π)) * exp(-x * x))
ChainRulesCore.@scalar_rule(erfc(x), -(2 / sqrt(π)) * exp(-x * x))
ChainRulesCore.@scalar_rule(erfcinv(x), -(sqrt(π) / 2) * exp(Ω^2))
ChainRulesCore.@scalar_rule(erfcx(x), (2 * x * Ω) - (2 / sqrt(π)))
ChainRulesCore.@scalar_rule(erfi(x), (2 / sqrt(π)) * exp(x * x))
ChainRulesCore.@scalar_rule(erfinv(x), (sqrt(π) / 2) * exp(Ω^2))
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

# binary
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
    hankelh1(ν, x),
    (
        ChainRulesCore.@not_implemented(BESSEL_ORDER_INFO),
        (hankelh1(ν - 1, x) - hankelh1(ν + 1, x)) / 2,
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
