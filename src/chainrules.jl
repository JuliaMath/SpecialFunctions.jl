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
    invdigamma(x),
    inv(trigamma(invdigamma(x))),
)
ChainRulesCore.@scalar_rule(trigamma(x), polygamma(2, x))

# binary
ChainRulesCore.@scalar_rule(
    besselj(ν, x),
    (
        ChainRulesCore.@thunk(error("not implemented")),
        (besselj(ν - 1, x) - besselj(ν + 1, x)) / 2
    ),
)
ChainRulesCore.@scalar_rule(
    besseli(ν, x),
    (
        ChainRulesCore.@thunk(error("not implemented")),
        (besseli(ν - 1, x) + besseli(ν + 1, x)) / 2,
    ),
)
ChainRulesCore.@scalar_rule(
    bessely(ν, x),
    (
        ChainRulesCore.@thunk(error("not implemented")),
        (bessely(ν - 1, x) - bessely(ν + 1, x)) / 2,
    ),
)
ChainRulesCore.@scalar_rule(
    besselk(ν, x),
    (
        ChainRulesCore.@thunk(error("not implemented")),
        -(besselk(ν - 1, x) + besselk(ν + 1, x)) / 2,
    ),
)
ChainRulesCore.@scalar_rule(
    hankelh1(ν, x),
    (
        ChainRulesCore.@thunk(error("not implemented")),
        (hankelh1(ν - 1, x) - hankelh1(ν + 1, x)) / 2,
    ),
)
ChainRulesCore.@scalar_rule(
    hankelh2(ν, x),
    (
        ChainRulesCore.@thunk(error("not implemented")),
        (hankelh2(ν - 1, x) - hankelh2(ν + 1, x)) / 2,
    ),
)
ChainRulesCore.@scalar_rule(
    polygamma(m, x),
    (
        ChainRulesCore.@thunk(error("not implemented")),
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
ChainRulesCore.@scalar_rule(logabsgamma(x), digamma(x), ChainRulesCore.Zero())

ChainRulesCore.@scalar_rule(loggamma(x), digamma(x))
