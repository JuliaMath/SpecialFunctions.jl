module SpecialFunctions

using IrrationalConstants:
    twoπ,
    halfπ,
    sqrtπ,
    sqrt2π,
    invπ,
    inv2π,
    invsqrt2,
    invsqrt2π,
    logtwo,
    logπ,
    log2π

import LogExpFunctions

using OpenLibm_jll: libopenlibm
using OpenSpecFun_jll: libopenspecfun

export
    airyai,
    airyaiprime,
    airybi,
    airybiprime,
    airyaix,
    airyaiprimex,
    airybix,
    airybiprimex,
    besselh,
    besselhx,
    besseli,
    besselix,
    besselj,
    sphericalbesselj,
    besselj0,
    besselj1,
    besseljx,
    besselk,
    besselkx,
    bessely,
    sphericalbessely,
    bessely0,
    bessely1,
    besselyx,
    jinc,
    dawson,
    ellipk,
    ellipe,
    erf,
    erfc,
    erfcinv,
    erfcx,
    erfi,
    erfinv,
    logerf,
    logerfc,
    logerfcx,
    faddeeva,
    eta,

    # Gamma functions
    gamma,
    loggamma,
    logabsgamma,
    logfactorial,
    digamma,
    invdigamma,
    trigamma,
    polygamma,
    gamma_inc,
    gamma_inc_inv,
    # beta functions
    beta,
    logbeta,
    logabsbeta,
    logabsbinomial,
    beta_inc,
    beta_inc_inv,

    ncbeta,
    ncF,
    hankelh1,
    hankelh1x,
    hankelh2,
    hankelh2x,
    zeta,
    expint,
    expinti,
    expintx,
    sinint,
    cosint,
    lbinomial

include("bessel.jl")
include("erf.jl")
include("ellip.jl")
include("expint.jl")
include("sincosint.jl")
include("logabsgamma/e_lgammaf_r.jl")
include("logabsgamma/e_lgamma_r.jl")
include("gamma.jl")
include("gamma_inc.jl")
include("betanc.jl")
include("beta_inc.jl")
include("deprecated.jl")

for f in (:digamma, :erf, :erfc, :erfcinv, :erfcx, :erfi, :erfinv, :logerfc, :logerfcx,
          :eta, :gamma, :invdigamma, :logfactorial, :lgamma, :trigamma, :ellipk, :ellipe)
    @eval $(f)(::Missing) = missing
end
for f in (:beta, :lbeta)
    @eval $(f)(::Number, ::Missing) = missing
    @eval $(f)(::Missing, ::Number) = missing
    @eval $(f)(::Missing, ::Missing) = missing
end
polygamma(m::Integer, x::Missing) = missing

# `fastabs` is identical to `Base.Math.fastabs` which is not used here since it is not public
# https://github.com/JuliaLang/julia/blob/93fb785831dcfcc442f82fab8746f0244c5274ae/base/special/trig.jl#L1057
fastabs(x::Number) = abs(x)
fastabs(x::Complex) = abs(real(x)) + abs(imag(x))

end # module
