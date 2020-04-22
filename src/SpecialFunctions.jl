module SpecialFunctions

using OpenSpecFun_jll

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
    dawson,
    ellipk,
    ellipe,
    erf,
    erfc,
    erfcinv,
    erfcx,
    erfi,
    erfinv,
    logerfc,
    logerfcx,
    eta,
    digamma,
    invdigamma,
    polygamma,
    trigamma,
    gamma_inc,
    beta_inc,
    beta_inc_inv,
    gamma_inc_inv,
    ncbeta,
    ncF,
    hankelh1,
    hankelh1x,
    hankelh2,
    hankelh2x,
    zeta,
    sinint,
    cosint,
    chebyshevt,
    chebyshevu,
    hermiteh,
    laguerrel,
    legendrep,
    legendreq,
    lbinomial

include("bessel.jl")
include("erf.jl")
include("ellip.jl")
include("sincosint.jl")
include("gamma.jl")
include("gamma_inc.jl")
include("betanc.jl")
include("beta_inc.jl")
include("deprecated.jl")
include("14_legendre.jl")


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

end # module
