__precompile__()

module SpecialFunctions

if isdefined(Base, :airyai) && VERSION < v"0.7.0-DEV.986" #22763
    import Base: airyai, airyaix, airyaiprime, airyaiprimex,
                 airybi, airybix, airybiprime, airybiprimex,
                 besselh, besselhx, besseli, besselix, besselj, besselj0, besselj1,
                 besseljx, besselk, besselkx, bessely, bessely0, bessely1, besselyx,
                 hankelh1, hankelh1x, hankelh2, hankelh2x,
                 dawson, erf, erfc, erfcinv, erfcx, erfi, erfinv,
                 eta, digamma, invdigamma, polygamma, trigamma, zeta,
                 # deprecated
                 airy, airyx, airyprime
else
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
        besselj0,
        besselj1,
        besseljx,
        besselk,
        besselkx,
        bessely,
        bessely0,
        bessely1,
        besselyx,
        dawson,
        erf,
        erfc,
        erfcinv,
        erfcx,
        erfi,
        erfinv,
        eta,
        digamma,
        invdigamma,
        polygamma,
        trigamma,
        hankelh1,
        hankelh1x,
        hankelh2,
        hankelh2x,
        zeta
end

export sinint,
       cosint

if isdefined(Base.Math, :openspecfun)
    const openspecfun = Base.Math.openspecfun
else
    const openspecfun = "libopenspecfun"
end

# Avoids a deprecation warning for the 0-arg case on 0.7
if VERSION < v"0.7.0-DEV.924"
    const DomainErrorNoArgs = DomainError()
else
    const DomainErrorNoArgs = DomainError(nothing)
end

include("bessel.jl")
include("erf.jl")
include("sincosint.jl")
include("gamma.jl")
include("deprecated.jl")

end # module
