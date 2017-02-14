__precompile__()

module SpecialFunctions

using Compat

if VERSION >= v"0.6.0-dev.2767"
    if isdefined(Base, :airyai)
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
end

if isdefined(Base.Math, :openspecfun)
    const openspecfun = Base.Math.openspecfun
else
    const openspecfun = "libopenspecfun"
end

include("bessel.jl")
include("erf.jl")
include("gamma.jl")
include("deprecated.jl")

end # module
