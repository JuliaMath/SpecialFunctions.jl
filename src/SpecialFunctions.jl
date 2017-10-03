__precompile__()

module SpecialFunctions

using Compat

if VERSION >= v"0.7.0-DEV.1760"
    depsfile = joinpath(dirname(@__FILE__), "..", "deps", "deps.jl")
    if isfile(depsfile)
        include(depsfile)
    else
        error("SpecialFunctions is not properly installed. Please run " *
              "Pkg.build(\"SpecialFunctions\") and restart Julia.")
    end
else
    using Base.Math: openspecfun
end

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

include("bessel.jl")
include("erf.jl")
include("sincosint.jl")
include("gamma.jl")
include("deprecated.jl")

end # module
