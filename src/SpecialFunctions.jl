__precompile__()

module SpecialFunctions

# On older versions these are exported from Base
# TODO: Uncomment me once these have been removed!
#if VERSION >= v"0.6.0-dev.XXXX"
#    export
#        airyai,
#        airyaiprime,
#        airybi,
#        airybiprime,
#        airyaix,
#        airyaiprimex,
#        airybix,
#        airybiprimex,
#        besselh,
#        besselhx,
#        besseli,
#        besselix,
#        besselj,
#        besselj0,
#        besselj1,
#        besseljx,
#        besselk,
#        besselkx,
#        bessely,
#        bessely0,
#        bessely1,
#        besselyx,
#        dawson,
#        erf,
#        erfc,
#        erfcinv,
#        erfcx,
#        erfi,
#        erfinv,
#        eta,
#        digamma,
#        invdigamma,
#        polygamma,
#        trigamma,
#        hankelh1,
#        hankelh1x,
#        hankelh2,
#        hankelh2x,
#        zeta
#end

# On older versions libopenspecfun is included in Base
# TODO: Uncomment me once Base no longer ships with it!
#if VERSION >= v"0.6.0-dev.XXXX"
#    const openspecfun = "libopenspecfun"
#else
    const openspecfun = Base.Math.openspecfun
#end

include("bessel.jl")
include("erf.jl")
include("gamma.jl")
include("deprecated.jl")

end # module
