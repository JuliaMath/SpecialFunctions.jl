module SpecialFunctions

# Load openspecfun libraries from our deps.jl
let depsjl_path = joinpath(@__DIR__, "..", "deps", "deps.jl")
    if !isfile(depsjl_path)
        error("SpecialFunctions is not installed properly, run `Pkg.build(\"SpecialFunctions\")`," *
              "restart Julia and try again")
    end
    include(depsjl_path)
end

__init__() = check_deps()

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
    zeta,
    sinint,
    cosint,
    lbinomial

include("bessel.jl")
include("erf.jl")
include("sincosint.jl")
include("gamma.jl")
include("deprecated.jl")

for f in (:digamma, :erf, :erfc, :erfcinv, :erfcx, :erfi, :erfinv,
          :eta, :gamma, :invdigamma, :lfactorial, :lgamma, :trigamma)
    @eval $(f)(::Missing) = missing
end
for f in (:beta, :lbeta)
    @eval $(f)(::Number, ::Missing) = missing
    @eval $(f)(::Missing, ::Number) = missing
    @eval $(f)(::Missing, ::Missing) = missing
end
polygamma(m::Integer, x::Missing) = missing

end # module
