__precompile__()

module SpecialFunctions

# On other versions these are exported from Base
# if VERSION >= v"0.6.0-dev.XXXX"
#     export
#         airyai,
#         airyaiprime,
#         airybi,
#         airybiprime,
#         airyaix,
#         airyaiprimex,
#         airybix,
#         airybiprimex,
#         besselh,
#         besselhx,
#         besseli,
#         besselix,
#         besselj,
#         besselj0,
#         besselj1,
#         besseljx,
#         besselk,
#         besselkx,
#         bessely,
#         bessely0,
#         bessely1,
#         besselyx,
#         eta,
#         hankelh1,
#         hankelh1x,
#         hankelh2,
#         hankelh2x,
#         polygamma,
#         zeta
# end

include("bessel.jl")
include("erf.jl")
include("gamma.jl")

end # module
