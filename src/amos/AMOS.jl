module AMOS
include("const.jl")
include("helper.jl")
include("warp.jl")


# internal subroutines
_subroutine_names = [
# No deps, leaf Functions
    "uchk",
    "gammaln",
    "s1s2",
# Only deps on leaf functions

# Othsers
] # _subroutine_names

for fname in _subroutine_names
    include(joinpath("subroutines", "$(fname).jl"))
end

end # AMOS
