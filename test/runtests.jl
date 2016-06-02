using SpecialFunctions
using Base.Test

# override Main
for f in names(SpecialFunctions)
    if isdefined(Main,f)
        @eval global const $f = SpecialFunctions.$f
    end
end

include("bessel.jl")
