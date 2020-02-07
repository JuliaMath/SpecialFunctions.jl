#prelude.jl

#some common imports used by almost all scripts:
using Base.MPFR: ROUNDING_MODE, big_ln2
using Base.Math: @horner, libm, nan_dom_err