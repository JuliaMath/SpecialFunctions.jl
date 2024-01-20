# SPDX-License-Identifier: BSD-3-Clause OR MIT

"""
    uchk(y::ComplexF64, ascle::Float64, tol::Float64)

`y` enters as a scaled quantity whose magnitude is greater than
`Exp(-alim) = ascle = 1.0e+3*d1mach(1)/tol`.
The test is made to see if the magnitude of the real or imaginary part would
Underflow when `y` is scaled (by `tol`) to its proper value.
`y` is accepted if the underflow is at least one precision below the magnitude of the largest component;
Otherwise the phase angle does not have absolute accuracy and an underflow is assumed.

# Impl Ref
- [`openspecfun/amos/zuchk.f`](https://github.com/JuliaMath/openspecfun/blob/v0.5.6/amos/zuchk.f)
"""
function uchk(y::ComplexF64, ascle::Float64, tol::Float64)
    yr = abs(real(y))
    yi = abs(imag(y))
    ss = max(yr, yi)
    st = min(yr, yi)
    
    if st > ascle
        # TODO: ret ::Bool
        return 0
    else
        st /= tol
        if ss < st
            return 1
        else
            return 0
        end
    end
end
