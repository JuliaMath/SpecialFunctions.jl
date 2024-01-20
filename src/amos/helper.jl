# SPDX-License-Identifier: BSD-3-Clause OR MIT

"""
    uchk(y::ComplexF64, ascle::Float64, tol::Float64)

Y ENTERS AS A SCALED QUANTITY WHOSE MAGNITUDE IS GREATER THAN
EXP(-ALIM)=ASCLE=1.0E+3*D1MACH(1)/TOL. THE TEST IS MADE TO SEE
IF THE MAGNITUDE OF THE REAL OR IMAGINARY PART WOULD UNDERFLOW
WHEN Y IS SCALED (BY TOL) TO ITS PROPER VALUE. Y IS ACCEPTED
IF THE UNDERFLOW IS AT LEAST ONE PRECISION BELOW THE MAGNITUDE
OF THE LARGEST COMPONENT; OTHERWISE THE PHASE ANGLE DOES NOT HAVE
ABSOLUTE ACCURACY AND AN UNDERFLOW IS ASSUMED.

# Impl Ref
- [`openspecfun/amos/zuchk.f`](https://github.com/JuliaMath/openspecfun/blob/v0.5.6/amos/zuchk.f)
"""
function uchk(y::ComplexF64, ascle::Float64, tol::Float64)
    yr = abs(real(y))
    yi = abs(imag(y))
    ss = max(yr, yi)
    st = min(yr, yi)
    
    if st > ascle
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
