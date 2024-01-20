# SPDX-License-Identifier: BSD-3-Clause OR MIT

"""
    uchk(y::ComplexF64, ascle::Float64, tol::Float64)

`y` enters as a scaled quantity whose magnitude is greater than
`Exp(-alim) = ascle = 1.0e+3*D1_MACH(1)/tol`.
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

"""
    gammaln(z::Float64)

Compute the logarithm of the gamma function.

# fortran comments
DGAMLN COMPUTES THE NATURAL LOG OF THE GAMMA FUNCTION FOR
Z.GT.0.  THE ASYMPTOTIC EXPANSION IS USED TO GENERATE VALUES
GREATER THAN ZMIN WHICH ARE ADJUSTED BY THE RECURSION
G(Z+1)=Z*G(Z) FOR Z.LE.ZMIN.  THE FUNCTION WAS MADE AS
PORTABLE AS POSSIBLE BY COMPUTIMG ZMIN FROM THE NUMBER OF BASE
10 DIGITS IN A WORD, RLN=AMAX1(-ALOG10(R1MACH(4)),0.5E-18)
LIMITED TO 18 DIGITS OF (RELATIVE) ACCURACY.

SINCE INTEGER ARGUMENTS ARE COMMON, A TABLE LOOK UP ON 100
VALUES IS USED FOR SPEED OF EXECUTION.

# Impl Ref
- [`openspecfun/amos/dgamln.f`](https://github.com/JuliaMath/openspecfun/blob/v0.5.6/amos/dgamln.f)
"""
function gammaln(z::Float64)
    con::Float64 = 1.83787706640934548  # ln(2*pi)
    @assert con == log(big"2" * pi) |> Float64

    # TODO:  handle z <= 0 || isnan(z)  first
    if z > 0.0
        if z <= 101.0
            nz = Int(trunc(z))
            fz::Float64 = z - nz
            if fz <= 0.0
                if nz <= 100
                    return GAMMALN_GLN[nz]
                end
            end
        end #= 10 =#

        wdtol = max(D1_MACH[3+1], 1e-18)
        i1m = I1_MACH[13+1]
        rln = D1_MACH[4+1] * i1m
        fln = max(min(rln, 20.0), 3.0) - 3.0
        zm = 1.8 + 0.3875 * fln
        mz = trunc(zm) + 1
        zmin = mz
        zdmy = z
        zinc = 0.0
        if z < zmin
            zinc = zmin - nz
            zdmy = z + zinc
        end #= 20 =#

        zp = 1.0 / zdmy
        t1 = GAMMALN_CF[0+1] * zp
        s = t1
        if zp >= wdtol
            zsq = zp * zp
            tst = t1 * wdtol
            # `i=1` is outside `if`
            for i = 2:22
                zp *= zsq
                trm = GAMMALN_CF[i] * zp
                if abs(trm) < tst
                    break  # Goto 40
                end

                s += trm
            end #= 30 =#
        end #= 40 =#

        if 0.0 == zinc
            tlg = log(z)
            return z * (tlg - 1.0) + 0.5 * (con - tlg) + s
        end #= 50 =#

        zp = 1.0
        nz = trunc(zinc)
        for i = 0:(nz-1)
            zp *= z + i
        end #= 60 =#

        tlg = log(zdmy)
        return zdmy * (tlg - 1.0) - log(zp) + 0.5 * (con - tlg) + s
    end #= 70 =#

    # Zero or negative argument
    return NaN
end
