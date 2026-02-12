# SPDX-License-Identifier: BSD-3-Clause OR MIT

"""

# fortran comments
ZACAI APPLIES THE ANALYTIC CONTINUATION FORMULA
    
    K(FNU,ZN*EXP(MP))=K(FNU,ZN)*EXP(-MP*FNU) - MP*I(FNU,ZN)
            MP=PI*MR*CMPLX(0.0,1.0)

TO CONTINUE THE K FUNCTION FROM THE RIGHT HALF TO THE LEFT
HALF Z PLANE FOR USE WITH ZAIRY WHERE FNU=1/3 OR 2/3 AND N=1.
ZACAI IS THE SAME AS ZACON WITH THE PARTS FOR LARGER ORDERS AND
RECURRENCE REMOVED. A RECURSIVE CALL TO ZACON CAN RESULT IF ZACON
IS CALLED FROM ZAIRY.

# Impl Ref
- [`openspecfun/amos/zacai.f`](https://github.com/JuliaMath/openspecfun/blob/v0.5.6/amos/zacai.f)
- [`scipy/scipy/special/_amos.c:amos_acai`](https://github.com/scipy/scipy/blob/b882f1b7ebe55e534f29a8d68a54e4ecd30aeb1a/scipy/special/_amos.c#L383-L485)

NOTE: Current implementation uses `@GOTO`!
"""
function acai!(
    y::Vector{ComplexF64},
    z::ComplexF64,
    fnu::Float64,
    kode::Int,
    mr::Int,
    n::Int,
    rl::Float64,
    tol::Float64,
    elim::Float64, 
    alim::Float64
)
    int_nu = Int(trunc(fnu))

    cy = ComplexF64[ 0.0, 0.0 ]
    nz::Int = 0
    zn = -z
    az = abs(z)
    dfnu = fnu + (n - 1)
    if (az > 2.0) && ((az * az * 0.25) > (dfnu + 1.0))
        #= 20 =#
        if az >= rl
            #
            # Asymptotic expansion for large z for the I function
            #
            nw = asyi!(y, zn, fnu, kode, n, rl, tol, elim, alim)
        else
            #= 30 =#
            #
            # Miller algorithm normalized by the series for the I function
            #
            nw = mlri!(y, zn, fnu, kode, n, tol)
        end
        if nw < 0
            #= 80 =#
            nz = -1
            if nw == -2
                nz = -2
            end
            return nz
        end
    else
        #= 10 =#
        #
        # Power series for the I function
        #
        seri!(y, zn, fnu, kode, n, tol, elim, alim)
    end

    #= 40 =#
    #
    # Analytic continuation to the left half plane for the K function
    #
    nw = bknu!(cy, zn, fnu, kode, 1, tol, elim, alim)
    if nw != 0
        #= 80 =#
        nz = -1
        if nw == -2
            nz = -2
        end
        return nz
    end

    fmr = mr
    sgn = fmr < 0.0 ? pi : -pi
    csgn = complex(0.0, sgn)

    if kode != 1
        yy = -imag(zn)
        cpn = cos(yy)
        spn = sin(yy)
        csgn *= complex(cpn, spn)
    end

    #= 50 =#
    #
    # Calculate cspn = exp(fnu*pi*i) to minimize losses of significance
    # when fnu is large
    #
    arg = (fnu - int_nu) * sgn
    cpn = cos(arg)
    spn = sin(arg)
    cspn = complex(cpn, spn)
    if isodd(int_nu)
        cspn = -cspn
    end

    #= 60 =#
    c1 = Ref(cy[1])
    c2 = Ref(y[1])

    if kode != 1
        iuf = 0
        ascle = 1e3 * D1_MACH[1] / tol
        nw, iuf = s1s2!(zn, c1, c2, ascle, alim, iuf)
        nz += nw
    end

    #= 70 =#
    y[1] = cspn * c1[] + csgn * c2[]
    return nz
end
