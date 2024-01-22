# SPDX-License-Identifier: BSD-3-Clause OR MIT

"""

# fortran comments
zseri computes the i bessel function for real(z).ge.0.0 by
means of the power series for large cabs(z) in the
region cabs(z).le.2*sqrt(fnu+1). nz=0 is a normal return.
nz.gt.0 means that the last nz components were set to zero
due to underflow. nz.lt.0 means underflow occurred, but the
condition cabs(z).le.2*sqrt(fnu+1) was violated and the
computation must be completed in another routine with n=n-abs(nz).

# Impl Ref
- [`openspecfun/amos/zseri.f`](https://github.com/JuliaMath/openspecfun/blob/v0.5.6/amos/zseri.f)
- [`scipy/scipy/special/_amos.c:amos_seri`](https://github.com/scipy/scipy/blob/b882f1b7ebe55e534f29a8d68a54e4ecd30aeb1a/scipy/special/_amos.c#L4178-L4345)

NOTE: Current implementation uses `@GOTO`!
"""
function seri!(
    y::Vector{ComplexF64},
    z::ComplexF64,
    fnu::Float64,
    kode::Int,
    n::Int,
    tol::Float64,
    elim::Float64, 
    alim::Float64
)
    @assert length(y) >= n

    w = ComplexF64[ 0.0, 0.0 ]
    nz::Int = 0
    az = abs(z)
    if az == 0.0
        #= 160 =#
        y[1] = 0.0 + 0.0im
        if fnu == 0.0
            y[1] = 1.0 + 0.0im
        end
        #= 170 =#
        if n == 1
            return nz
        end

        for i in 2:n
            y[i] = 0.0 + 0.0im
        end
        #= 180 =#
        return nz
    end

    arm = 1.0e3 * D1_MACH[1]
    rtr1 = sqrt(arm)
    crsc = 1.0
    iflag::Int = 0
    if az >= arm
        half_z = 0.5 * z
        cz = 0.0 + 0.0im
        if az > rtr1
            cz = half_z * half_z
        end #= 10 =#
        acz = abs(cz)
        nn::Int = n
        ck = log(half_z)

@label _line20
        #= 20 =#
        dfnu = fnu + (nn - 1)
        fnup = dfnu + 1.0
        #
        # UNDERFLOW TEST
        #
        ak1 = ck * dfnu
        ak = gammaln(fnup)
        ak1 -= ak
        if kode == 2
            ak1 -= real(z)
        end
        rak1 = real(ak1)
        if rak1 > (-elim)
            @goto _line40
        end

@label _line30
        #= 30 =#
        nz += 1
        y[nn] = 0.0
        if acz > dfnu
            #= 190 =#
            #
            # RETURN WITH NZ.LT.0 IF CABS(Z*Z/4).GT.FNU+N-NZ-1 COMPLETE
            # THE CALCULATION IN CBINU WITH N=N-IABS(NZ)
            #
            return -nz
        end
        nn -= 1
        if nn == 0
            return nz
        end
        @goto _line20

@label _line40
        #= 40 =#
        if rak1 <= -alim
            iflag = 1
            ss = 1.0 / tol
            crsc = tol
            ascle = arm * ss
        end
        #= 50 =#
        ak = imag(ak1)
        aa = exp(rak1)
        if iflag == 1
            aa *= ss
        end
        coef = aa * complex(cos(ak), sin(ak))
        atol = tol * acz / fnup
        il = min(2, nn)
        for i in 1:il
            dfnu = fnu + (nn - i)
            fnup = dfnu + 1.0
            s1 = 1.0 + 0.0im
            if acz >= (tol * fnup)
                ak1 = 1.0 + 0.0im
                ak = fnup + 2.0
                s = fnup
                aa = 2.0
                #= 60 =#
                while true
                    rs = 1.0 / s
                    ak1 *= cz
                    ak1 *= rs
                    s1 += ak1
                    s += ak
                    ak += 2.0
                    aa *= acz
                    aa *= rs
                    if aa <= atol
                        break
                    end
                end
            end
            #= 70 =#
            s2 = s1 * coef
            w[i] = s2
            if iflag != 0
                if amos_uchk(s2, ascle, tol)
                    @goto _line30
                end
            end
            #= 80 =#
            m = nn - i + 1
            y[m] = s2 * crsc
            if i != il
                coef *= dfnu / half_z
            end
        end
        #= 90 =#
        if nn <= 2
            return nz
        end
        
        k = nn - 2
        ak = Float64(k)
        rz = 2.0 / z
        if iflag == 1
            @goto _line120
        end
        ib = 3

@label _line100
        #= 100 =#
        for i in ib:nn
            y[k] = (ak + fnu) * rz * y[k + 1] + y[k + 2]
            ak -= 1.0
            k -= 1
        end
        #= 110 =#
        return nz

        #
        # RECUR BACKWARD WITH SCALED VALUES
        #
@label _line120
        #= 120 =#
        #
        # EXP(-ALIM)=EXP(-ELIM)/TOL=APPROX. ONE PRECISION ABOVE THE
        # UNDERFLOW LIMIT = ASCLE = D1MACH(1)*SS*1.0D+3
        #
        s1 = w[1]
        s2 = w[2]
        for l in 3:nn
            ck = s2
            s2 = s1 + (ak + fnu) * rz * s2
            s1 = ck
            ck = s2 * crsc
            y[k] = ck
            ak -= 1.0
            k -= 1
            if abs(ck) > ascle
                @goto _line140
            end
        end
        #= 130 =#
        return nz

@label _line140
        #= 140 =#
        ib = l + 1
        if ib > nn
            return nz
        end
        @goto _line100
    end

    
    #= 150 =#
    nz = n
    if fnu == 0.0
        nz -= 1
    end

    #= 160 =#
    y[1] = 0.0 + 0.0im
    if fnu == 0.0
        y[1] = 1.0 + 0.0im
    end
    #= 170 =#
    if n == 1
        return nz
    end

    for i in 2:n
        y[i] = 0.0 + 0.0im
    end
    #= 180 =#
    return nz
end
