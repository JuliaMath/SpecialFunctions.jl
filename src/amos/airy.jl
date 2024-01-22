# SPDX-License-Identifier: BSD-3-Clause OR MIT

"""

COMPUTE AIRY FUNCTIONS AI(Z) AND DAI(Z) FOR COMPLEX Z

# parameters
- `z` Complex Number
- `id`: ORDER OF DERIVATIVE
    id=0, or id=1
- `kode`: INDICATE THE SCALING OPTION
    kode=1, 
        AI=AI(Z)                ON ID=0 OR
        AI=DAI(Z)/DZ            ON ID=1
    kode=2, 
        AI=CEXP(ZTA)*AI(Z)      ON ID=0 OR
        AI=CEXP(ZTA)*DAI(Z)/DZ  ON ID=1 WHERE
        ZTA=(2/3)*Z*CSQRT(Z)

# Impl Ref
- [`openspecfun/amos/zairy.f`](https://github.com/JuliaMath/openspecfun/blob/v0.5.6/amos/zairy.f)
- [`scipy/scipy/special/_amos.c:amos_airy`](https://github.com/scipy/scipy/blob/b882f1b7ebe55e534f29a8d68a54e4ecd30aeb1a/scipy/special/_amos.c#L650-L994)
"""
function airy(z::ComplexF64, id::Int, kode::Int)
    # aa, ad, ak, alim, atrm, az, az3, bk, ck, dig, dk, d1, d2, elim, fid, fnu, rl, r1m5, sfac, tol,  bb, alaz = 0.0
    TTH = 2.0 / 3.0
    GAMMA_C1 = 0.35502805388781723926  # 1/(Gamma(2/3) * 3**(2/3))
    GAMMA_C2 = 0.25881940379280679841  # 1/(Gamma(1/3) * 3**(1/3))
    COEF = 0.18377629847393068317  # 1 / (sqrt(3) * PI)

    cy = [ 0.0 + 0.0im ]
    zta = 0.0 + 0.0im
    csq = 0.0 + 0.0im
    
    # TODO: return or throw
    ierr = 0
    nz = 0
    ai = 0.0 + 0.0im
    if (id < 0) || (id > 1)
        ierr = 1
    end
    if (kode < 1) || (kode > 2)
        ierr = 1
    end
    if ierr != 0
        # TODO: throw(AmosException(ierr)
        return 0.0 + 0.0im
    end

    az = abs(z)
    tol = max(D1_MACH[4], 1e-18)
    fid = float(id)

    if az <= 1.0
        #
        # POWER SERIES FOR ABS(Z) <= 1.
        #
        s1 = 1.0 + 0.0im
        s2 = 1.0 + 0.0im
        if az < tol
            #= 170 =#
            aa = 1.0e3 * D1_MACH[1]
            s1 = 0.0 + 0.0im
            if id != 1
                if az > aa
                    s1 = GAMMA_C2 * z
                end
                #= 180 =#
                ai = GAMMA_C1 - s1
                # @info "180: ai" nz ierr
                # @show csq s1 ai
                return ai
            end
        
            #= 190 =#
            ai = complex(-GAMMA_C2)
            aa = sqrt(aa)
            if az > aa
                s1 = z * z * 0.5
            end
            #= 200 =#
            ai += s1 * GAMMA_C1
            # @info "200: ai" nz ierr
            # @show csq s1 ai
            return ai
        end
    
        aa = az * az
        if aa >= (tol / az)
            trm1 = 1.0 + 0.0im
            trm2 = 1.0 + 0.0im
            atrm = 1.0
            z3 = z * z * z
            az3 = az * aa
            ak = 2.0 + fid
            bk = 3.0 - fid - fid
            ck = 4.0 - fid
            dk = 3.0 + fid + fid
            d1 = ak * dk
            d2 = bk * ck
            ad = min(d1, d2)
            ak = 24.0 + 9.0 * fid
            bk = 30.0 - 9.0 * fid
            for k = 1:25
                trm1 *= z3 / d1
                s1 += trm1
                trm2 *= z3 / d2
                s2 += trm2
                atrm *= az3 / ad
                d1 += ak
                d2 += bk
                ad = (d1 > d2 ? d2 : d1)
                if atrm < tol * ad
                    break
                end
                ak += 18.0
                bk += 18.0
            end
        end #= 30 =#
        #= 40 =#
        if id != 1
            ai = s1*GAMMA_C1 - z*s2*GAMMA_C2
            if kode == 1
                return ai
            end
        
            zta = z * sqrt(z) * TTH
            ai *= exp(zta)
            # @info "40: ai" nz ierr
            # @show zta
            return ai
        end
        #= 50 =#
        ai = -s2 * GAMMA_C2
        if az > tol
            ai += z * z * s1 * GAMMA_C1 / (1.0 + fid)
        end
        if kode == 1
            return ai
        end
    
        zta = z * sqrt(z) * TTH
        # @info "50: ai * exp(zta)" nz ierr
        # @show csq s1 ai
        return ai * exp(zta)
    end

    #= 70 =#
    #
    # CASE FOR CABS(Z).GT.1.0
    #
    fnu = (1.0 + fid) / 3.0
    #
    # SET PARAMETERS RELATED TO MACHINE CONSTANTS.
    # TOL IS THE APPROXIMATE UNIT ROUNDOFF LIMITED TO 1.0E-18.
    # ELIM IS THE APPROXIMATE EXPONENTIAL OVER- AND UNDERFLOW LIMIT.
    # EXP(-ELIM) < EXP(-ALIM)=EXP(-ELIM)/TOL    AND
    # EXP(ELIM) > EXP(ALIM)=EXP(ELIM)*TOL       ARE INTERVALS NEAR
    # UNDERFLOW AND OVERFLOW LIMITS WHERE SCALED ARITHMETIC IS DONE.
    # RL IS THE LOWER BOUNDARY OF THE ASYMPTOTIC EXPANSION FOR LARGE Z.
    # DIG = NUMBER OF BASE 10 DIGITS IN TOL = 10**(-DIG).
    #
    k1 = trunc(Int, I1_MACH[15])
    k2 = trunc(Int, I1_MACH[16])
    r1m5 = D1_MACH[5]
    k = min(abs(k1), abs(k2))
    elim = 2.303 * (k * r1m5 - 3.0)
    k1 = trunc(Int, I1_MACH[14]) - 1
    aa = r1m5 * k1
    dig = min(aa, 18.0)
    aa *= 2.303
    alim = elim + max(-aa, -41.45)
    rl = 1.2 * dig + 3.0
    alaz = log(az)

    #
    # TEST FOR PROPER RANGE
    #
    aa = 0.5 / tol
    bb = I1_MACH[9] * 0.5
    aa = min(aa, bb)
    aa = aa^(TTH)
    if az > aa
        #= 260 =#
        ierr = 4
        nz = 0
        # TODO: throw(AmosException(ierr)
        # @info "260: 0.0im" nz ierr
        # @show az aa
        return 0.0 + 0.0im
    end

    aa = sqrt(aa)
    if az > aa
        ierr = 3
    end
    csq = sqrt(z)
    zta = z * csq * TTH

    #
    # RE(ZTA).LE.0 WHEN RE(Z).LT.0, ESPECIALLY WHEN IM(Z) IS SMALL
    #
    iflag = 0
    sfac = 1.0
    ak = imag(zta)
    if real(z) < 0.0
        bk = real(zta)
        ck = -abs(bk)
        zta = complex(ck, ak)
    end
    #= 80 =#
    if (imag(z) == 0.0) && (real(z) <= 0.0)
        zta = complex(0.0, ak)
    end
    #= 90 =#
    aa = real(zta)
    if (aa < 0.0) || (real(z) <= 0.0)
        if kode != 2
            #
            # OVERFLOW TEST
            #
            if aa <= -alim
                aa = -aa + 0.25 * alaz
                iflag = 1
                sfac = tol
                if aa > elim
                    #= 270 =#
                    nz = 0
                    ierr = 2
                    # TODO: throw(AmosException(ierr)
                    # @info "90: ai" nz ierr iflag
                    # @show aa sfac
                    return ai
                end
            end
        end

        #= 100 =#
        #
        # CBKNU AND CACON RETURN EXP(ZTA)*K(FNU,ZTA) ON KODE=2
        #
        mr = 1
        if imag(z) < 0.0
            mr = -1
        end
        nn = acai!(cy, zta, fnu, kode, mr, 1, rl, tol, elim, alim)
        @assert length(cy) == 1
        if nn < 0
            #= 280 =#
            if nn == -1
                nz = 1
                # @info "280: 0.0" nz ierr iflag
                return 0.0
            else
                nz = 0
                ierr = 5
                # TODO: throw(AmosException(ierr)
                # @info "208-else: 0.0" nz ierr iflag
                return 0.0
            end
        end
    
        nz += nn
    else
        #= 110 =#
        if kode != 2
            #
            # OVERFLOW TEST
            #
            if aa >= alim
                aa = -aa - 0.25 * alaz
                iflag = 2
                sfac = 1.0 / tol
                if aa < -elim
                    nz = 1
                    # @info "110: 0.0" nz ierr iflag
                    return 0.0
                end
            end
        end
        #= 120 =#
        nz = bknu!(cy, zta, fnu, kode, 1, tol, elim, alim)
        @assert length(cy) == 1
    end

    #= 130 =#
    s1 = cy[1] * COEF

    # 0: normal;  3: underflow
    @assert ierr==0 || ierr==3
    if iflag == 0
        if id != 1
            # @info "130: csq * s1" nz ierr iflag
            # @show csq s1
            return csq * s1
        end
        #= 140 =#
        # @info "140: csq * s1" nz ierr iflag
        # @show csq s1
        return (-z * s1)
    end

    #= 150 =#
    s1 *= sfac
    if id != 1
        s1 *= csq
        # @info "150: s1 / sfac" nz ierr iflag
        # @show csq s1
        return (s1 / sfac)
    end
    
    #= 160 =#
    s1 *= -z
    # @info "160: s1 / sfac" nz ierr iflag
    # @show csq s1
    return (s1 / sfac)
end
