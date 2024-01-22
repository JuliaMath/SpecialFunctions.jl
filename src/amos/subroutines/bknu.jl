# SPDX-License-Identifier: BSD-3-Clause OR MIT

"""

ZBKNU COMPUTES THE K BESSEL FUNCTION IN THE RIGHT HALF Z PLANE.

# Impl Ref
- [`openspecfun/amos/zbknu.f`](https://github.com/JuliaMath/openspecfun/blob/v0.5.6/amos/zbknu.f)
- [`scipy/scipy/special/_amos.c:amos_bknu`](https://github.com/scipy/scipy/blob/b882f1b7ebe55e534f29a8d68a54e4ecd30aeb1a/scipy/special/_amos.c#L3005-L3467)

NOTE: Current implementation uses `@GOTO`!
"""
function bknu!(
    y::Vector{ComplexF64},
    z::ComplexF64,
    fnu::Float64,
    kode::Int,
    n::Int,
    tol::Float64,
    elim::Float64, 
    alim::Float64
)
    s1 = 0.0 + 0.0im
    s2 = 0.0 + 0.0im
    ck = 0.0 + 0.0im
    dnu2 = 0.0

    r1 = 2.0
    tth = 2.0 / 3.0
    cc = Float64[
        5.77215664901532861e-01, -4.20026350340952355e-02,
        -4.21977345555443367e-02, 7.21894324666309954e-03,
        -2.15241674114950973e-04, -2.01348547807882387e-05,
        1.13302723198169588e-06, 6.11609510448141582e-09
    ]

    caz = abs(z)
    cscl = 1.0 / tol
    crsc = tol
    css = Float64[cscl, 1.0, crsc]
    csr = Float64[crsc, 1.0, cscl]
    bry = [
        1e3 * D1_MACH[1] / tol,
        tol / (1e3 * D1_MACH[1]),
        D1_MACH[2]
    ]
    nz = 0
    iflag = 0
    koded = kode
    rz = 2.0 / z
    inu = Int(trunc(fnu + 0.5))
    dnu = fnu - inu
    if abs(dnu) != 0.5
        dnu2 = 0.0
        if abs(dnu) > tol
            dnu2 = dnu * dnu
        end
        if caz <= r1
            #
            # SERIES FOR CABS(Z).LE.R1
            #
            fc = 1.0
            smu = log(rz)
            fmu = smu * dnu
            csh = sinh(fmu)
            cch = cosh(fmu)
            if dnu != 0.0
                fc = dnu * pi
                fc *= 1.0 / sin(fc)
                smu = csh / dnu
            end
            #= 10 =#
            a2 = 1.0 + dnu
            #
            # GAM(1-Z)*GAM(1+Z)=PI*Z/SIN(PI*Z), T1=1/GAM(1-DNU), T2=1/GAM(1+DNU)
            #
            t2 = exp(-gammaln(a2))
            t1 = 1.0 / (t2 * fc)
            if abs(dnu) <= 0.1
                #
                # SERIES FOR F0 TO RESOLVE INDETERMINACY FOR SMALL ABS(DNU)
                #
                ak = 1.0
                s = cc[1]
                for k = 2:8
                    ak *= dnu2
                    tm = cc[k] * ak
                    s += tm
                    if abs(tm) < tol
                        break
                    end
                end #= 20 =#
                #= 30 =#
                g1 = -s
            else
                #= 40 =#
                g1 = (t1 - t2) / (dnu + dnu)
            end
            #= 50 =#
            g2 = 0.5 * (t1 + t2)
            f = fc * (g1 * cch + smu * g2)
            pt = exp(fmu)
            p = (0.5 / t2) * pt
            q = (0.5 / t1) / pt
            s1 = f
            s2 = p
            ak = 1.0
            a1 = 1.0
            ck = 1.0 + 0.0im
            bk = 1.0 - dnu2
            if (inu <= 0) && (n <= 1)
                #
                # GENERATE K(FNU,Z), 0.0D0 .LE. FNU .LT. 0.5D0 AND N=1
                #
                if caz >= tol
                    cz = z * z * 0.25
                    t1 = 0.25 * caz * caz
                    #= 60 =#
                    while true
                        f = (f * ak + p + q) / bk
                        p = p / (ak - dnu)
                        q = q / (ak + dnu)
                        rk = 1.0 / ak
                        ck *= cz * rk
                        s1 += ck * f
                        a1 *= t1 * rk
                        bk += ak + ak + 1.0
                        ak += 1.0
                        if a1 <= tol
                            break
                        end
                    end
                end
                #= 70 =#
                y[1] = s1
                if koded == 1
                    return nz
                end
                
                y[1] = s1 * exp(z)
                return nz
            end
    
            #= 80 =#
            #
            # GENERATE K(DNU,Z) AND K(DNU+1,Z) FOR FORWARD RECURRENCE
            #
            if caz >= tol
                cz = z * z * 0.25
                t1 = 0.25 * caz * caz
                #= 90 =#
                while true
                    f = (f * ak + p + q) / bk
                    p *= 1.0 / (ak - dnu)
                    q *= 1.0 / (ak + dnu)
                    rk = 1.0 / ak
                    ck *= cz * rk
                    s1 += ck * f
                    s2 += ck * (p - f * ak)
                    a1 *= t1 * rk
                    bk += ak + ak + 1.0
                    ak += 1.0
                    if a1 <= tol
                        break
                    end
                end
            end
            #= 100 =#
            kflag = 2
            a1 = fnu + 1.0
            ak = a1 * abs(real(smu))
            if ak > alim
                kflag = 3
            end
            p2 = s2 * css[kflag]
            s2 = p2 * rz
            s1 *= css[kflag]
            if koded != 1
                f = exp(z)
                s1 *= f
                s2 *= f
            end
            @goto _line210
        end
    end

    #= 110 =#
    #
    #   IFLAG=0 MEANS NO UNDERFLOW OCCURRED
    #   IFLAG=1 MEANS AN UNDERFLOW OCCURRED- COMPUTATION PROCEEDS WITH
    #   KODED=2 AND A TEST FOR ON SCALE VALUES IS MADE DURING FORWARD
    #   RECURSION
    #
    coef = rthpi / sqrt(z)
    kflag = 2
    if koded != 2
        if real(z) > alim
            #= 290 =#
            #
            # SCALE BY DEXP(Z), IFLAG = 1 CASES
            #
            koded = 2
            iflag = 1
            kflag = 2
        else
            a1 = exp(-real(z)) * real(css[kflag])
            pt = a1 * complex(cos(imag(z)), -sin(imag(z)))
            coef *= pt
        end
    end
    #= 120 =#
    if abs(dnu) == 0.5
        #= 300 =#
        #
        # FNU=HALF ODD INTEGER CASE, DNU=-0.5
        #
        s1 = coef
        s2 = coef
        @goto _line210
    end

    #
    # MILLER ALGORITHM FOR CABS(Z).GT.R1
    #
    ak = abs(cos(pi * dnu))
    if ak == 0.0
        #= 300 =#
        #
        # FNU=HALF ODD INTEGER CASE, DNU=-0.5
        #
        s1 = coef
        s2 = coef
        @goto _line210
    end

    fhs = abs(0.25 - dnu2)
    if fhs == 0.0
        #= 300 =#
        #
        # FNU=HALF ODD INTEGER CASE, DNU=-0.5
        #
        s1 = coef
        s2 = coef
        @goto _line210
    end

    #
    # COMPUTE R2=F(E). IF ABS(Z) >= R2, USE FORWARD RECURRENCE TO
    # DETERMINE THE BACKWARD INDEX K. R2=F(E) IS A STRAIGHT LINE ON
    # 12 <= E <= 60. E IS COMPUTED FROM 2**(-E)=B**(1-DIGITS(0.0_dp))=
    # TOL WHERE B IS THE BASE OF THE ARITHMETIC.
    #
    t1 = I1_MACH[14] - 1
    t1 *= D1_MACH[5] * (log(10) / log(2))
    t1 = min(max(t1, 12.0), 60.0)
    t2 = tth * t1 - 6.0
    if real(z) == 0.0
        t1 = hpi
    else #= 130 =#
        t1 = abs(atan(imag(z) / real(z)))
    end
    #= 140 =#
    if t2 <= caz
        #
        # FORWARD RECURRENCE LOOP WHEN CABS(Z).GE.R2
        #
        etest = ak / (pi * caz * tol)
        fk = 1.0
        if etest < 1.0
            @goto _line180
        end
        
        fks = 2.0
        rk = caz + caz + 2.0
        a1 = 0.0
        a2 = 1.0
        for i = 1:kmax
            ak = fhs / fks
            bk = rk / (fk + 1.0)
            tm = a2
            a2 = bk * a2 - ak * a1
            a1 = tm
            rk += 2.0
            fks += fk + fk + 2.0
            fhs += fk + fk
            fk += 1.0
            tm = abs(a2) * fk
            if etest < tm
                break
            end
            if i == kmax
                #= 310 =#
                return -2
            end
        end #= 150 =#

        #= 160 =#
        fk += spi * t1 * sqrt(t2 / caz)
        fhs = abs(0.25 - dnu2)
    else
        #= 170 =#
        #
        # COMPUTE BACKWARD INDEX K FOR CABS(Z).LT.R2
        #
        a2 = sqrt(caz)
        ak *= fpi / (tol * sqrt(a2))
        aa = 3.0 * t1 / (1.0 + caz)
        bb = 14.7 * t1 / (28.0 + caz)
        ak = (log(ak) + caz * cos(aa) / (1.0 + 0.008 * caz)) / cos(bb)
        fk = 0.12125 * ak * ak / caz + 1.5
    end

@label _line180
    #= 180 =#
    #
    # BACKWARD RECURRENCE LOOP FOR MILLER ALGORITHM
    #
    k = Int(trunc(fk))
    fk = float(k)
    fks = fk * fk
    p1 = 0.0 + 0.0im
    p2 = complex(tol)
    cs = p2
    for i = 1:k
        a1 = fks - fk
        a2 = (fks + fk) / (a1 + fhs)
        rk = 2.0 / (fk + 1.0)
        t1 = (fk + real(z)) * rk
        t2 = imag(z) * rk
        pt = p2
        p2 = (p2 * complex(t1, t2) - p1) * a2
        p1 = pt
        cs += p2
        fks = a1 - fk + 1.0
        fk -= 1.0
    end
    #= 190 =#

    #
    # COMPUTE (P2/CS)=(P2/CABS(CS))*(CONJG(CS)/CABS(CS))
    #   FOR BETTER SCALING
    #
    tm = abs(cs)
    pt = 1.0 / tm
    s1 = pt * p2
    cs = conj(cs) * pt
    s1 *= coef * cs
    if (inu <= 0) && (n <= 1)
        zd = z
        if iflag == 1
            @goto _line270
        end
        @goto _line240
    end
    #= 200 =#

    #
    # COMPUTE P1/P2=(P1/CABS(P2)*CONJG(P2)/CABS(P2) FOR SCALING
    #
    tm = abs(p2)
    pt = 1.0 / tm
    p1 = pt * p1
    p2 = conj(p2) * pt
    pt = p1 * p2
    s2 = s1 * (1.0 + (dnu + 0.5 - pt) / z)

@label _line210
    #= 210 =#
    #
    # FORWARD RECURSION ON THE THREE TERM RECURSION WITH RELATION WITH
    # SCALING NEAR EXPONENT EXTREMES ON KFLAG=1 OR KFLAG=3
    #
    ck = (dnu + 1.0) * rz
    if n == 1
        inu -= 1
    end
    if inu <= 0
        if n <= 1
            s1 = s2
        end
        #= 215 =#
        zd = z
        if iflag == 1
            @goto _line270
        end
        @goto _line240
    end

    #= 220 =#
    inub = 1
    if iflag == 1
        @goto _line261
    end

@label _line225
    #= 225 =#
    p1 = csr[kflag]
    ascle = bry[kflag]
    for i = inub:inu
        st = s2
        s2 = ck * s2 + s1
        s1 = st
        ck += rz
        if kflag < 3
            p2 = s2 * p1
            p2m = max(abs(real(p2)), abs(imag(p2)))
            if p2m > ascle
                kflag += 1
                ascle = bry[kflag]
                s1 *= p1
                s2 = p2
                s1 *= css[kflag]
                s2 *= css[kflag]
                p1 = csr[kflag]
            end
        end
    end
    #= 230 =#
    if n == 1
        s1 = s2
    end

@label _line240
    #= 240 =#
    y[1] = s1 * csr[kflag]
    if n == 1
        return nz
    end

    y[2] = s2 * csr[kflag]
    if n == 2
        return nz
    end

    kk = 2
@label _line250
    #= 250 =#
    kk += 1
    if kk > n
        return nz
    end

    p1 = csr[kflag]
    ascle = bry[kflag]
    for i = kk:n
        p2 = s2
        s2 = ck * s2 + s1
        s1 = p2
        ck += rz
        p2 = s2 * p1
        y[i] = p2
        if kflag < 3
            p2m = max(abs(real(p2)), abs(imag(p2)))
            if p2m > ascle
                kflag += 1
                ascle = bry[kflag]
                s1 *= p1
                s2 = p2
                s1 *= css[kflag]
                s2 *= css[kflag]
                p1 = csr[kflag]
            end
        end
    end
    #= 260 =#
    return nz

@label _line261
    #= 261 =#
    #
    # IFLAG=1 CASES, FORWARD RECURRENCE ON SCALED VALUES ON UNDERFLOW
    #
    elm = exp(-elim)
    ascle = bry[1]
    zd = z
    xd = real(z)
    yd = imag(z)
    ic = -1
    j = 2
    for i = 1:inu
        st = s2
        s2 = ck * s2 + s1
        s1 = st
        ck += rz
        as = abs(s2)
        alas = log(as)
        p2r = alas - xd
        if p2r >= -elim
            p2 = -zd + clog(s2)
            p2r = real(p2)
            p2i = imag(p2)
            p2m = exp(p2r) / tol
            p1 = p2m * complex(cos(p2i), sin(p2i))
            if uchk(p1, ascle, tol)
                j = 3 - j
                cy[j] = p1
                if ic == i - 1
                    @goto _line264
                end
                ic = i
                continue
            end
            if alas >= 0.5 * elim
                xd -= elim;
                s1 *= elm;
                s2 *= elm;
                zd = complex(xd, yd);
            end
        end
        #= 262 =#
    end
    if n == 1
        s1 = s2
    end
    @goto _line270
    
@label _line264
    #= 264 =#
    kflag = 1;
    inub = i + 1;
    s2 = cy[j];
    j = 3 - j;
    s1 = cy[j];
    if inub <= inu
        @goto _line225
    end
    if n == 1
        s1 = s2;
    end
    @goto _line240;

@label _line270
    #= 270 =#
    y[1] = s1;
    if n == 1
        y[2] = s2
    end
    #= 280 =#
    ascle = bry[1];
    nz = kscl!(y, zd, fnu, n, rz, ascle, tol, elim);
    inu = n - nz;
    if inu <= 0
        return nz;
    end

    kk = nz + 1;
    s1 = y[kk-1];
    y[kk-1] = s1 * csr[0];
    if inu == 1
        return nz;
    end

    kk = nz + 2;
    s2 = y[kk-1];
    y[kk-1] = s2 * csr[0];
    if inu == 2
        return nz;
    end

    t2 = fnu + (kk-1);
    ck = t2 * rz;
    kflag = 1;
    @goto _line250;
end
