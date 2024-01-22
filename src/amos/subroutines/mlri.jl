# SPDX-License-Identifier: BSD-3-Clause OR MIT

"""

# fortran comments
zmlri computes the i bessel function for re(z).ge.0.0 by the
miller algorithm normalized by a neumann series.

# Impl Ref
- [`openspecfun/amos/zmlri.f`](https://github.com/JuliaMath/openspecfun/blob/v0.5.6/amos/zmlri.f)
- [`scipy/scipy/special/_amos.c:amos_mlri`](https://github.com/scipy/scipy/blob/b882f1b7ebe55e534f29a8d68a54e4ecd30aeb1a/scipy/special/_amos.c#L3778-L3946)
"""
function mlri!(
    y::Vector{ComplexF64},
    z::ComplexF64,
    fnu::Float64,
    kode::Int,
    n::Int,
    tol::Float64,
)
    @assert length(y) >= n

    scle::Float64 = D1_MACH[1] / tol
    az = abs(z)
    iaz = Int(trunc(az))
    ifnu = Int(trunc(fnu))
    inu::Int = ifnu + n - 1
    at = iaz + 1
    ck::ComplexF64 = at / z
    rz = 2.0 / z
    p1 = 0.0 + 0.0im
    p2 = 1.0 + 0.0im
    ack = (at + 1.0) / az
    rho = ack + sqrt(ack*ack - 1.0)
    rho2 = rho * rho
    tst = (rho2 + rho2) / ((rho2 - 1.0) * (rho - 1.0))
    tst /= tol
    
    # COMPUTE RELATIVE TRUNCATION ERROR INDEX FOR SERIES
    ak = at
    i::Int = 1
    converged = false
    for _i = 1:80
        # NOTE: will use i later
        i = _i
        pt = p2
        p2 = p1 - ck * p2
        p1 = pt
        ck += rz
        ap = abs(p2)
        if ap > (tst*ak*ak)
            converged = true
            break
        end
        ak += 1.0
    end #= 10 =#
    if !converged
        return -2
    end

    i += 1
    k::Int = 0
    if inu >= iaz
        # COMPUTE RELATIVE TRUNCATION ERROR FOR RATIOS
        p1 = 0.0 + 0.0im
        p2 = 1.0 + 0.0im
        at = inu + 1
        ck = at / z
        ack = at / az
        tst = sqrt(ack / tol)
        itime = 1
        k = 1
        converged = false
        for _k = 1:80
            # NOTE: will use k later
            k = _k
            pt = p2
            p2 = p1 - ck * p2
            p1 = pt
            ck += rz
            ap = abs(p2)
            if ap >= tst
                if itime == 2
                    converged = true
                    break
                end
                
                ack = abs(ck)
                flam = ack + sqrt(ack*ack - 1.0)
                fkap = ap / abs(p1)
                rho = min(flam, fkap)
                tst *= sqrt(rho / (rho*rho - 1.0))
                itime = 2
            end
        end #= 30 =#
        if !converged
            return -2
        end
    end #= 40 =#

    #
    # BACKWARD RECURRENCE AND SUM NORMALIZING RELATION
    #
    k += 1
    kk = max(i + iaz, k + inu)
    fkk = kk
    p1 = 0.0 + 0.0im
    #
    # SCALE P2 AND SUM BY SCLE
    #
    p2 = complex(scle)
    fnf = fnu - ifnu
    tfnf = fnf + fnf
    bk = gammaln(fkk + tfnf + 1.0) - gammaln(fkk + 1.0) - gammaln(tfnf + 1.0)
    bk = exp(bk)
    sum = 0.0 + 0.0im
    km = kk - inu
    for i = 1:km
        pt = p2
        p2 = p1 + (fkk + fnf) * rz * p2
        p1 = pt
        ak = 1.0 - tfnf / (fkk + tfnf)
        ack = bk * ak
        sum += (ack + bk) * p1
        bk = ack
        fkk -= 1.0
    end #= 50 =#

    y[n] = p2
    if n != 1
        for i = 2:n
            pt = p2
            p2 = p1 + (fkk + fnf) * rz * p2
            p1 = pt
            ak = 1.0 - tfnf / (fkk + tfnf)
            ack = bk * ak
            sum += (ack + bk) * p1
            bk = ack
            fkk -= 1.0
            
            m = n - i + 1
            y[m] = p2
        end
    end #= 60 =#

    #= 70 =#
    if ifnu > 0
        for i = 1:ifnu
            pt = p2
            p2 = p1 + (fkk + fnf) * rz * p2
            p1 = pt
            ak = 1.0 - tfnf / (fkk + tfnf)
            ack = bk * ak
            sum += (ack + bk) * p1
            bk = ack
            fkk -= 1.0
        end
    end #= 80 =#

    #= 90 =#
    pt::ComplexF64 = z
    if kode == 2
        # pt.real = 0.0
        pt -= real(z)
    end
    p1 = -fnf * log(rz) + pt
    ap = gammaln(1.0 + fnf)
    pt = p1 - ap
    # 
    # THE DIVISION EXP(PT)/(SUM+P2) IS ALTERED TO AVOID OVERFLOW
    # IN THE DENOMINATOR BY SQUARING LARGE QUANTITIES
    # 
    p2 += sum
    ap = abs(p2)
    p1 = 1.0 / ap
    ck = exp(pt) * p1
    pt = conj(p2) * p1
    cnorm = ck * pt
    for i = 1:n
        y[i] *= cnorm
    end

    return 0
end
