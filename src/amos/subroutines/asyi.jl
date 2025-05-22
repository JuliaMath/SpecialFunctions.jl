# SPDX-License-Identifier: BSD-3-Clause OR MIT

"""

# fortran comments
zasyi computes the i bessel function for real(z).ge.0.0 by
means of the asymptotic expansion for large cabs(z) in the
region cabs(z).gt.max(rl,fnu*fnu/2). nz=0 is a normal return.
nz.lt.0 indicates an overflow on kode=1.

# Impl Ref
- [`openspecfun/amos/zasyi.f`](https://github.com/JuliaMath/openspecfun/blob/v0.5.6/amos/zasyi.f)
- [`scipy/scipy/special/_amos.c:amos_asyi`](https://github.com/scipy/scipy/blob/b882f1b7ebe55e534f29a8d68a54e4ecd30aeb1a/scipy/special/_amos.c#L997-L1120)
"""
function asyi!(
    y::Vector{ComplexF64},
    z::ComplexF64, 
    fnu::Float64, 
    kode::Int, 
    n::Int, 
    rl::Float64, 
    tol::Float64,
    elim::Float64, 
    alim::Float64
)
    "1 / (2*pi)"
    rpi = INV_2PI

    az = abs(z)
    arm = 1.0e3 * D1_MACH[1]
    rtr1 = sqrt(arm)
    il = min(2, n)
    dfnu = fnu + (n - il)
    
    # OVERFLOW TEST
    ak1 = sqrt(rpi / z)
    cz = z
    if kode == 2
        cz = complex(0.0, imag(z))
    end #= 10 =#

    czr = real(cz)
    if abs(czr) <= elim
        dnu2 = dfnu + dfnu
        koded = 1
        if !(abs(czr) > alim && n > 2)
            koded = 0
            ak1 *= exp(cz)
        end
        #= 20 =#
        fdn = 0.0
        if dnu2 > rtr1
            fdn = dnu2 * dnu2
        end
        ez = z * 8.0
        
        # WHEN Z IS IMAGINARY, THE ERROR TEST MUST BE MADE
        # RELATIVE TO THE FIRST RECIPROCAL POWER SINCE THIS
        # IS THE LEADING TERM OF THE EXPANSION FOR THE
        # IMAGINARY PART.
        aez = 8.0 * az
        s = tol / aez
        jl = Int(trunc(rl + rl)) + 2
        yy = imag(z)
        p1 = 0.0 + 0.0im
        if yy != 0.0
            inu = Int(trunc(fnu))
            arg = (fnu - inu) * pi
            inu += n - il
            ak = -sin(arg)
            bk = cos(arg)
            if yy < 0.0
                bk = -bk
            end
            p1 = complex(ak, bk)
            if isodd(inu)
                p1 = -p1
            end
        end
        #= 30 =#
        for k in 1:il
            sqk = fdn - 1.0
            atol = s * abs(sqk)
            sgn = 1.0
            cs1 = 1.0 + 0.0im
            cs2 = 1.0 + 0.0im
            ck = 1.0 + 0.0im
            ak = 0.0
            aa = 1.0
            bb = aez
            dk = ez
            
            converged = false
            for _ in 1:jl
                ck *= sqk / dk
                cs2 += ck
                sgn = -sgn
                cs1 += ck * sgn
                dk += ez
                aa *= abs(sqk) / bb
                bb += aez
                ak += 8.0
                sqk -= ak
                if aa <= atol
                    converged = true
                    break
                end
            end
            #= 40 =#
            if !converged
                #= 110 =#
                return -2
            end
            
            #= 50 =#
            s2 = cs1
            zr = real(z)
            if (zr + zr) < elim
                s2 += p1 * cs2 * exp(-z - z)
            end
            #= 60 =#
            fdn += 8.0 * dfnu + 4.0
            p1 = -p1
            m = n - il + k
            y[m] = s2 * ak1
        end
        #= 70 =#
        if n <= 2
            return 0
        end
        
        nn = n
        k = nn - 2
        ak = k
        rz = 2.0 / z
        ib = 3
        @assert length(y) >= k
        for _ in ib:nn
            y[k] = (ak + fnu) * rz * y[k+1] + y[k+2]
            ak -= 1.0
            k -= 1
        end
        #= 80 =#
        if koded == 0
            return 0
        end
        
        ck = exp(cz)
        for i in 1:nn
            # ck isa ComplexF64
            y[i] *= ck
        end

        #= 90 =#
        return 0
    end

    #= 100 =#
    return -1
end
