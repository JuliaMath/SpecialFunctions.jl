# SPDX-License-Identifier: BSD-3-Clause OR MIT

"""
    kscl!(
        y::Vector{ComplexF64},
        zr::ComplexF64,
        fnu::Float64,
        n::Float64,
        rz::ComplexF64,
        ascle::Float64,
        tol::Float64,
        elim::Float64,
    )

# fortran comments
Set k functions to zero on underflow,
continue recurrence on scaled functions until two members come on scale,
then return with `min(nz+2,n)` values scaled by `1/tol`.

# Impl Ref
- [`openspecfun/amos/zkscl.f`](https://github.com/JuliaMath/openspecfun/blob/v0.5.6/amos/zkscl.f)
- [`scipy/scipy/special/_amos.c:amos_kscl`](https://github.com/scipy/scipy/blob/b882f1b7ebe55e534f29a8d68a54e4ecd30aeb1a/scipy/special/_amos.c#L3949-L4060)
"""
function kscl!(
    y::Vector{ComplexF64},
    zr::ComplexF64,
    fnu::Float64,
    n::Int,
    rz::ComplexF64,
    ascle::Float64,
    tol::Float64,
    elim::Float64,
)
    cy = ComplexF64[0.0im, 0.0im]
    nz::Int = 0
    ic::Int = 0
    nn::Int = ( n > 2 ? 2 : n )
    @assert length(y) >= nn

    kk::Int = 0
    elm = exp(-elim)
    zdr = real(zr)

    for i = 1:nn
        s1 = y[i]
        cy[i] = s1
        as = abs(s1)
        acs = -real(zr) + log(as)
        nz += 1
        y[i] = 0.0im
        if acs < -elim
            continue
        end
        cs = log(s1)
        cs -= zr
        str = exp(real(cs)) / tol
        cs = str * complex(cos(imag(cs)), sin(imag(cs)))
        if 0==uchk(cs, ascle, tol)
            y[i] = cs
            nz -= 1
            ic = i
        end
    end #= 10 =#

    if n == 1
        # @info "ret n==1: " nz
        return nz
    end

    if ic <= 1
        @assert n >= 1
        # TODO: may throw BoundsError here
        y[1] = 0.0im
        nz = 2
    end #= 20 =#
    if n == 2
        # @info "ret n==2: " nz
        return nz
    end
    if nz == 0
        # @info "ret nz==0: " n
        return nz
    end

    fn = fnu + 1.0
    ck = fn * rz
    s1 = cy[1]
    s2 = cy[2]
    zri = imag(zr)
    zd = zr

#   FIND TWO CONSECUTIVE Y VALUES ON SCALE. SCALE RECURRENCE IF
#   S2 GETS LARGER THAN EXP(ELIM/2)
    for i = 3:n
        kk = i
        cs = s2
        s2 *= ck
        s2 += s1
        s1 = cs
        ck += rz
        as = abs(s2)
        alas = log(as)
        acs = alas - zdr
        nz += 1
        y[i] = 0.0im
        if acs >= -elim
            cs = log(s2)
            cs -= zd
            str = exp(real(cs)) / tol
            cs = str * complex(cos(imag(cs)), sin(imag(cs)))
            if 0==uchk(cs, ascle, tol)
                y[i] = cs
                nz -= 1
                if ic == (kk - 1)
                    #= 40 =#
                    nz = kk - 2
                    #= 45 =#
                    for i in 1:nz
                        y[i] = 0.0im
                    end
                    # @info "ret i==(kk-1): " n nz
                    return nz
                end

                ic = kk
                continue
            end
        end #= 25 =#

        if alas >= (0.5 * elim)
            zdr -= elim
            zd = complex(zdr, zri)
            s1 *= elm
            s2 *= elm
        end
    end #= 30 =#

    nz = n
    if ic == n
        nz = n - 1
    end #= 45 =#

    for i in 1:nz
        y[i] = 0.0im
    end
    # @info "ret end kscl!: " n nz
    return nz
end
