# SPDX-License-Identifier: BSD-3-Clause OR MIT

"""
    s1s2!(
        zr::ComplexF64,
        s1_ref::Ref{ComplexF64},
        s2_ref::Ref{ComplexF64},
        ascle::Float64,
        alim::Float64,
        iuf::Int
    )

# fortran comments
Zs1s2 tests for a possible underflow resulting from the
addition of the I and k functions in the analytic con-
tinuation formula where s1=k function and s2=i function.
On kode=1 the I and k functions are different orders of
magnitude, but for kode=2 they can be of the same order
of magnitude and the maximum must be at least one
precision above the underflow limit.

# Impl Ref
- [`openspecfun/amos/zs1s2.f`](https://github.com/JuliaMath/openspecfun/blob/v0.5.6/amos/zs1s2.f)
- [`scipy/scipy/special/_amos.c:amos_s1s2`](https://github.com/scipy/scipy/blob/b882f1b7ebe55e534f29a8d68a54e4ecd30aeb1a/scipy/special/_amos.c#L4348-L4402)
"""
function s1s2!(
    zr::ComplexF64,
    s1_ref::Ref{ComplexF64},
    s2_ref::Ref{ComplexF64},
    ascle::Float64,
    alim::Float64,
    iuf::Int
)
    s1 = s1_ref[]
    s2 = s2_ref[]

    nz::Int = 0
    as1 = abs(s1)
    as2 = abs(s2)

    if (real(s1) != 0.0 || imag(s1) != 0.0)
        if as1 != 0.0
            xx = real(zr)
            aln = -xx - xx + log(as1)
            s1d = s1
            s1_ref[] = 0.0
            as1 = 0.0
            if aln >= -alim
                c1 = log(s1d)
                c1 -= zr
                c1 -= zr
                s1_ref[] = exp(c1)
                as1 = abs(s1_ref[])
                iuf += 1
            end
        end
    end #= 10 =#

    aa = max(as1, as2)
    if aa > ascle
        return nz, iuf
    end

    s1_ref[] = 0.0
    s2_ref[] = 0.0
    nz = 1
    iuf = 0
    return nz, iuf
end
