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

end
