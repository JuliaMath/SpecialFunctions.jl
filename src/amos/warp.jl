# SPDX-License-Identifier: MIT
# This file wraps the exported functions of `libopenspecfun`
#   and is not part of AMOS or its C translation.
using OpenSpecFun_jll

function _uchk(y::ComplexF64, ascle::Float64, tol::Float64)
    nz = Ref{Int32}()

    ccall((:zuchk_, libopenspecfun), Cvoid,
          (Ref{Float64}, Ref{Float64}, Ref{Int32},
           Ref{Float64}, Ref{Float64}),
           real(y), imag(y), nz,
           ascle, tol)

    nz[]
end

function _gammaln(z::Float64)
    isnan(z) && return NaN

    ierr = Ref{Int32}()
    gamln = ccall((:dgamln_, libopenspecfun), Cdouble,
                (Ref{Float64}, Ref{Int32}),
                 z, ierr)

    if 0==ierr[]
        gamln
    else
        NaN
    end
end
