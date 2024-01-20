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
