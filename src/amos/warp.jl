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

function update_arr!(y::Vector{ComplexF64}, y_new::Vector{ComplexF64})
    @assert size(y) == size(y_new)
    for idx in 1:length(y)
        y[idx] = y_new[idx]
    end
end

function _kscl!(
    y::Vector{ComplexF64},
    zr::ComplexF64,
    fnu::Float64,
    n::Int,
    rz::ComplexF64,
    ascle::Float64,
    tol::Float64,
    elim::Float64,
)
    y_len = length(y)
    nz = Ref{Int32}()
    yr = real.(y)
    yi = imag.(y)
    # @info "[_kscl!] before call" yr, yi
    # @show yr
    # @show yi

    ccall((:zkscl_, libopenspecfun), Cvoid,
            # ZRR,ZRI,FNU,N,
            (Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Int32},
            # YR,YI,NZ,
            Ptr{Float64}, Ptr{Float64}, Ref{Int32},
            # RZR,RZI,
             Ref{Float64}, Ref{Float64},
            # ASCLE,TOL,ELIM
             Ref{Float64}, Ref{Float64}, Ref{Float64},
            ),
             real(zr), imag(zr), fnu, Int32(n),
             yr, yi, nz,
             real(rz), imag(rz),
             ascle, tol, elim)
    # @info "[_kscl!] after call" yr, yi
    # @show yr
    # @show yi

    # Update `y`
    y_new = complex.(yr, yi)
    update_arr!(y, y_new)
    @assert length(y) == y_len

    nz[]
end
