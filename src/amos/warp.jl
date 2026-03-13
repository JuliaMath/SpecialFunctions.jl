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

"""
`_s1s2!` will modify `s1_ref` and `s2_ref`.

return (nz, iuf)
"""
function _s1s2!(
    zr::ComplexF64,
    s1_ref::Ref{ComplexF64},
    s2_ref::Ref{ComplexF64},
    ascle::Float64,
    alim::Float64,
    iuf::Int
)
    s1r = Ref{Float64}(real(s1_ref[]))
    s1i = Ref{Float64}(imag(s1_ref[]))
    s2r = Ref{Float64}(real(s2_ref[]))
    s2i = Ref{Float64}(imag(s2_ref[]))
    nz = Ref{Int32}()
    iuf_ref = Ref{Int32}(iuf)

    ccall((:zs1s2_, libopenspecfun), Cvoid,
            (Ref{Float64}, Ref{Float64},  # (ZRR,ZRI) 
             Ref{Float64}, Ref{Float64},  # (S1R,S1I)
             Ref{Float64}, Ref{Float64},  # (S2R,S2I)
            #    NZ,         ASCLE,        ALIM,         IUF
             Ref{Int32}, Ref{Float64}, Ref{Float64}, Ref{Int32}
            ),
             real(zr), imag(zr), 
             s1r, s1i,
             s2r, s2i,
             nz, ascle, alim, iuf_ref)
    
    s1_ref[] = complex(s1r[], s1i[])
    s2_ref[] = complex(s2r[], s2i[])
    nz[], Int(iuf_ref[])
end

"""
Will modify `y[]`.

Return `nz`.
"""
function _asyi!(
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
    y_len = length(y)
    yr = real.(y)
    yi = imag.(y)
    nz = Ref{Int32}()

    ccall((:zasyi_, libopenspecfun),
        Cvoid,
        (
            # (ZR,ZI), 
            Ref{Float64}, Ref{Float64}, 
            # FNU, KODE, N, 
            Ref{Float64}, Ref{Int32}, Ref{Int32},
            # (YR,YI), 
            Ptr{Float64}, Ptr{Float64},
            # NZ, RL, TOL, ELIM, ALIM
            Ref{Int32}, Ref{Float64}, Ref{Float64},
            Ref{Float64}, Ref{Float64}
        ),
        real(z), imag(z),
        fnu, kode, n,
        yr, yi,
        nz, rl, tol, elim, alim
    )
    
    # Update `y`
    y_new = complex.(yr, yi)
    update_arr!(y, y_new)
    @assert length(y) == y_len
    
    nz[]
end


"""
Will modify `y[]`.

Return `nz`.
"""
function _mlri!(
    y::Vector{ComplexF64},
    z::ComplexF64,
    fnu::Float64,
    kode::Int,
    n::Int,
    tol::Float64,
)
    y_len = length(y)
    yr = real.(y)
    yi = imag.(y)
    nz = Ref{Int32}()

    ccall((:zmlri_, libopenspecfun),
        Cvoid,
        (
            # (ZR,ZI),
            Ref{Float64}, Ref{Float64}, 
            # FNU, KODE, N,
            Ref{Float64}, Ref{Int32}, Ref{Int32},
            # (YR,YI),
            Ptr{Float64}, Ptr{Float64},
            # NZ, TOL
            Ref{Int32}, Ref{Float64},
        ),
        real(z), imag(z),
        fnu, kode, n,
        yr, yi,
        nz, tol,
    )

    # Update `y`
    y_new = complex.(yr, yi)
    update_arr!(y, y_new)
    @assert length(y) == y_len
    
    nz[]
end


"""
Will modify `y[]`.

Return `nz`.
"""
function _seri!(
    y::Vector{ComplexF64},
    z::ComplexF64,
    fnu::Float64,
    kode::Int,
    n::Int,
    tol::Float64,
    elim::Float64, 
    alim::Float64
)
    y_len = length(y)
    yr = real.(y)
    yi = imag.(y)
    nz = Ref{Int32}()

    ccall((:zseri_, libopenspecfun),
        Cvoid,
        (
            # (ZR,ZI),
            Ref{Float64}, Ref{Float64}, 
            # FNU, KODE, N,
            Ref{Float64}, Ref{Int32}, Ref{Int32},
            # (YR,YI),
            Ptr{Float64}, Ptr{Float64},
            # NZ, TOL, ELIM, ALIM
            Ref{Int32}, Ref{Float64}, Ref{Float64}, Ref{Float64},
        ),
        real(z), imag(z),
        fnu, kode, n,
        yr, yi,
        nz, tol, elim, alim
    )

    # Update `y`
    y_new = complex.(yr, yi)
    update_arr!(y, y_new)
    @assert length(y) == y_len
    
    nz[]
end


"""
Will modify `y[]`.

Return `nz`.
"""
function _bknu!(
    y::Vector{ComplexF64},
    z::ComplexF64,
    fnu::Float64,
    kode::Int,
    n::Int,
    tol::Float64,
    elim::Float64, 
    alim::Float64
)
    y_len = length(y)
    yr = real.(y)
    yi = imag.(y)
    nz = Ref{Int32}()

    ccall((:zbknu_, libopenspecfun),
        Cvoid,
        (
            # (ZR,ZI),
            Ref{Float64}, Ref{Float64}, 
            # FNU, KODE, N,
            Ref{Float64}, Ref{Int32}, Ref{Int32},
            # (YR,YI),
            Ptr{Float64}, Ptr{Float64},
            # NZ, TOL, ELIM, ALIM
            Ref{Int32}, Ref{Float64}, Ref{Float64}, Ref{Float64},
        ),
        real(z), imag(z),
        fnu, kode, n,
        yr, yi,
        nz, tol, elim, alim
    )

    # Update `y`
    y_new = complex.(yr, yi)
    update_arr!(y, y_new)
    @assert length(y) == y_len
    
    nz[]
end


"""
Will modify `y[]`.

Return `nz`.
"""
function _acai!(
    y::Vector{ComplexF64},
    z::ComplexF64, 
    fnu::Float64, 
    kode::Int, 
    mr::Int, 
    n::Int, 
    rl::Float64, 
    tol::Float64,
    elim::Float64, 
    alim::Float64
)
    y_len = length(y)
    yr = real.(y)
    yi = imag.(y)
    nz = Ref{Int32}()

    ccall((:zacai_, libopenspecfun),
        Cvoid,
        (
            # (ZR,ZI), 
            Ref{Float64}, Ref{Float64}, 
            # FNU, KODE, MR, N, 
            Ref{Float64}, Ref{Int32}, Ref{Int32}, Ref{Int32},
            # (YR,YI), 
            Ptr{Float64}, Ptr{Float64},
            # NZ, RL, TOL, ELIM, ALIM
            Ref{Int32}, Ref{Float64}, Ref{Float64},
            Ref{Float64}, Ref{Float64}
        ),
        real(z), imag(z),
        fnu, kode, mr, n,
        yr, yi,
        nz, rl, tol, elim, alim
    )
    
    # Update `y`
    y_new = complex.(yr, yi)
    update_arr!(y, y_new)
    @assert length(y) == y_len
    
    nz[]
end

#
# `_airy(z::Complex{Float64}, id::Int32, kode::Int32)`
#   defined in SpecialFunctions
#
