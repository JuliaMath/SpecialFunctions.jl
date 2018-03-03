# References:
#  [1] Abramowitz, Stegun: Handbook of Mathematical Functions (1965)
#  [2] ellipjc.m in Toby Driscoll's Schwarz-Christoffel Toolbox


#------------------------------------------------
# Descending and ascending Landen Transformation


descstep(m) = m/(1+sqrt(1-m))^2

@generated function shrinkm(m,::Val{N}) where {N}
    # [1, Sec 16.12]
    quote
        f = one(m)
        Base.Cartesian.@nexprs $N i->begin
            k_i = descstep(m)
            m = k_i^2
            f *= 1+k_i
        end
        return (Base.Cartesian.@ntuple $N k), f, m
    end
end

function ellipj_smallm(u,m)
    # [1, Sec 16.13]
    sn = sin(u) - m*(u-sin(u)*cos(u))*cos(u)/4
    cn = cos(u) + m*(u-sin(u)*cos(u))*sin(u)/4
    dn = 1 - m*sin(u)^2/2;
    return sn,cn,dn
end
function ellipj_largem(u,m1)
    # [1, Sec 16.15]
    sn = tanh(u) + m1*(sinh(u)*cosh(u)-u)*sech(u)^2/4
    cn = sech(u) - m1*(sinh(u)*cosh(u)-u)*tanh(u)*sech(u)/4
    dn = sech(u) + m1*(sinh(u)*cosh(u)+u)*tanh(u)*sech(u)/4
    return sn,cn,dn
end

@generated function ellipj_growm(sn,cn,dn, k::NTuple{N,<:Any}) where {N}
    # [1, Sec 16.12]
    quote
        Base.Cartesian.@nexprs $N i->begin
            kk = k[end-i+1]
            sn,cn,dn = (1+kk)*sn/(1+kk*sn^2),
                       cn*dn/(1+kk*sn^2),
                       (1-kk*sn^2)/(1+kk*sn^2)
                       # ^ Use [1, 16.9.1]. Idea taken from [2]
        end
        return sn,cn,dn
    end
end
@generated function ellipj_shrinkm(sn,cn,dn, k::NTuple{N,<:Any}) where {N}
    # [1, Sec 16.14]
    quote
        Base.Cartesian.@nexprs $N i->begin
            kk = k[end-i+1]
            sn,cn,dn = (1+kk)*sn*cn/dn,
                       (cn^2-kk*sn^2)/dn, # Use [1, 16.9.1]
                       (cn^2+kk*sn^2)/dn  # Use [1, 16.9.1]
        end
        return sn,cn,dn
    end
end

function ellipj_viasmallm(u,m,::Val{N}) where {N}
    k,f,m = shrinkm(m,Val{N}())
    sn,cn,dn = ellipj_smallm(u/f,m)
    return ellipj_growm(sn,cn,dn,k)
end
function ellipj_vialargem(u,m,::Val{N}) where {N}
    k,f,m1 = shrinkm(1-m,Val{N}())
    sn,cn,dn = ellipj_largem(u/f,m1)
    return ellipj_shrinkm(sn,cn,dn,k)
end


#----------------
# Pick algorithm

Base.@pure function nsteps(m,ε)
    i = 0
    while abs(m) > ε
        m = descstep(m)^2
        i += 1
    end
    return i
end
@generated function ellipj_dispatch(u,m)
    T = real(m)
    ε = sqrt(eps(T))
    N = max(
        nsteps(one(T)/2,ε),
        nsteps(sqrt(Complex{T}(0,1)),ε)
    )
    quote
        if abs(m) <= 1 && real(m) <= 0.5
            return ellipj_viasmallm(u,m, Val{$N}())
        elseif abs(1-m) <= 1
            return ellipj_vialargem(u,m, Val{$N}())
        elseif imag(m) == 0 && real(m) < 0
            # [1, Sec 16.10]
            sn,cn,dn = ellipj_dispatch(u*sqrt(1-m),-m/(1-m))
            return sn/(dn*sqrt(1-m)), cn/dn, 1/dn
        else
            # [1, Sec 16.11]
            sn,cn,dn = ellipj_dispatch(u*sqrt(m),1/m)
            return sn/sqrt(m), dn, cn
        end
    end
end


#-----------------------------------
# Type promotion and special values

function ellipj_check(u,m)
    if isfinite(u) && isfinite(m)
        return ellipj_dispatch(u,m)
    else
        T = promote_type(typeof(u),typeof(m))
        return (T(NaN),T(NaN),T(NaN))
    end
end

ellipj(u::Real,m::Real) = ellipj_check(promote(float(u),float(m))...)
ellipj(u::Complex,m::Complex) = ellipj_check(promote(float(u),float(m))...)
function ellipj(u::Complex,m::Real)
    T = promote_type(real.(typeof.(float.((u,m))))...)
    return ellipj_check(convert(Complex{T},u), convert(T,m))
end
function ellipj(u::Real,m::Complex)
    T = promote_type(real.(typeof.(float.((u,m))))...)
    return ellipj_check(convert(T,u), convert(Complex{T},m))
end
