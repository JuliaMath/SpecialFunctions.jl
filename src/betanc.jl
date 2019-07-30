const errmax = 1e-15

#Compute tail of noncentral Beta distribution
"""
	ncbeta_tail(x,a,b,lambda)

Compute tail of the noncentral beta distribution.
"""
function ncbeta_tail(x::Float64, a::Float64, b::Float64, lambda::Float64)
    if x <= 0.0
        return 0.0
    elseif x >= 1.0
        return 1.0
    end
    
    c = 0.5*lambda
    #Init series

    beta = logabsbeta(a,b)[1]
    temp = beta_inc(a,b,x)[1]
    gx = (SpecialFunctions.beta_integrand(a,b,x,1.0-x))/a
    q = exp(-c)
    xj = 0.0
    ax = q*temp
    sumq = 1.0 - q
    ans = ax

    while true
        xj += 1.0
        temp -= gx
        gx *= x*(a+b+xj-1.0)/(a+xj)
        q *= c/xj
        sumq -= q
        ax = temp*q
        ans += ax

        #Check convergence
        errbd = abs((temp-gx)*sumq)
        if xj > 1000 || errbd < 1e-10
            break
        end
    end
    return ans
end

#R Chattamvelli, R Shanmugam, Algorithm AS 310: Computing the Non-central Beta Distribution #Function, Applied Statistics, Volume 46, Number 1, 1997, pages 146-156

"""
	ncbeta(a,b,lambda,x)

Compute the CDF of the noncentral beta distribution.
"""
function ncbeta(a::Float64, b::Float64, lambda::Float64, x::Float64)
    ans = x
    if x <= 0.0
        return 0.0
    elseif x >= 1.0
        return 1.0
    end
    c = 0.5*lambda
    xj = 0.0

    if lambda < 54.0
        return ncbeta_tail(x,a,b,lambda)
    else
        m = trunc(Int, c+0.5)
        mr = float(m)
        iterlo = m - trunc(Int, 5.0*sqrt(mr))
        iterhi = m + trunc(Int, 5.0*sqrt(mr))
        t = -c + mr*log(c) - logabsgamma(mr + 1.0)[1]
        q = exp(t)
        r = q
        psum = q

        beta = logabsbeta(a+mr,b)[1]
        s1 = (a+mr)*log(x) + b*log(1.0-x) - log(a+mr) - beta
        gx = exp(s1)
        fx = gx
        temp = beta_inc(a+mr,b,x)[1]
        ftemp = temp
        xj += 1.0

        sm = q*temp
        iter1 = m

        #Iterations start from M and goes downwards

        while true
            if iter1 < iterlo || q < errmax
                break
            end

            q *= iter1/c
            xj += 1.0
            gx *= (a + iter1)/(x*(a+b+iter1-1.0))
            iter1 -= 1
            temp += gx
            psum += q
            sm += q*temp
        end

        t0 = logabsgamma(a+b)[1] - logabsgamma(a+1.0)[1] - logabsgamma(b)[1]
        s0 = a*log(x) + b*log(1.0-x)

        s = 0.0
        for j = 0:iter1-1
            s += exp(t0+s0+j*log(x))
            t1 = log(a+b+j) - log(a+j+1.0) + t0
            t0 = t1
        end
        #Compute first part of error bound

        errbd = (gamma_inc(float(iter1),c,0)[2])*(temp+s)
        q = r
        temp = ftemp
        gx = fx
        iter2 = m

        while true
            ebd = errbd + (1.0 - psum)*temp
            if ebd < errmax || iterhi <= iter2
                return sm
            end
            iter2 += 1
            xj += 1.0
            q *= c/iter2
            psum += q
            temp -= gx
            gx *= x*(a+b+iter2-1.0)/(a+iter2)
            sm += q*temp
        end
    end
end

function ncbeta(a::T,b::T,lambda::T,x::T) where {T<:Union{Float16,Float32}}
	T.(ncbeta(Float64(a),Float64(b),Float64(lambda),Float64(x)))
end

ncbeta(a::Real,b::Real,lambda::Real,x::Real) = ncbeta(promote(float(a),float(b),float(lambda),float(x))...)
ncbeta(a::T,b::T,lambda::T,x::T) where {T<:AbstractFloat} = throw(MethodError(ncbeta,(a,b,lambda,x,"")))

