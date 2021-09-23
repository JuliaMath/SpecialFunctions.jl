const errmax = 1e-15

#Compute tail of noncentral Beta distribution
#Russell Lenth, Algorithm AS 226: Computing Noncentral Beta Probabilities,
#Applied Statistics,Volume 36, Number 2, 1987, pages 241-244

"""
	ncbeta_tail(x,a,b,lambda)

Compute tail of the noncentral beta distribution.
Uses the recursive relation
```math
I_{x}(a,b+1;0) = I_{x}(a,b;0) - \\Gamma(a+b)/\\Gamma(a+1)\\Gamma(b)x^{a}(1-x)^{b}
```
and ``\\Gamma(a+1) = a\\Gamma(a)`` given in https://dlmf.nist.gov/8.17.21.
"""
function ncbeta_tail(a::Float64, b::Float64, lambda::Float64, x::Float64)
    if x <= 0.0
        return 0.0
    elseif x >= 1.0
        return 1.0
    end

    c = 0.5*lambda
    #Init series

    beta = logabsbeta(a,b)[1]
    temp = beta_inc(a,b,x)[1]
    gx = (beta_integrand(a,b,x,1.0-x))/a
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

"""
    ncbeta_poisson(a,b,lambda,x)

Compute CDF of noncentral beta if lambda >= 54 using:
First ``\\lambda/2`` is calculated and the Poisson term is calculated using ``P(j-1)=j/\\lambda P(j)`` and ``P(j+1) = \\lambda/(j+1) P(j)``.
Then backward recurrences are used until either the Poisson weights fall below `errmax` or `iterlo` is reached.
```math
I_{x}(a+j-1,b) = I_{x}(a+j,b) + \\Gamma(a+b+j-1)/\\Gamma(a+j)\\Gamma(b)x^{a+j-1}(1-x)^{b}
```
Then forward recurrences are used until error bound falls below `errmax`.
```math
I_{x}(a+j+1,b) = I_{x}(a+j,b) - \\Gamma(a+b+j)/\\Gamma(a+j)\\Gamma(b)x^{a+j}(1-x)^{b}
```
"""
function ncbeta_poisson(a::Float64, b::Float64, lambda::Float64, x::Float64)
    c = 0.5*lambda
    xj = 0.0
    m = round(Int, c)
    mr = float(m)
    iterlo = m - trunc(Int, 5.0*sqrt(mr))
    iterhi = m + trunc(Int, 5.0*sqrt(mr))
    t = -c + mr*log(c) - logabsgamma(mr + 1.0)[1]
    q = exp(t)
    r = q
    psum = q

    beta = logabsbeta(a+mr,b)[1]
    gx = beta_integrand(a+mr,b,x,1.0-x)/(a + mr)
    fx = gx
    temp = beta_inc(a+mr,b,x)[1]
    ftemp = temp
    xj += 1.0

    sm = q*temp
    iter1 = m

    #Iterations start from M and goes downwards

    for iter1 = m:-1:iterlo
        if q < errmax
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
    s0 = a*log(x) + b*log1p(-x)

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
    #Iterations for the higher part

    for iter2 = m:iterhi-1
        ebd = errbd + (1.0 - psum)*temp
        if ebd < errmax
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
    return sm
end

#R Chattamvelli, R Shanmugam, Algorithm AS 310: Computing the Non-central Beta Distribution Function,
#Applied Statistics, Volume 46, Number 1, 1997, pages 146-156

"""
	ncbeta(a,b,lambda,x)

Compute the CDF of the noncentral beta distribution given by
```math
I_{x}(a,b;\\lambda ) = \\sum_{j=0}^{\\infty}q(\\lambda/2,j)I_{x}(a+j,b;0)
```
For ``\\lambda < 54`` : algorithm suggested by Lenth(1987) in `ncbeta_tail(a,b,lambda,x)`.
Else for ``\\lambda >= 54`` : modification in Chattamvelli(1997) in `ncbeta_poisson(a,b,lambda,x)` by using both forward and backward recurrences.
"""
function ncbeta(a::Float64, b::Float64, lambda::Float64, x::Float64)
    ans = x
    if x <= 0.0
        return 0.0
    elseif x >= 1.0
        return 1.0
    end

    if lambda < 54.0
        return ncbeta_tail(a,b,lambda,x)
    else
        return ncbeta_poisson(a,b,lambda,x)
    end
end

"""
    ncF(x,v1,v2,lambda)

Compute CDF of noncentral F distribution given by:
```math
F(x, v1, v2; lambda) = I_{v1*x/(v1*x + v2)}(v1/2, v2/2; \\lambda)
```
where ``I_{x}(a,b; lambda)`` is the noncentral beta function computed above.

Wikipedia: https://en.wikipedia.org/wiki/Noncentral_F-distribution
"""
function ncF(x::Float64, v1::Float64, v2::Float64, lambda::Float64)
    return ncbeta(v1/2, v2/2, lambda, (v1*x)/(v1*x + v2))
end

function ncbeta(a::T,b::T,lambda::T,x::T) where {T<:Union{Float16,Float32}}
	T.(ncbeta(Float64(a),Float64(b),Float64(lambda),Float64(x)))
end

function ncF(x::T,v1::T,v2::T,lambda::T) where {T<:Union{Float16,Float32}}
	T.(ncF(Float64(x),Float64(v1),Float64(v2),Float64(lambda)))
end

ncbeta(a::Real,b::Real,lambda::Real,x::Real) = ncbeta(promote(float(a),float(b),float(lambda),float(x))...)
ncbeta(a::T,b::T,lambda::T,x::T) where {T<:AbstractFloat} = throw(MethodError(ncbeta,(a,b,lambda,x,"")))
ncF(x::Real,v1::Real,v2::Real,lambda::Real) = ncF(promote(float(x),float(v1),float(v2),float(lambda))...)
ncF(x::T,v1::T,v2::T,lambda::T) where {T<:AbstractFloat} = throw(MethodError(ncF,(x,v1,v2,lambda,"")))
