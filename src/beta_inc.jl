using Base.Math: @horner
using Base.MPFR: ROUNDING_MODE
const exparg_n = log(nextfloat(floatmin(Float64)))
const exparg_p =  log(prevfloat(floatmax(Float64)))

#COMPUTE log(gamma(b)/gamma(a+b)) when b >= 8
"""
    loggammadiv(a,b)

Computes ``log(\\Gamma(b)/\\Gamma(a+b))`` when b >= 8
"""
function loggammadiv(a::Float64, b::Float64)
    if a > b
        h = b/a
        c = 1.0/(1.0 + h)
        x = h/(1.0 + h)
        d = a + (b - 0.5)
    else
        h = a/b
        c = h/(1.0 + h)
        x = 1.0/(1.0 + h)
        d = b + a - 0.5
    end
    x² = x*x
    s₃ = 1.0 + (x + x²)
    s₅ = 1.0 + (x + x²*s₃)
    s₇ = 1.0 + (x + x²*s₅)
    s₉ = 1.0 + (x + x²*s₇)
    s₁₁ = 1.0 + (x + x²*s₉)

    # SET W = stirling(b) - stirling(a+b)
    t = inv(b)^2
    w = @horner(t, .833333333333333E-01, -.277777777760991E-02*s₃, .793650666825390E-03*s₅, -.595202931351870E-03*s₇, .837308034031215E-03*s₉, -.165322962780713E-02*s₁₁)
    w *= c/b

    #COMBINING
    u = d*log1p(a/b)
    v = a*(log(b) - 1.0)
    return u <= v ? w - u - v : w - v - u
end

"""
    stirling_corr(a0,b0)

Compute stirling(a0) + stirling(b0) - stirling(a0 + b0)
for a0,b0 >= 8
"""
function stirling_corr(a0::Float64, b0::Float64)
    a = min(a0,b0)
    b = max(a0,b0)

    h = a/b
    c = h/(1.0 + h)
    x = 1.0/(1.0 + h)
    x² = x*x
    #SET SN = (1-X^N)/(1-X)
    s₃ = 1.0 + (x + x²)
    s₅ = 1.0 + (x + x²*s₃)
    s₇ = 1.0 + (x + x²*s₅)
    s₉ = 1.0 + (x + x²*s₇)
    s₁₁ = 1.0 + (x + x²*s₉)
    t = inv(b)^2
    w = @horner(t, .833333333333333E-01, -.277777777760991E-02*s₃, .793650666825390E-03*s₅, -.595202931351870E-03*s₇, .837308034031215E-03*s₉, -.165322962780713E-02*s₁₁)
    w *= c/b
    # COMPUTE stirling(a) + w
    t = inv(a)^2
    return @horner(t, .833333333333333E-01, -.277777777760991E-02, .793650666825390E-03, -.595202931351870E-03, .837308034031215E-03, -.165322962780713E-02)/a + w
end

"""
    esum(mu,x)

Compute ``e^{mu+x}``
"""
function esum(mu::Float64, x::Float64)
    if x > 0.0
        if mu > 0.0 || mu + x < 0.0
            return exp(mu)*exp(x)
        else
            return exp(mu + x)
        end
    elseif mu < 0.0 || mu + x > 0.0
        return exp(mu)*exp(x)
    else
        return exp(mu + x)
    end
end

"""
    beta_integrand(a,b,x,y,mu=0.0)

Compute ``e^{mu} * x^{a}y^{b}/B(a,b)``
"""
function beta_integrand(a::Float64,b::Float64,x::Float64,y::Float64,mu::Float64=0.0)
    a0, b0 = minmax(a,b)
    if a0 >= 8.0
        if a > b
            h = b/a
            x0 = 1.0/(1.0 + h)
            y0 = h/(1.0 + h)
            lambda = (a+b)*y - b
        else
             h = a/b
             x0 = h/(1.0 + h)
             y0 = 1.0/(1.0 + h)
             lambda = a - (a+b)*x
        end
        e = -lambda/a
        u = abs(e) > 0.6 ? u = e - log(x/x0) : - log1pmx(e)
        e = lambda/b
        v = abs(e) > 0.6 ? e - log(y/y0) : - log1pmx(e)
        z = esum(mu, -(a*u + b*v))
        return (1.0/sqrt(2*pi))*sqrt(b*x0)*z*exp(-stirling_corr(a,b))
    elseif x > 0.375
        if y > 0.375
            lnx = log(x)
            lny = log(y)
        else
            lnx = log1p(-y)
            lny = log(y)
        end
    else
        lnx = log(x)
        lny = log1p(-x)
    end
    z = a*lnx + b*lny
    if a0 < 1.0
        b0 = max(a,b)
        if b0 >= 8.0
            u = loggamma1p(a0) + loggammadiv(a0,b0)
            return a0*(esum(mu, z-u))
        elseif b0 > 1.0
            u = loggamma1p(a0)
            n = trunc(Int,b0 - 1.0)
            if n >= 1
                c = 1.0
                for i = 1:n
                    b0 -= 1.0
                    c *= (b0/(a0+b0))
                end
                u += log(c)
            end
            z -= u
            b0 -= 1.0
            apb = a0 + b0
            if apb > 1.0
                u = a0 + b0 - 1.0
                t = (1.0 + rgamma1pm1(u))/apb
            else
                t = 1.0 + rgamma1pm1(apb)
            end
            return a0*(esum(mu,z))*(1.0 + rgamma1pm1(b0))/t
        end
    else
        z -= logbeta(a,b)
        ans = esum(mu,z)
        return ans
    end
    if ans == 0.0
        return 0.0
    end
    apb = a + b
    if apb > 1.0
        z = (1.0 + rgamma1pm1(apb - 1.0))/apb
    else
        z = 1.0 + rgamma1pm1(apb)
    end
    c = (1.0 + rgamma1pm1(a))*(1.0 + rgamma1pm1(b))/z
    return ans*(a0*c)/(1.0 + a0/b0)
end

"""
    beta_inc_cont_fraction(a,b,x,y,lambda,epps)

Compute ``I_{x}(a,b)`` using continued fraction expansion when a,b > 1.
It is assumed that ``\\lambda = (a+b)*y - b``
DLMF : https://dlmf.nist.gov/8.17#E22
BFRAC(A,B,X,Y,LAMBDA,EPS) from Didonato and Morris (1982)
"""
function beta_inc_cont_fraction(a::Float64, b::Float64, x::Float64, y::Float64, lambda::Float64, epps::Float64)
    ans = beta_integrand(a,b,x,y)
    if ans == 0.0
        return 0.0
    end
    c = 1.0 + lambda
    c0 = b/a
    c1 = 1.0 + 1.0/a
    yp1 = y + 1.0

    n = 0.0
    p = 1.0
    s = a + 1.0
    an = 0.0
    bn = 1.0
    anp1 = 1.0
    bnp1 = c/c1
    r = c1/c
    #CONT FRACTION

    while true
     n += 1.0
     t = n/a
     w = n*(b - n)*x
     e = a/s
     alpha = (p*(p+c0)*e*e)*(w*x)
     e = (1.0 + t)/(c1 + 2*t)
     beta = n + w/s +e*(c + n*yp1)
     p = 1.0 + t
     s += 2.0

     #update an, bn, anp1, bnp1
     t = alpha*an  + beta*anp1
     an = anp1
     anp1 = t
     t = alpha*bn + beta*bnp1
     bn = bnp1
     bnp1 = t

     r0 = r
     r = anp1/bnp1
     if abs(r - r0) <= epps*r
        break
     end
     #rescale
     an /= bnp1
     bn /= bnp1
     anp1 = r
     bnp1 = 1.0
    end
    return ans*r
end

"""
    beta_inc_asymptotic_symmetric(a,b,lambda,epps)

Compute ``I_{x}(a,b)`` using asymptotic expansion for a,b >= 15.
It is assumed that ``\\lambda = (a+b)*y - b``
BASYM(A,B,LAMBDA,EPS) from Didonato and Morris (1982)
"""
function beta_inc_asymptotic_symmetric(a::Float64, b::Float64, lambda::Float64, epps::Float64)
    a0 =zeros(22)
    b0 = zeros(22)
    c = zeros(22)
    d = zeros(22)
    e0 = 2/sqrt(pi)
    e1 = 2^(-1.5)
    sm = 0.0
    ans = 0.0
    if a > b
        h = b/a
        r0 = 1.0/(1.0 + h)
        r1 = (b-a)/a
        w0 = 1.0/sqrt(b*(1.0+h))
    else
        h = a/b
        r0 = 1.0/(1.0 + h)
        r1 = (b-a)/b
        w0 = 1.0/sqrt(a*(1.0+h))
    end
    f = -a*log1pmx(-(lambda/a)) - b*log1pmx((lambda/b))
    t = exp(-f)
    if t == 0.0
        return ans
    end
    z0 = sqrt(f)
    z = 0.5*(z0/e1)
    z² = 2.0*f

    a0[1] = (2.0/3.0)*r1
    c[1] = -0.5*a0[1]
    d[1] = - c[1]
    j0 = (0.5/e0)*erfcx(z0)
    j1 = e1
    sm = j0 + d[1]*w0*j1

    s = 1.0
    h² = h*h
    hn = 1.0
    w = w0
    znm1 = z
    zn = z²

    for n = 2: 2: 20
        hn *= h²
        a0[n] = 2.0*r0*(1.0 + h*hn)/(n + 2.0)
        s += hn
        a0[n+1] = 2.0*r1*s/(n+3.0)

        for i = n: n+1
            r = -0.5*(i + 1.0)
            b0[1] = r*a0[1]
            for m = 2:i
                bsum = 0.0
                for j =1: m-1
                    bsum += (j*r - (m-j))*a0[j]*b0[m-j]
                end
                b0[m] = r*a0[m] + bsum/m
            end
            c[i] = b0[i]/(i+1.0)
            dsum = 0.0
            for j = 1: i-1
                imj = i - j
                dsum += d[imj]*c[j]
            end
            d[i] = -(dsum + c[i])
        end

        j0 = e1*znm1 + (n - 1)*j0
        j1 = e1*zn + n*j1
        znm1 *= z²
        zn *= z²
        w *= w0
        t0 = d[n]*w*j0
        w *= w0
        t1 = d[n+1]*w*j1
        sm += (t0 + t1)
        if (abs(t0) + abs(t1)) <= epps*sm
            break
        end
    end

    u = exp(-stirling_corr(a,b))
    return e0*t*u*sm
end

"""
    beta_inc_asymptotic_asymmetric(a,b,x,y,w,epps)

Evaluation of I_{x}(a,b) when b < min(epps,epps*a) and x <= 0.5 using asymptotic expansion.
It is assumed a >= 15 and b <= 1, and epps is tolerance used.
"""
function beta_inc_asymptotic_asymmetric(a::Float64, b::Float64, x::Float64, y::Float64, w::Float64, epps::Float64)
    c = zeros(31)
    d = zeros(31)
    bm1 = b - 1.0
    nu = a + 0.5*bm1
    if y > 0.375
        lnx = log(x)
    else
        lnx = log1p(-y)
    end
    z = -nu*lnx
    if b*z == 0.0
        return error("expansion can't be computed")
    end

    # COMPUTATION OF THE EXPANSION
    #SET R = EXP(-Z)*Z**B/GAMMA(B)
    r = b*(1.0 + rgamma1pm1(b))*exp(b*log(z))
    r *= exp(a*lnx)*exp(0.5*bm1*lnx)
    u = loggammadiv(b,a) + b*log(nu)
    u = r*exp(-u)
    if u == 0.0
        return error("expansion can't be computed")
    end
    (p, q) = gamma_inc(b,z,0)
    v = inv(nu)^2/4
    t2 = lnx^2/4
    l = w/u
    j = q/r
    sm = j
    t = 1.0
    cn = 1.0
    n2 = 0.0
    for n = 1:30
        bp2n = b + n2
        j = (bp2n*(bp2n + 1.0)*j + (z + bp2n + 1.0)*t)*v
        n2 += 2.0
        t *= t2
        cn /= n2*(n2 + 1.0)
        c[n] = cn
        s = 0.0
        if n != 1
            nm1 = n -1
            coef = b - n
            for i = 1:nm1
                s += coef*c[i]*d[n-i]
                coef += b
            end
        end
        d[n] = bm1*cn + s/n
        dj = d[n] * j
        sm += dj
        if sm <= 0.0
            return error("expansion can't be computed")
        end
        if abs(dj) <= epps*(sm+l)
            break
        end
    end
    return w + u*sm
end


#For b < min(eps, eps*a) and x <= 0.5
"""
    beta_inc_power_series2(a,b,x,epps)

Variant of BPSER(A,B,X,EPS).
FPSER(A,B,X,EPS) from Didonato and Morris (1982)
"""
function beta_inc_power_series2(a::Float64, b::Float64, x::Float64, epps::Float64)
    ans = 1.0
    if a > 1.0e-3 * epps
        ans = 0.0
        t = a*log(x)
        if t < exparg_n
            return ans
        end
        ans = exp(t)
    end
    ans *= b/a
    tol = epps/a
    an = a + 1.0
    t = x
    s = t/an
    an += 1.0
    t *= x
    c = t/an
    s += c
    while abs(c) > tol
        an += 1.0
        t *= x
        c = t/an
        s += c
    end
    ans *= (1.0 + a*s)
    return ans
end

#A <= MIN(EPS,EPS*B), B*X <= 1, AND X <= 0.5.,  A is small
"""
    beta_inc_power_series1(a,b,x,epps)

Another variant of BPSER(A,B,X,EPS)
APSER(A,B,X,EPS) from Didonato and Morris (1982)
"""
function beta_inc_power_series1(a::Float64, b::Float64, x::Float64, epps::Float64)
    g = Base.MathConstants.eulergamma
    bx = b*x
    t = x - bx
    if b*epps > 2e-2
        c = log(bx) + g + t
    else
        c = log(x) + digamma(b) + g + t
    end
    tol = 5.0*epps*abs(c)
    j = 1.0
    s = 0.0
    j += 1.0
    t *= (x - bx/j)
    aj = t/j
    s += aj
    while abs(aj) > tol
       j += 1.0
       t *= (x - bx/j)
       aj = t/j
       s += aj
    end
    return -a*(c + s)
end

#B .LE. 1 OR B*X .LE. 0.7
"""
    beta_inc_power_series(a,b,x,epps)

Computes Ix(a,b) using power series :
```math
I_{x}(a,b) = G(a,b)x^{a}/a (1 + a\\sum_{1}^{\\infty}((1-b)(2-b)...(j-b)/j!(a+j)) x^{j})
```
BPSER(A,B,X,EPS) from Didonato and Morris (1982)
"""
function beta_inc_power_series(a::Float64, b::Float64, x::Float64, epps::Float64)
    ans = 0.0
    if x == 0.0
        return 0.0
    end
    a0 = min(a,b)
    b0 = max(a,b)
    if a0 >= 1.0
        z = a*log(x) - logbeta(a,b)
        ans = exp(z)/a
    else

        if b0 >= 8.0
            u = loggamma1p(a0) + loggammadiv(a0,b0)
            z = a*log(x) - u
            ans = (a0/a)*exp(z)
            if ans == 0.0 || a <= 0.1*epps
                return ans
            end
        elseif b0 > 1.0
            u = loggamma1p(a0)
            m = b0 - 1.0
            if m >= 1.0
                c = 1.0
                for i = 1:m
                    b0 -= 1.0
                    c *= (b0/(a0+b0))
                end
                u += log(c)
            end
            z = a*log(x) - u
            b0 -= 1.0
            apb = a0 + b0
            if apb > 1.0
                u = a0 + b0 - 1.0
                t = (1.0 + rgamma1pm1(u))/apb
            else
                t = 1.0 + rgamma1pm1(apb)
            end
            ans = exp(z)*(a0/a)*(1.0 + rgamma1pm1(b0))/t
            if ans == 0.0 || a <= 0.1*epps
                return ans
            end
        else
        #PROCEDURE FOR A0 < 1 && B0 < 1
            ans = x^a
            if ans == 0.0
                return ans
            end
            apb = a + b
            if apb > 1.0
                u = a + b - 1.0
                z = (1.0 + rgamma1pm1(u))/apb
            else
                z = 1.0 + rgamma1pm1(apb)
            end
            c = (1.0 + rgamma1pm1(a))*(1.0 + rgamma1pm1(b))/z
            ans *= c*(b/apb)
            #label l70 start
            if ans == 0.0 || a <= 0.1*epps
                return ans
            end
        end
    end
    if ans == 0.0 || a <= 0.1*epps
        return ans
    end
    # COMPUTE THE SERIES

    sm = 0.0
    n = 0.0
    c = 1.0
    tol = epps/a
    n += 1.0
    c *= x*(1.0 - b/n)
    w = c/(a + n)
    sm += w
    while abs(w) > tol
        n += 1.0
        c *= x*(1.0 - b/n)
        w = c/(a+n)
        sm += w
    end
    return ans*(1.0 + a*sm)
end

"""
    beta_inc_diff(a,b,x,y,n,epps)

Compute ``I_{x}(a,b) - I_{x}(a+n,b)`` where n is positive integer and epps is tolerance.
A more generalised version of https://dlmf.nist.gov/8.17#E20
"""
function beta_inc_diff(a::Float64, b::Float64, x::Float64, y::Float64, n::Integer, epps::Float64)
    apb = a + b
    ap1 = a + 1.0
    mu = 0.0
    d = 1.0
    if n != 1 && a >= 1.0 && apb >= 1.1*ap1
        mu = abs(exparg_n)
        k = exparg_p
        if k < mu
            mu = k
        end
        t = mu
        d = exp(-t)
    end

    ans = beta_integrand(a,b,x,y,mu)/a
    if n == 1 || ans == 0.0
        return ans
    end
    nm1 = n -1
    w = d

    k = 0
    if b <= 1.0
        kp1 = k + 1
        for i = kp1:nm1
            l = i - 1
            d *= ((apb + l)/(ap1 + l))*x
            w += d
            if d <= epps*w
                break
            end
        end
        return ans*w
    elseif y > 1.0e-4
        r = trunc(Int,(b - 1.0)*x/y - a)
        if r < 1.0
            kp1 = k + 1
            for i = kp1:nm1
                l = i - 1
                d *= ((apb + l)/(ap1 + l))*x
                w += d
                if d <= epps*w
                    break
                end
            end
            return ans*w
        end
        k = t = nm1
        if r < t
            k = r
        end
        # ADD INC TERMS OF SERIES
        for i = 1:k
            l = i -1
            d *= ((apb + l)/(ap1 + l))*x
            w += d
        end
        if k == nm1
            return ans*w
        end
    else
        k = nm1
        for i = 1:k
            l = i -1
            d *= ((apb + l)/(ap1 + l))*x
            w += d
        end
        if k == nm1
            return ans*w
        end
    end

    kp1 = k + 1
    for i = kp1:nm1
        l = i - 1
        d *= ((apb + l)/(ap1 + l))*x
        w += d
        if d <= epps*w
            break
        end
   end
   return ans*w
end


#SIGNIFICANT DIGIT COMPUTATION OF INCOMPLETE BETA FUNCTION RATIOS
#by ARMIDO R. DIDONATO AND ALFRED H. MORRIS, JR.
#ACM Transactions on Mathematical Software. Vol18, No3, September1992, Pages360-373
#DLMF : https://dlmf.nist.gov/8.17#E1
#Wikipedia : https://en.wikipedia.org/wiki/Beta_function#Incomplete_beta_function

"""
    beta_inc(a,b,x,y)

Computes Incomplete Beta Function Ratios given by:
```math
I_{x}(a,b) = G(a,b) \\int_{0}^{x} t^{a-1}(1-t)^{b-1} dt,
```
and ``I_{y}(a,b) = 1.0 - I_{x}(a,b)``.
given
``B(a,b) = 1/G(a,b) = \\Gamma(a)\\Gamma(b)/\\Gamma(a+b)`` and ``x+y = 1``.
"""
function beta_inc(a::Float64, b::Float64, x::Float64, y::Float64)
    p = 0.0
    q = 0.0
   # lambda = a - (a+b)*x
    if a < 0.0 || b < 0.0
        return error("a or b is negative")
    elseif a == 0.0 && b == 0.0
        return error("a and b are 0.0")
    elseif x < 0.0 || x > 1.0
        return error("x < 0 or x > 1")
    elseif y < 0.0 || y > 1.0
        return error("y < 0  or y > 1")
    else
        z = x + y - 1.0
        if abs(z) > 3.0*eps()
            return error("x + y != 1.0")         # ERROR HANDLING
        end
    end

    if x == 0.0
        return (0.0,1.0)
    elseif y == 0.0
        return (1.0,0.0)
    elseif a == 0.0
        return (1.0,0.0)
    elseif b == 0.0
        return (0.0,1.0)
    end
#EVALUATION OF ALGOS FOR PROPER SUB-DOMAINS ABOVE
    epps = max(eps(), 1.0e-15)
    if max(a,b) < 1.0E-3 * epps
        return (b/(a+b), a/(a+b))
    end
    ind = false
    a0 = a
    b0 = b
    x0 = x
    y0 = y

    if min(a0,b0) > 1.0
        #PROCEDURE FOR A0>1 AND B0>1
        lambda = a > b ? (a+b)*y - b : a - (a+b)*x
        if lambda < 0.0
            ind = true
            a0 = b
            b0 = a
            x0 = y
            y0 = x
            lambda = abs(lambda)
        end
        if b0 < 40.0 && b0*x0 <= 0.7
            p = beta_inc_power_series(a0,b0,x0,epps)
            q = 1.0 - p
        elseif b0 < 40.0
                n = trunc(Int, b0)
                b0 -= float(n)
                if b0 == 0.0
                    n-=1
                    b0=1.0
                end
                p = beta_inc_diff(b0,a0,y0,x0,n,epps)
                if x0 <= 0.7
                    p += beta_inc_power_series(a0,b0,x0,epps)
                    q = 1.0 - p
                else
                    if a0 <= 15.0
                        n = 20
                        p += beta_inc_diff(a0, b0, x0, y0, n, epps)
                        a0 += n
                    end
                    p = beta_inc_asymptotic_asymmetric(a0,b0,x0,y0,p,15.0*eps())
                    q = 1.0 - p
                end
        elseif a0 > b0
            if b0 <= 100.0 || lambda > 0.03*b0
                p = beta_inc_cont_fraction(a0,b0,x0,y0,lambda,15.0*eps())
                q = 1.0 - p
            else
                p = beta_inc_asymptotic_symmetric(a0,b0,lambda,100.0*eps())
                q = 1.0 - p
            end
        elseif a0 <= 100.0 || lambda > 0.03*a0
            p = beta_inc_cont_fraction(a0,b0,x0,y0,lambda,15.0*eps())
            q = 1.0 - p
        else
            p = beta_inc_asymptotic_symmetric(a0,b0,lambda,100.0*eps())
            q = 1.0 - p
        end
        return ind ? (q, p) : (p, q)
    end
#PROCEDURE FOR A0<=1 OR B0<=1
    if x > 0.5
        ind = true
        a0 = b
        b0 = a
        y0 = x
        x0 = y
    end

    if b0 < min(epps, epps*a0)
        p = beta_inc_power_series2(a0,b0,x0,epps)
        q = 1.0 - p
    elseif a0 < min(epps, epps*b0) && b0*x0 <= 1.0
        q = beta_inc_power_series1(a0,b0,x0,epps)
        p = 1.0 - q
    elseif max(a0,b0) > 1.0
        if b0 <= 1.0
            p = beta_inc_power_series(a0,b0,x0,epps)
            q = 1.0 - p
        elseif x0 >= 0.3
            q = beta_inc_power_series(b0,a0,y0,epps)
            p = 1.0 - q
        elseif x0 >= 0.1
            if b0 > 15.0
                q = beta_inc_asymptotic_asymmetric(b0,a0,y0,x0,q,15.0*eps())
                p = 1.0 - q
            else
                n = 20
                q = beta_inc_diff(b0,a0,y0,x0,n,epps)
                b0 += n
                q = beta_inc_asymptotic_asymmetric(b0,a0,y0,x0,q,15.0*eps())
                p = 1.0 - q
            end
        elseif (x0*b0)^(a0) <= 0.7
            p = beta_inc_power_series(a0,b0,x0,epps)
            q = 1.0 - p
        else
            n = 20
            q = beta_inc_diff(b0,a0,y0,x0,n,epps)
            b0 += n
            q = beta_inc_asymptotic_asymmetric(b0,a0,y0,x0,q,15.0*eps())
            p = 1.0 - q
        end
    elseif a0 >= min(0.2, b0)
        p = beta_inc_power_series(a0,b0,x0,epps)
        q = 1.0 - p
    elseif x0^a0 <= 0.9
        p = beta_inc_power_series(a0,b0,x0,epps)
        q = 1.0 - p
    elseif x0 >= 0.3
        q = beta_inc_power_series(b0,a0,y0,epps)
        p = 1.0 - q
    else
        n = 20
        q = beta_inc_diff(b0,a0,y0,x0,n,epps)
        b0 += n
        q = beta_inc_asymptotic_asymmetric(b0,a0,y0,x0,q,15.0*eps())
        p = 1.0 - q
    end

#TERMINATION
    return ind ? (q, p) : (p, q)
end

beta_inc(a::Float64, b::Float64, x::Float64) = beta_inc(a, b, x, 1.0 - x)
function beta_inc(a::T, b::T, x::T, y::T) where {T<:Union{Float16, Float32}}
    T.(beta_inc(Float64(a), Float64(b), Float64(x), Float64(y)))
end
beta_inc(a::Real, b::Real, x::Real, y::Real) = beta_inc(promote(float(a), float(b), float(x), float(y))...)
beta_inc(a::T, b::T, x::T, y::T) where {T<:AbstractFloat} = throw(MethodError(beta_inc,(a, b, x, y,"")))

#GW Cran, KJ Martin, GE Thomas, Remark AS R19 and Algorithm AS 109: A Remark on Algorithms AS 63: The Incomplete Beta Integral and AS 64: Inverse of the Incomplete Beta Integeral,
#Applied Statistics,
#Volume 26, Number 1, 1977, pages 111-114.

"""
    beta_inc_inv(a,b,p,q,lb=logbeta(a,b)[1])

Computes inverse of incomplete beta function. Given `a`,`b` and ``I_{x}(a,b) = p`` find `x` and return tuple (x,y).
"""
function beta_inc_inv(a::Float64, b::Float64, p::Float64, q::Float64; lb = logbeta(a,b)[1])
    fpu = 1e-30
    x = p
    if p == 0.0
        return (0.0, 1.0)
    elseif p == 1.0
        return (1.0, 0.0)
    end

    #change tail if necessary

    if p > 0.5
        pp = q
        aa = b
        bb = a
        indx = true
    else
        pp = p
        aa = a
        bb = b
        indx = false
    end

    #Initial approx

    r = sqrt(-log(pp^2))
    pp_approx = r - @horner(r, 2.30753e+00, 0.27061e+00) / @horner(r, 1.0, .99229e+00, .04481e+00)

    if a > 1.0 && b > 1.0
        r = (pp_approx^2 - 3.0)/6.0
        s = 1.0/(2*aa - 1.0)
        t = 1.0/(2*bb - 1.0)
        h = 2.0/(s+t)
        w = pp_approx*sqrt(h+r)/h - (t-s)*(r + 5.0/6.0 - 2.0/(3.0*h))
        x = aa/ (aa+bb*exp(w^2))
    else
        r = 2.0*bb
        t = 1.0/(9.0*bb)
        t = r*(1.0-t+pp_approx*sqrt(t))^3
        if t <= 0.0
            x = -expm1((log((1.0-pp)*bb)+lb)/bb)
        else
            t = (4.0*aa+r-2.0)/t
            if t <= 1.0
                x = exp((log(pp*aa)+lb)/aa)
            else
                x = 1.0 - 2.0/(t+1.0)
            end
        end
    end

    #solve x using modified newton-raphson iteration

    r = 1.0 - aa
    t = 1.0 - bb
    pp_approx_prev = 0.0
    sq = 1.0
    prev = 1.0

    if x < 0.0001
        x = 0.0001
    end
    if x > .9999
        x = .9999
    end

    iex = max(-5.0/aa^2 - 1.0/pp^0.2 - 13.0, -30.0)
    acu = 10.0^iex

    #iterate
    while true
        pp_approx = beta_inc(aa,bb,x)[1]
        xin = x
        pp_approx = (pp_approx-pp)*exp(lb+r*log(xin)+t*log1p(-xin))
        if pp_approx * pp_approx_prev <= 0.0
            prev = max(sq, fpu)
        end
        g = 1.0

        tx = x - g*pp_approx
        while true
            adj = g*pp_approx
            sq = adj^2
            tx = x - adj
            if (prev > sq && tx >= 0.0 && tx <= 1.0)
                break
            end
            g /= 3.0
        end

        #check if current estimate is acceptable

        if prev <= acu || pp_approx^2 <= acu
            x = tx
            return indx ? (1.0 - x, x) : (x, 1.0-x)
        end

        if tx == x
            return indx ? (1.0 - x, x) : (x, 1.0-x)
        end

        x = tx
        pp_approx_prev = pp_approx
    end
end

beta_inc_inv(a::Float64, b::Float64, p::Float64) = beta_inc_inv(a, b, p, 1.0-p)
