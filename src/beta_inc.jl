using Base.Math: @horner
using Base.MPFR: ROUNDING_MODE
const exparg_n = -745.133
const exparg_p =  707.5

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
        x² = x*x
        s₃ = 1.0 + (x + x²)
        s₅ = 1.0 + (x + x²*s₃)
        s₇ = 1.0 + (x + x²*s₅)
        s₉ = 1.0 + (x + x²*s₇)
        s₁₁ = 1.0 + (x + x²*s₉)
    else
        h = a/b
        c = h/(1.0 + h)
        x = 1.0/(1.0 + h)
        d = b + a - 0.5
        x² = x*x
        s₃ = 1.0 + (x + x²)
        s₅ = 1.0 + (x + x²*s₃)
        s₇ = 1.0 + (x + x²*s₅)
        s₉ = 1.0 + (x + x²*s₇)
        s₁₁ = 1.0 + (x + x²*s₉)
    end

    # SET W = del(b) - del(a+b)
    t = (1.0/b)^2.0
    w = @horner(t, .833333333333333E-01, -.277777777760991E-02*s₃, .793650666825390E-03*s₅, -.595202931351870E-03*s₇, .837308034031215E-03*s₉, -.165322962780713E-02*s₁₁)
    w *= c/b

    #COMBINING 
    u = d*log1p(a/b)
    v = a*(log(b) - 1.0)
    if u <= v
        return w - u - v
    end
    return w - v - u
end 

"""
    bcorr(a0,b0)

Compute del(a0) + del(b0) - del(a0 + b0) where del(a) is given by:
```math
log(\\Gamma(a)) = (a - 0.5)*log(a) - a + 0.5*log(2\\pi) + del(a)
```
for a0,b0 >= 8
"""
function bcorr(a0::Float64, b0::Float64)
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
    t = (1.0/b)^2.0
    w = @horner(t, .833333333333333E-01, -.277777777760991E-02*s₃, .793650666825390E-03*s₅, -.595202931351870E-03*s₇, .837308034031215E-03*s₉, -.165322962780713E-02*s₁₁)
    w *= c/b
    # COMPUTE DEL(A) + W
    t = inv(a)^2.0
    return @horner(t, .833333333333333E-01, -.277777777760991E-02, .793650666825390E-03, -.595202931351870E-03, .837308034031215E-03, -.165322962780713E-02)/a + w
end

"""
    esum(mu,x)

Compute ``e^{mu+x}``
"""
function esum(mu::Float64, x::Float64)
    if x > 0.0
        if mu > 0.0 || mu + x < 0.0
            w = mu
            return exp(w)*exp(x)
        else
            return exp(mu + x)
        end
    elseif mu < 0.0 || mu + x > 0.0
        w = mu
        return exp(w)*exp(x)
    else
        return exp(mu + x)
    end
end

"""
    brcmp1(mu,a,b,x,y,case)

If case == 1: compute ``e^{mu} * x^{a}y^{b}/B(a,b)``
else : compute above with mu = 0.0
"""
function brcmp1(mu::Float64, a::Float64, b::Float64, x::Float64, y::Float64, case::Bool)
    a0 = min(a,b)
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
        if abs(e) > 0.6
            u = e - log(x/x0)
        else
            u = -SpecialFunctions.log1pmx(e)
        end
        e = lambda/b
        if abs(e) > 0.6
            v = e - log(y/y0)
        else
            v = -SpecialFunctions.log1pmx(e)
        end
        z = case ? esum(mu, -(a*u + b*v)) : exp(-(a*u+b*v))
        return (1.0/sqrt(2*pi))*sqrt(b*x0)*z*exp(-bcorr(a,b))
    elseif x > 0.375
        if y > 0.375
            lnx = log(x)
        else
            lnx = log1p(-y)
        end
        lny = log(y)
    else
        lnx = log(x)
        lny = log1p(-x)
    end
    z = a*lnx + b*lny
    if a0 >= 1.0
        z -= logbeta(a,b)
        return case ? esum(mu,z) : exp(z)
    end

   #PROCEDURE A < 1 OR B < 1
    b0 = max(a,b)
    if b0 >= 8.0
        u = SpecialFunctions.loggamma1p(a0) + loggammadiv(a0,b0)
        return a0*(case ? esum(mu, z-u) : exp(z-u))
    elseif b0 > 1.0
        u = SpecialFunctions.loggamma1p(a0)
        n = b0 - 1.0
        if n >= 1
            c = 1.0
            for i = 1:n
                b0 -= 1.0
                c *= (b0/(a0+b0))
            end
            u += log(c)
        end

    end

     # FOR B0 <= 1
    ans = case ? esum(mu,z) : exp(z)
    if ans == 0.0
        return 0.0
    end
    apb = a + b
    if apb > 1.0
        u = a + b - 1.0
        z = (1.0 + SpecialFunctions.rgamma1pm1(u))/apb
    else
        z = 1.0 + SpecialFunctions.rgamma1pm1(apb)
    end

    c = (1.0 + SpecialFunctions.rgamma1pm1(a))*(1.0 + SpecialFunctions.rgamma1pm1(b))/z
    return ans*(a0*c)/(1.0 + a0/b0)

    z -= u
    b0 -= 1.0
    apb = a0 + b0
    if apb > 1.0
        u = a0 + b0 - 1.0
        t = (1.0 + SpecialFunctions.rgamma1pm1(u))/apb
    else
        t = 1.0 + SpecialFunctions.rgamma1pm1(apb)
    end
    return a0*(case ? esum(mu,z) : exp(z))*(1.0 + SpecialFunctions.rgamma1pm1(b0))/t
end   

"""
    bfrac(a,b,x,y,lambda,epps)

Compute ``I_{x}(a,b)`` using continued fraction expansion when a,b > 1.
It is assumed that ``\\lambda = (a+b)*y - b``
"""
function bfrac(a::Float64, b::Float64, x::Float64, y::Float64, lambda::Float64, epps::Float64)
    ans = brcmp1(0.0,a,b,x,y,false)
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
    basym(a,b,lambda,epps)

Compute ``I_{x}(a,b)`` using asymptotic expansion for a,b >= 15.
It is assumed that ``\\lambda = (a+b)*y - b``
"""
function basym(a::Float64, b::Float64, lambda::Float64, epps::Float64)
    a0 = b0 = c = d = zeros(21)
    e0 = 2/sqrt(pi)
    e1 = 2^(-1.5)

    ans = 0.0
    if a >= b
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
    f = -a*SpecialFunctions.logmxp1(-lambda/a) + b*SpecialFunctions.logmxp1(lambda/b)
    t = exp(-f)
    if t == 0.0
    return ans
    end
    z0 = sqrt(f)
    z = 0.5*(z0/e1)
    z² = 2*f

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

    for n = 2:2:20
        hn *= h²
        a0[n] = 2.0*r0*(1.0 + h*hn)/(n + 2.0)
        np1 = n + 1
        s += hn
        a0[np1] = 2.0*r1*s/(n+3.0)

        for i = n:np1
            r = -0.5*(i + 1.0)
            b0[1] = r*a0[1]
            for m = 2:i
                bsum = 0.0
                mm1 = m - 1
                for j =1:mm1
                    mmj = m - j
                    bsum += (j*r - mmj)*a0[j]*b0[mmj]
                end
                b0[m] = r*a0[m] + bsum/m
            end
            c[i] = b0[i]/(i + 1.0)

            dsum = 0.0
            im1 = i - 1
            for j = 1:im1
                imj = i - j
                dsum += d[imj]*c[j]
            end
            d[i] = -(dsum + c[i])
        end

        j0 = e1*znm1 + (n - 1.0)*j0
        j1 = e1*zn + n*j1
        znm1 *= z²
        zn *= z²
        w *= w0
        t0 = d[n]*w*j0
        w *= w0
        t1 = d[np1]*w*j1
        sm += t0 + t1
        if (abs(t0) + abs(t1)) <= epps*sm
            break
        end
    end
    
    u = exp(-bcorr(a,b))
    return e0*t*u*sm
end        

"""
    bgrat(a,b,x,y,w,epps)

Evaluation of I_{x}(a,b) when b < min(epps,epps*a) and x <= 0.5 using asymptotic expansion.
It is assumed a >= 15 and b <= 1, and epps is tolerance used.
"""
function bgrat(a::Float64, b::Float64, x::Float64, y::Float64, w::Float64, epps::Float64)
    c = zeros(30)
    d = zeros(30)
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
    r = b*(1.0 + SpecialFunctions.rgamma1pm1(b))*exp(b*log(z))
    r *= exp(a*lnx)*exp(0.5*bm1*lnx)
    u = loggammadiv(b,a) + b*log(nu)
    u = r*exp(-u)
    if u == 0.0
        return error("expansion can't be computed")
    end
    (p, q) = gamma_inc(b,z,0)
    v = 0.25*(1.0/nu)^2.0
    t2 = 0.25*lnx*lnx
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
            return w+u*sm
        end
    end
end 



function fpser(a::Float64, b::Float64, x::Float64, epps::Float64)
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

#A .LE. MIN(EPS,EPS*B), B*X .LE. 1, AND X .LE. 0.5.,  A is small

function apser(a::Float64, b::Float64, x::Float64, epps::Float64)         
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

function bpser(a::Float64, b::Float64, x::Float64, epps::Float64)
    ans = 0.0
    if x == 0.0
        return 0.0
    end
    a0 = min(a,b)
    if a0 >= 1.0
        z = a*log(x) - logbeta(a,b)
        ans = exp(z)/a
    end


    b0 = max(a,b)
    if b0 >= 8.0
        u = SpecialFunctions.rgamma1pm1(a0) + loggammadiv(a0,b0)
        z = a*log(x) - u
        ans = (a0/a)*exp(z)
        if ans == 0.0 || a <= 0.1*epps
            return ans
        end
    elseif b0 > 1.0
        u = SpecialFunctions.loggamma1p(a0)
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
            t = (1.0 + SpecialFunctions.rgamma1pm1(u))/apb
        else
            t = 1.0 + SpecialFunctions.rgamma1pm1(apb)
        end
        ans = exp(z)*(a0/a)*(1.0 + SpecialFunctions.rgamma1pm1(b0))/t
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
            z = (1.0 + SpecialFunctions.rgamma1pm1(u))/apb
        else
            z = 1.0 + SpecialFunctions.rgamma1pm1(apb)
        end
        c = (1.0 + SpecialFunctions.rgamma1pm1(a))*(1.0 + SpecialFunctions.rgamma1pm1(b))/z
        ans *= c*(b/apb)
        #label l70 start
        if ans == 0.0 || a <= 0.1*epps
            return ans
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
    bup(a,b,x,y,n,epps)

Compute ``I_{x}(a,b) - I_{x}(a+n,b)`` where n is positive integer and epps is tolerance.
"""
function bup(a::Float64, b::Float64, x::Float64, y::Float64, n::Integer, epps::Float64)
    apb = a + b
    ap1 = a + 1.0
    mu = 0.0
    d = 1.0
    if n != 1 && a >= 1.0
        mu = abs(exparg_n)
        k = exparg_p
        if k < mu 
            mu = k
        end
        t = mu
        d = exp(-t)
    end

    ans = brcmp1(mu,a,b,x,y,true)/a
    if n == 1 || ans == 0.0
        return ans
    end
    nm1 = n -1
    w = d

    k = 0
    if b > 1.0
        if y > 1e-4
            r = (b - 1.0)*x/y - a
            if r < 1.0
                @goto l40
            end
            k = t = nm1
            if r < t
                k = r
            end
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
        end
        for i = 1:k
            l = i -1
            d *= ((apb + l)/(ap1 + l))*x
            w += d
        end
        if k == nm1
            return ans*w
        end
    end
    @label l40
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

"""
    beta_inc(a,b,x,y)

Computes Incomplete Beta Function Ratios given by:
```math
I_{x}(a,b) = G(a,b) \\int_{0}^{x} t^{a-1}(1-t)^{b-1} dt,
```
and ``I_{y}(a,b) = 1.0 - I_{x}(a,b)``.
given
``B(a,b) = 1/G(a,b) = \\Gamma(a)\\Gamma(b)/\\Gamma(a+b)`` and y = 1.0 - x
"""
function beta_inc(a::Float64, b::Float64, x::Float64, y::Float64)
    ans_x = 0.0
    ans_y = 0.0
    lambda = a - (a+b)*x
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
        @goto l30
    end
#PROCEDURE FOR A0<=1 OR B0<=1
    if x <= 0.5
        @goto l10
    end
    ind = true
    a0 = b
    b0 = a
    y0 = x
    x0 = y

    @label l10
     if b0 < min(epps, epps*a0)
        @goto l80
     elseif a0 < min(epps, epps*b0)
        @goto l90
     elseif max(a0,b0) > 1.0
        if b0 <= 1.0 || (x0*b0)^a0 <= 0.7
            @goto l100
        elseif x0 >= 0.3
            @goto l110
        elseif x0 >= 0.1
            if b0 > 15.0
                @goto l131
            else
                n = 20
                @goto l130
            end
        end
     elseif a0 >= min(0.2, b0) || x0^a0 <= 0.9
        @goto l100
     elseif x0 >= 0.3
        @goto l110
     else
        n = 20
        @goto l130
     end

#PROCEDURE FOR A0>1 AND B0>1
    @label l30
     if a > b
        lambda = (a+b)*y - b
     else
        lambda = a - (a+b)*x
     end
     if lambda >= 0.0
        @goto l40
     end
    
    ind = true
    a0 = b
    b0 = a
    x0 = y
    y0 = x

    lambda = abs(lambda)
    @label l40
     if b0 < 40.0 && b0*x0 <= 0.7

        @goto l100
     elseif b0 < 40.0
        @goto l140
     elseif a0 > b0
        if b0 <= 100.0 || lambda > 0.03*b0
            @goto l120
        else
            @goto l180
        end
     elseif a0 <= 100.0 || lambda > 0.03*a0
        @goto l120
     else
        @goto l180
     end

#EVALUATION OF ALGOS FOR PROPER SUB-DOMAINS

    @label l80
     ans_x = fpser(a0,b0,x0,epps)
     ans_y = 1.0 - ans_x
     @goto l220
     
    @label l90
     ans_y = apser(a0,b0,x0,epps)
     ans_x = 1.0 - ans_y
     @goto l220

    @label l100
     ans_x = bpser(a0,b0,x0,epps)
     ans_y = 1.0 - ans_x
     @goto l220

    @label l110
     ans_y = bpser(b0,a0,y0,epps)
     ans_x = 1.0 - ans_y
     @goto l220

    @label l120
     ans_x = bfrac(a0,b0,x0,y0,lambda,15.0*eps())
     ans_y = 1.0 - ans_x
     @goto l220

    @label l130
     ans_y = bup(b0,a0,y0,x0,n,epps)
     b0 += n
    @label l131
     ans_y = bgrat(b0,a0,y0,x0,ans_y,15.0*eps())
     ans_x = 1.0 - ans_y
     @goto l220
    
    @label l140
     n = b0
     b0 -= n
     if b0 != 0.0
        @goto l141
     end
     n-=1
     b0=1.0

    @label l141
     ans_x = bup(b0,a0,y0,x0,n,epps)
     if x0 > 0.7
        @goto l150
     end
     ans_x += bpser(a0,b0,x0,epps)
     ans_y = 1.0 - ans_x
     @goto l220

    @label l150
     if a0 > 15.0
        @goto l151
     end
     n = 20
     ans_x += bup(a0, b0, x0, y0, n, epps)
     a0 += n

    @label l151
     ans_x = bgrat(a0,b0,x0,y0,15.0*eps())
     ans_y = 1.0 - ans_x
     @goto l220

    @label l180
     ans_x = basym(a0,b0,lambda,100.0*eps())
     ans_y = 1.0 - ans_x
     @goto l220

#TERMINATION
    @label l220
     if !ind
        return (ans_x,ans_y)
     end
     #ans_x, ans_y = ans_y, ans_x
     return (ans_y, ans_x)
end
