using Base.Math: @horner
using Base.MPFR: ROUNDING_MODE
#useful constants
const acc0 = [5.0e-15 , 5.0e-7 , 5.0e-4] #accuracy options
const big1 = [25.0 , 14.0, 10.0]
const e0 = [0.25e-3 , 0.25e-1 , 0.14]
const x0=[31.0 , 17.0 , 9.7]
const alog10 = log(10)
const rt2pin = 1.0/sqrt(2*pi)
const rtpi = sqrt(pi) 
const exparg = -745.1


#----------------COEFFICIENTS FOR TEMME EXPANSION------------------

const d00 = -.333333333333333E+00
const d0=[.833333333333333E-01 , -.148148148148148E-01 , .115740740740741E-02 , .352733686067019E-03 , -.178755144032922E-03 , .391926317852244E-04]
const d10 = -.185185185185185E-02
const d1=[-.347222222222222E-02 , .264550264550265E-02 , -.990226337448560E-03 , .205761316872428E-03]
const d20 = .413359788359788E-02
const d2=[-.268132716049383E-02 , .771604938271605E-03]
const d30 = .649434156378601E-03
const d3=[.229472093621399E-03 , -.469189494395256E-03]
const d40 = -.861888290916712E-03
const d4=[.784039221720067E-03]
const d50 = -.336798553366358E-03
const d5=[-.697281375836586E-04]
const d60 = .531307936463992E-03
const d6=[-.592166437353694E-03]
const d70 = .344367606892378E-03
const d80 = -.652623918595309E-03
#Source of logmxp1(x): https://github.com/JuliaStats/StatsFuns.jl/blob/master/src/basicfuns.jl
# The kernel of log1pmx
# Accuracy within ~2ulps for -0.227 < x < 0.315
function _log1pmx_ker(x::Float64)
    r = x/(x+2.0)
    t = r*r
    w = @horner(t,
                6.66666666666666667e-1, # 2/3
                4.00000000000000000e-1, # 2/5
                2.85714285714285714e-1, # 2/7
                2.22222222222222222e-1, # 2/9
                1.81818181818181818e-1, # 2/11
                1.53846153846153846e-1, # 2/13
                1.33333333333333333e-1, # 2/15
                1.17647058823529412e-1) # 2/17
    hxsq = 0.5*x*x
    r*(hxsq+w*t)-hxsq
end
"""
    log1pmx(x::Float64)
Return `log(1 + x) - x`
Use naive calculation or range reduction outside kernel range.  Accurate ~2ulps for all `x`.
"""
function log1pmx(x::Float64)
    if !(-0.7 < x < 0.9)
        return log1p(x) - x
    elseif x > 0.315
        u = (x-0.5)/1.5
        return _log1pmx_ker(u) - 9.45348918918356180e-2 - 0.5*u
    elseif x > -0.227
        return _log1pmx_ker(x)
    elseif x > -0.4
        u = (x+0.25)/0.75
        return _log1pmx_ker(u) - 3.76820724517809274e-2 + 0.25*u
    elseif x > -0.6
        u = (x+0.5)*2.0
        return _log1pmx_ker(u) - 1.93147180559945309e-1 + 0.5*u
    else
        u = (x+0.625)/0.375
        return _log1pmx_ker(u) - 3.55829253011726237e-1 + 0.625*u
    end
end
"""
    logmxp1(x::Float64)
Return `log(x) - x + 1` carefully evaluated.
"""
function logmxp1(x::Float64)
    if x <= 0.3
        return (log(x) + 1.0) - x
    elseif x <= 0.4
        u = (x-0.375)/0.375
        return _log1pmx_ker(u) - 3.55829253011726237e-1 + 0.625*u
    elseif x <= 0.6
        u = 2.0*(x-0.5)
        return _log1pmx_ker(u) - 1.93147180559945309e-1 + 0.5*u
    else
        return log1pmx(x - 1.0)
    end
end

"""
   rgamma1pm1(a)

   Computation of 1/Gamma(a+1) - 1 for -0.5<=a<=1.5 : ``1/\\Gamma (a+1) - 1``
   Uses the relation gamma(a+1) = a*gamma(a)
"""
function rgamma1pm1(a::Float64)
    t=a
    rangereduce = a > 0.5
    t = rangereduce ? a-1 : a #-0.5<= t <= 0.5
    if t == 0.0
        return 0.0
    elseif t < 0.0        
        top = @horner(t , -.422784335098468E+00 , -.771330383816272E+00 , -.244757765222226E+00 , .118378989872749E+00 , .930357293360349E-03 , -.118290993445146E-01 , .223047661158249E-02 , .266505979058923E-03 , -.132674909766242E-03)
        bot = @horner(t , 1.0 , .273076135303957E+00 , .559398236957378E-01)
        w = top/bot
        return rangereduce ? t*w/a : a*(w+1)
    else
        top = @horner(t , .577215664901533E+00 , -.409078193005776E+00 , -.230975380857675E+00 , .597275330452234E-01 , .766968181649490E-02 , -.514889771323592E-02 , .589597428611429E-03)
        bot = @horner(t , 1.0 , .427569613095214E+00 , .158451672430138E+00 , .261132021441447E-01 , .423244297896961E-02) 
        w = top/bot
        return rangereduce ? (t/a)*(w - 1.0) : a*w
    end
end
"""
    rgammax(a,x)

Evaluation of exp(-x)*x^a/gamma(a) : ``1/\\Gamma(a) e^{-x} x^{a}``
"""
function rgammax(a::Float64,x::Float64)
    if x == 0.0
        return 0.0
    elseif a >= 20.0
        u =x/a
        if u == 0.0
            return 0.0
        end
        t = (1.0/a)^2
        t1 = (((0.75*t - 1.0)*t + 3.5)*t - 105.0)/(a*1260.0) #Using Stirling Series : https://dlmf.nist.gov/5.11.1
        t1 = t1 + a*logmxp1(u)
        if t1 >= exparg
            return rt2pin*sqrt(a)*exp(t1)
        end
    else
        t = a*log(x) - x
        if t < exparg
            return 0.0
        elseif a >= 1.0
            return exp(t)/gamma(a)
        else
            return (a*exp(t))*(1.0 + rgamma1pm1(a))
        end
    end
end   
"""
    gamma_p_cf(a, x, ind)

Computes P(a,x) by continued fraction expansion given by : 
```math
R(a,x) * \\frac{1}{1-\\frac{z}{a+1+\\frac{z}{a+2-\\frac{(a+1)z}{a+3+\\frac{2z}{a+4-\\frac{(a+2)z}{a+5+\\frac{3z}{a+6-\\dots}}}}}}}
```
Used when 1 <= a <= BIG and x < x0.
DLMF : https://dlmf.nist.gov/8.9#E2
"""
function gamma_p_cf(a::Float64, x::Float64, ind::Integer)
    acc = acc0[ind + 1]
    tol = 4.0*acc
    a2nm1 = 1.0
    a2n = 1.0
    b2nm1 = x
    b2n = x + (1.0 - a)
    c = 1.0
    while true
       a2nm1 = x*a2n + c*a2nm1
       b2nm1 = x*b2n + c*b2nm1
       c = c + 1.0
       t = c - a
       a2n = a2nm1 + t*a2n
       b2n = b2nm1 + t*b2n
       a2nm1 = a2nm1/b2n
       b2nm1 = b2nm1/b2n
       a2n = a2n/b2n
       b2n = 1.0
       if abs(a2n - a2nm1/b2nm1) < tol*a2n
           break
       end
    end
    q_ans = rgammax(a,x)*a2n
    return (1.0 - q_ans, q_ans)
end
"""
    gamma_p_taylor(a, x, ind)

Compute P(a,x) using Taylor Series for P/R given by : 
```math
R(a,x)/a * (1 + \\sum_{1}^{\\infty} x^{n}/((a+1)(a+2)...(a+n)))
```
Used when 1 <= a <= BIG and x <= max{a, ln 10}.
DLMF : https://dlmf.nist.gov/8.11#E2
"""
function gamma_p_taylor(a::Float64, x::Float64, ind::Integer)
    acc = acc0[ind + 1]
    wk = zeros(30)
    flag = false
    apn = a + 1.0
    t = x/apn
    wk[1] = t
    loop=2
    for indx = 2:20
       apn = apn + 1.0
       t = t*(x/apn)
       if t <= 1.0e-3
           loop = indx
           flag = true
           break
       end
       wk[indx] = t
    end
    if !flag
        loop = 20
    end
    sm = t
    tol = 0.5*acc #tolerance
    while true
       apn = apn+1.0
       t = t*(x/apn)
       sm = sm + t
       if t <= tol
           break
       end
    end
    for j = loop-1:-1:1
       sm += wk[j]
    end
    p_ans = (rgammax(a,x)/a)*(1.0 + sm)
    return (p_ans, 1.0-p_ans)
end
"""
    gamma_p_asym(a, x, ind)

Compute P(a,x) using asymptotic expansion given by : 
```math
R(a,x)/a * (1 + \\sum_{1}^{N-1}(a_{n}/x^{n} + \\Theta _{n}a_{n}/x^{n}))
```
where R(a,x) = rgammax(a,x). Used when 1 <= a <= BIG and x >= x0.
DLMF : https://dlmf.nist.gov/8.11#E2
"""
function gamma_p_asym(a::Float64, x::Float64, ind::Integer)
    wk = zeros(30)
    flag = false
    acc = acc0[ind + 1]
    amn = a-1.0
    t=amn/x
    wk[1]=t
    loop=2
    for indx = 2 : 20
       amn = amn-1.0
       t=t*(amn/x)
       if abs(t) <= 1.0e-3
           loop = indx
           flag = true
           break
       end
       wk[indx]=t
    end
    if !flag
        loop = 20
    end
    sm = t
    while true
       if abs(t) < acc
           break
       end
       amn=amn-1.0
       t=t*(amn/x)
       sm=sm+t
    end
    for j = loop-1:-1:1
       sm += wk[j]
    end
    q_ans = (rgammax(a,x)/x)*(1.0 + sm)
    return (1.0 - q_ans, q_ans) 
end
"""
    gamma_p_taylor_x(a,x,ind)

Computes P(a,x) based on Taylor expansion of P(a,x)/x**a given by:
```math
J = -a * \\sum_{1}^{\\infty} (-x)^{n}/((a+n)n!)
``` and P(a,x)/x**a is given by :
```math
(1 - J)/ \\Gamma(a+1)
``` resulting from term-by-term integration of gamma_p(a,x,ind).
This is used when a < 1 and x < 1.1 - Refer Eqn (9) in the paper.
"""
function gamma_p_taylor_x(a::Float64, x::Float64, ind::Integer)
    acc = acc0[ind + 1]
    l=3.0
    c=x
    sm= x/(a + 3.0)
    tol = 3.0*acc/(a + 1.0)
    while true
       l=l+1.0
       c=-c*(x/l)
       t=c/(a+l)
       sm=sm+t
       if abs(t) <= tol 
           break
       end
    end
    temp = a*x*((sm/6.0 - 0.5/(a + 2.0))*x + 1.0/(a + 1.0))
    z = a*log(x)
    h = rgamma1pm1(a)
    g = 1.0 + h
    if (x < 0.25 && z > -.13394) || a < x/2.59
       l = expm1(z)
       w = 1.0+l
       rangered = ((w*temp - l)*g - h < 0.0)
       return rangered ? (1.0,0.0) : ((1.0 - ((w*temp - l)*g - h)),((w*temp - l)*g - h))
    else
       w = exp(z)
       return (w*g*(1.0 - temp), 1.0 - w*g*(1.0-temp))
    end
end
"""
    gamma_p_minimax(a,x,z)

Compute P(a,x) using minimax approximations given by : 
```math
1/2 * erfc(\\sqrt{y}) - e^{-y}/\\sqrt{2\\pi*a}* T(a,\\lambda)
``` where 
```math
T(a,\\lambda) = \\sum_{0}^{N} c_{k}(z)a^{-k}
```
DLMF : https://dlmf.nist.gov/8.12#E8
This is a higher accuracy approximation of Temme expansion, which deals with the region near a ≈ x with a large.
Refer Appendix F in the paper for the extensive set of coefficients calculated using Brent's multiple precision arithmetic(set at 50 digits) in BRENT, R. P. A FORTRAN multiple-precision arithmetic package, ACM Trans. Math. Softw. 4(1978), 57-70 .
"""
function gamma_p_minimax(a::Float64, x::Float64, z::Float64)
    l = x/a
    s = 1.0 - l
    y = -a*logmxp1(l)
    c = exp(-y)
    w = 0.5*erfcx(sqrt(y))

    if abs(s) <= 1.0e-3
        c0 = @horner(z , d00 , d0[1] , d0[2] , d0[3])
        c1 = @horner(z , d10 , d1[1] , d1[2] , d1[3])
        c2 = @horner(z , d20 , d2[1] , d2[2])
        c3 = @horner(z , d30 , d3[1] , d3[2])
        c4 = @horner(z , d40 , d4[1])
        c5 = @horner(z , d50 , d5[1])
        c6 = @horner(z , d60 , d6[1])

        t = @horner(u , c0 , c1 , c2 , c3 , c4 , c5 , c6 , d70 , d80)
        if l < 1.0
            return (c*(w - rt2pin*t/sqrt(a)) , 1.0 - c*(w - rt2pin*t/sqrt(a)))
        end
        return (1.0 - c*(w + rt2pin*t/sqrt(a)) , c*(w + rt2pin*t/sqrt(a)))
    end
    #---USING THE MINIMAX APPROXIMATIONS---
    c0 = @horner(z , -.333333333333333E+00 , -.159840143443990E+00 , -.335378520024220E-01 , -.231272501940775E-02)/(@horner(z , 1.0 , .729520430331981E+00 , .238549219145773E+00, .376245718289389E-01 , .239521354917408E-02 , -.939001940478355E-05 , .633763414209504E-06) )
    c1 = @horner(z , -.185185185184291E-02 , -.491687131726920E-02 , -.587926036018402E-03 , -.398783924370770E-05)/(@horner(z , 1.0 , .780110511677243E+00 , .283344278023803E+00 , .506042559238939E-01 , .386325038602125E-02))
    c2 = @horner(z , .413359788442192E-02 , .669564126155663E-03)/(@horner(z , 1.0 , .810647620703045E+00 , .339173452092224E+00 , .682034997401259E-01 , .650837693041777E-02 , -.421924263980656E-03))
    c3 = @horner(z , .649434157619770E-03 , .810586158563431E-03)/(@horner(z , 1.0 , .894800593794972E+00, .406288930253881E+00 , .906610359762969E-01 , .905375887385478E-02 , -.632276587352120E-03))
    c4 = @horner(z , -.861888301199388E-03 , -.105014537920131E-03)/(@horner(z , 1.0 , .103151890792185E+01 , .591353097931237E+00 , .178295773562970E+00 , .322609381345173E-01))
    c5 = @horner(z , -.336806989710598E-03, -.435211415445014E-03)/(@horner(z , 1.0 , .108515217314415E+01 , .600380376956324E+00 , .178716720452422E+00))
    c6 = @horner(z , .531279816209452E-03 , -.182503596367782E-03)/(@horner(z , 1.0 , .770341682526774E+00 , .345608222411837E+00))
    c7 = @horner(z , .344430064306926E-03, .443219646726422E-03)/(@horner(z , 1.0 , .115029088777769E+01 , .821824741357866E+00))
    c8 = @horner(z , -.686013280418038E-03 , .878371203603888E-03)

    t = @horner(1.0/a , c0 , c1 , c2 , c3 , c4 , c5 , c6 , c7 , c8)
    if l < 1.0
       return (c*(w - rt2pin*t/sqrt(a)) , 1.0 - c*(w - rt2pin*t/sqrt(a)))
    end
    return (1.0 - c*(w + rt2pin*t/sqrt(a)) , c*(w + rt2pin*t/sqrt(a)))
end
"""
    gamma_p_temme(a, x, z)

Compute P(a,x) using Temme's expansion given by : 
```math
1/2 * erfc(\\sqrt{y}) - e^{-y}/\\sqrt{2\\pi*a}* T(a,\\lambda)
``` where 
```math
T(a,\\lambda) = \\sum_{0}^{N} c_{k}(z)a^{-k}
```
DLMF : https://dlmf.nist.gov/8.12#E8
This mainly solves the problem near the region when a ≈ x with a large, and is a lower accuracy version of the minimax approximation.
"""
function gamma_p_temme(a::Float64, x::Float64, z::Float64)
    l = x/a
    y = -a*logmxp1(x/a)
    c = exp(-y)
    w = 0.5*erfcx(sqrt(y))
    c0 = @horner(z , d00 , d0[1] , d0[2] , d0[3] , d0[4] , d0[5] , d0[6]) 
    c1 = @horner(z , d10 , d1[1] , d1[2] , d1[3] , d1[4]) 
    c2 = @horner(z , d20 , d2[1])
    t = @horner(1.0/a , c0 , c1 , c2)
    if l < 1.0
       return (c*(w - rt2pin*t/sqrt(a)) , 1.0 - c*(w - rt2pin*t/sqrt(a)))
    end
    return (1.0 - c*(w + rt2pin*t/sqrt(a)) , c*(w + rt2pin*t/sqrt(a)))
end
"""
    gamma_p_temme_1(a, x, z, ind)

Computes P(a,x) using simplified Temme expansion near y=0 by : 
```math
E(y) - (1 - y)/\\sqrt{2\\pi*a} * T(a,\\lambda)
``` where 
```math
E(y) = 1/2 - (1 - y/3)*(\\sqrt(y/\\pi))
```
Used instead of it's previous function when ``\\sigma <= e_{0}/\\sqrt{a}``.
DLMF : https://dlmf.nist.gov/8.12#E8
"""
function gamma_p_temme_1(a::Float64, x::Float64, z::Float64, ind::Integer)
    iop = ind + 1
    l = x/a
    y = -a * logmxp1(l)
    if a*eps()*eps() > 3.28e-3
        throw(DomainError((a,x,ind,"P(a,x) or Q(a,x) is computationally indeterminant in this case.")))
    end
    c = exp(-y)
    w = 0.5*erfcx(sqrt(y))
    u = 1.0/a
    if iop == 1
        c0 = @horner(z , d00 , d0[1] , d0[2] , d0[3])
        c1 = @horner(z , d10 , d1[1] , d1[2] , d1[3])
        c2 = @horner(z , d20 , d2[1] , d2[2])
        c3 = @horner(z , d30 , d3[1] , d3[2])
        c4 = @horner(z , d40 , d4[1])
        c5 = @horner(z , d50 , d5[1])
        c6 = @horner(z , d60 , d6[1])

        t = @horner(u , c0 , c1 , c2 , c3 , c4 , c5 , c6 , d70 , d80)
        
    elseif iop == 2
        c0 = @horner(d00 , d0[1] , d0[2])
        c1 = @horner(d10 , d1[1])
        t = @horner(u , c0 , c1 , d20)
    
    else
        t = @horner(z , d00 , d0[1])
        
    end
    if l < 1.0
        return (c*(w - rt2pin*t/sqrt(a)) , 1.0 - c*(w - rt2pin*t/sqrt(a)))
    end
    return (1.0 - c*(w + rt2pin*t/sqrt(a)) , c*(w + rt2pin*t/sqrt(a)))
end
"""
    gamma_p_fsum(a,x)

Compute using Finite Sums for Q(a,x) when a >= 1 && 2a is integer
Used when a <= x <= x0 and a = n/2.
Refer Eqn (14) in the paper.
"""
function gamma_p_fsum(a::Float64, x::Float64)
    if isinteger(a)           
        sm = exp(-x)
        t = sm
        N = 1
        c=0.0
        i = trunc(Int,a )
    else
        rtx = sqrt(x)
        sm = erfc(rtx)
        t = exp(-x)/(rtpi*rtx)
        N=0
        c=-0.5
        i = trunc(Int,a )
    end
    for lp = N : i-1
        if i == 0 
            break
        end
        c += 1.0
        t = (x*t)/c
        sm += t
    end
    return (1.0 - sm, sm)

end
# Reference : 'Computation of the incomplete gamma function ratios and their inverse' by Armido R DiDonato , Alfred H Morris.
# Published in Journal: ACM Transactions on Mathematical Software (TOMS)
# Volume 12 Issue 4, Dec. 1986 Pages 377-393
# doi>10.1145/22721.23109

"""
    gamma_inc(a,x,IND)
    
DLMF: https://dlmf.nist.gov/8.2#E4 , https://dlmf.nist.gov/8.2#E5
Wikipedia: https://en.wikipedia.org/wiki/Incomplete_gamma_function
IND --> Accuracy desired ; IND=0 means 14 significant digits accuracy , IND=1 means 6 significant digit and IND=2 means only 3 digit accuracy suffices.
gamma_p(a,x) or P(a,x) is the Incomplete gamma function ratio given by : ``1/\\Gamma (a) \\int_{0}^{x} e^{-t}t^{a-1} dt``
gamma_q(a,x) or Q(a,x) is the Incomplete gamma function ratio given by : 1 - P(a,x) ->  ``1/\\Gamma (a) \\int_{x}^{\\infty} e^{-t}t^{a-1} dt``
Returns a tuple (gamma_p, gamma_q) where gamma_p + gamma_q = 1.0
"""
function gamma_inc(a::Float64,x::Float64,ind::Integer)
    iop = ind + 1
    acc = acc0[iop]
    if a<0.0 || x<0.0
        throw(DomainError((a,x,ind,"`a` and `x` must be greater than 0 ---- Domain : (0,inf)")))
    elseif a==0.0 && x==0.0
        throw(DomainError((a,x,ind,"`a` and `x` must be greater than 0 ---- Domain : (0,inf)")))
    elseif a*x==0.0 
        if x<=a
            return (0.0,1.0)
        else
            return (1.0,0.0)
        end               
    end
    
    if a >= 1.0
        if a >= big1[iop]
            l = x/a
            if l == 0.0
                return (0.0,1.0)
            end
            s = 1.0 - l
            z = -logmxp1(l)
            if z >= 700.0/a
                return (0.0,1.0)
            end
            y = a*z
            rta = sqrt(a)
            if abs(s) <= e0[iop]/rta
                z = sqrt(z + z) 
                if l < 1.0 
                    z=-z
                return gamma_p_temme_1(a, x, z, ind)
                end
            end

            if abs(s) <= 0.4
                if abs(s) <= 2.0*eps() && a*eps()*eps() > 3.28e-3
                    throw(DomainError((a,x,ind,"P(a,x) or Q(a,x) is computationally indeterminant in this case.")))
                end
                c = exp(-y)
                w = 0.5*erfcx(sqrt(y))
                u = 1.0/a
                z = sqrt(z + z) 
                if l < 1.0 
                    z=-z
                end
                if iop == 1
                    return gamma_p_minimax(a,x,z)
                elseif iop == 2
                    return gamma_p_temme(a,x,z)
                else
                    t = @horner(z , d00 , d0[1] , d0[2] , d0[3])
                    return gamma_p_temme_1(a, x, z, ind)
                end
            end
        elseif a > x || x >= x0[iop] || !isinteger(2*a)  
            r = rgammax(a,x)
            if r == 0.0
                if x <= a
                    return (0.0,1.0)
                else
                    return (1.0,0.0)
                end
            end 
            if x <= max(a,alog10)
                return gamma_p_taylor(a, x, ind)
            elseif x < x0[iop]
                return gamma_p_cf(a, x, ind)
            else
                return gamma_p_asym(a, x, ind)
            end
        else
            return gamma_p_fsum(a,x)
            
        end
    elseif a == 0.5
        if x >= 0.25
            return ( 1.0 - erfc(sqrt(x)) , erfc(sqrt(x)) )
        end
        return ( erf(sqrt(x)) , 1.0 - erf(sqrt(x)) )
    elseif x < 1.1
        return gamma_p_taylor_x(a, x, ind)  
    end
    r = rgammax(a,x)
    if r == 0.0
        return (1.0, 0.0)
    else
        return gamma_p_cf(a, x, ind)    
    end

    
    
end

function gamma_inc(a::BigFloat,x::BigFloat,ind::Integer) #BigFloat version from GNU MPFR wrapped via ccall
    z = BigFloat()
    ccall((:mpfr_gamma_inc, :libmpfr), Int32 , (Ref{BigFloat} , Ref{BigFloat} , Ref{BigFloat} , Int32) , z , a , x , ROUNDING_MODE[])
    return (1.0 - z/gamma(a), z/gamma(a))
end
gamma_inc(a::Float32,x::Float32,ind::Integer) = ( Float32(gamma_inc(Float64(a),Float64(x),ind)[1]) , Float32(gamma_inc(Float64(a),Float64(x),ind)[2]) )
gamma_inc(a::Float16,x::Float16,ind::Integer) = ( Float16(gamma_inc(Float64(a),Float64(x),ind)[1]) , Float16(gamma_inc(Float64(a),Float64(x),ind)[2]) )
gamma_inc(a::Real,x::Real,ind::Integer) = (gamma_inc(float(a),float(x),ind))
gamma_inc(a::Integer,x::Integer,ind::Integer) = gamma_inc(Float64(a),Float64(x),ind)
gamma_inc(a::AbstractFloat,x::AbstractFloat,ind::Integer) = throw(MethodError(gamma_inc,(a,x,ind,"")))
