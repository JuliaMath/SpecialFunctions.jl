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

const eu = Base.MathConstants.eulergamma
const a0 = [3.31125922108741, 11.6616720288968, 4.28342155967104, .213623493715853]
const b0 = [6.61053765625462, 6.40691597760039, 1.27364489782223, .036117081018842]

const eps0 = [1.0E-10, 1.0E-08]
const amin = [500.0, 100.0]
const bmin = [1.0E-28, 1.0E-13]
const dmin = [1.0E-6, 1.0E-4]
const emin = [2.0E-03, 6.0E-03]
const tol = 1.0E-5
const xmin = floatmin(Float64)
const xmax = floatmax(Float64)

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
# #EVALUATE LOG(1+A)
# function loga1(a::Float64)
#     if abs(a) > .375
#         xx = 1.0 + a
#         if a<0.0
#             xx = a + 1.0
#         end
#         return log(xx)
#     else
#         t = a/(a + 2.0)
#         t2 = t*t
#         w = @horner(t2, 1.0, -.129418923021993E+01, .405303492862024E+00, -.178874546012214E-01) / @horner(t2, 1.0, -.162752256355323E+01, .747811014037616E+00, -.845104217945565E-01)
#         return 2.0*t*w
#     end
# end
#EVALUATE LOG(GAMMA(1+A)) FOR -0.2 <= A <= 1.25
function log1gamma1(a::Float64)
  
    if a >= 0.6
        x = a - 1.0
        w = @horner(x, .422784335098467E+00, .848044614534529E+00, .565221050691933E+00, .156513060486551E+00, .170502484022650E-01, .497958207639485E-03) / @horner(x, 1.0, .124313399877507E+01, .548042109832463E+00, .101552187439830E+00, .713309612391000E-02, .116165475989616E-03)
        return x*w
    else
        w = @horner(a, .577215664901533E+00, .844203922187225E+00, -.168860593646662E+00, -.780427615533591E+00, -.402055799310489E+00, -.673562214325671E-01, -.271935708322958E-02) / @horner(a, 1.0, .288743195473681E+01, .312755088914843E+01, .156875193295039E+01, .361951990101499E+00, .325038868253937E-01, .667465618796164E-03)
        return -a*w
    end
end
# #EVALUATE INVERSE NORMAL DISTRIBUTION
# function normal_inv(p::Float64, q::Float64, d::Float64)
#     t = min(p,q)
#     w = 1.0 - (p+q)
#     if abs(w) > 2.0*eps()
#         throw(DomainError())
#     end
#     u = abs(2*d)
#     v = 2*t
#     w = erfinv(u)
#     if w < 0.0
#         throw(DomainError)
#     end
#     w *= sqrt(2)
#     if d < 0.0 
#         w = -w
#     end
#     return w
# end
"""
    chepolsum(n,x,a)

Computes a series of Chebyshev Polynomials given by : a[0]/2 + a[1]T1(x) + .... + a[n]Tn(X)
"""
function chepolsum(n::Integer, x::Float64, a::Array{Float64,1})
  
    if n == 0
        return a[1]/2.0
    elseif n == 1
        return a[1]/2.0 + a[1]*x
    else
        tx = 2*x
        r = a[n]
        h = a[n-1] + r*tx
        for k = n-2:-1:1
            s=r
            r=h
            h=a[k]+r*tx-s
        end
        return a[0]/2.0 - r + h*x
    end
end
"""
    stirling(x)

Stirling Series
"""
function stirling(x::Float64)
  
    a = zeros(18) 
    if x < floatmin(Float64)*1000.0
        return floatmax(Float64)/1000.0
    elseif x < 1
        return log1gamma1(x) - (x+0.5)*log(x)+x - log(sqrt(2*pi))
    elseif x < 2
        return log1gamma1(x-1) - (x-0.5)*log(x)+x - log(sqrt(2*pi))
    elseif x < 3
        return log1gamma1(x-2) - (x-0.5*log(x)+x  - log(sqrt(2*pi))+log(x-1))
    elseif x < 12
        a[1]=1.996379051590076518221
        a[2]=-0.17971032528832887213e-2
        a[3]=0.131292857963846713e-4
        a[4]=-0.2340875228178749e-6
        a[5]=0.72291210671127e-8
        a[6]=-0.3280997607821e-9
        a[7]=0.198750709010e-10
        a[8]=-0.15092141830e-11
        a[9]=0.1375340084e-12
        a[10]=-0.145728923e-13
        a[11]=0.17532367e-14
        a[12]=-0.2351465e-15
        a[13]=0.346551e-16
        a[14]=-0.55471e-17
        a[15]=0.9548e-18
        a[16]=-0.1748e-18
        a[17]=0.332e-19
        a[18]=-0.58e-20
        z=18.0/(x*x)-1.0
        return chepolsum(17,z,a)/(12.0*x)
    else
        c = zeros(7)
        z = 1.0/(x*x)
        if x < 1000
            c[1]=0.25721014990011306473e-1
            c[2]=0.82475966166999631057e-1
            c[3]=-0.25328157302663562668e-2
            c[4]=0.60992926669463371e-3
            c[5]=-0.33543297638406e-3
            c[6]=0.250505279903e-3
            c[7]=0.30865217988013567769
            return @horner(z,c[1],c[2],c[3],c[4],c[5],c[6])/((c[7]+z)/x)
        else
            return (((-z/1680.0+1.0/1260.0)*z-1.0/360.0)*z+1.0/12.0)/x
        end
    end
end
function gamstar(x::Float64)
    if x>=3
        return exp(stirling(x))
    elseif x>0
        return gamma(x)/(exp(-x+(x-0.5)*log(x))*sqrt(2*pi))
    else
        return floatmax(Float64)/1000.0
    end
end
function lambdaeta(eta::Float64)
  
    ak = zeros(7)
    s = eta*eta*0.5
    if eta == 0.0
        la = 1
    elseif eta < -1.0
        r = exp(-1-s)
        ak[1]=1.0
        ak[2]=1.0
        ak[3]=1.5
        ak[4]=8.0/3.0
        ak[5]=125.0/24.0
        ak[6]=54.0/5.0
        la = @horner(r, 0.0, ak[1], ak[2], ak[3], ak[4], ak[5], ak[6])
    elseif eta < 1.0
        ak[1] = 1.0
        ak[2] = 1.0/3.0
        ak[3] = 1.0/36.0
        ak[4] = -1.0/270.0
        ak[5] = 1.0/4320.0
        ak[6] = 1.0/17010.0
        r = eta
        la = @horner(r, 1.0, ak[1], ak[2], ak[3], ak[4], ak[5], ak[6])
    else
        r = 11 + s
        l = log(r)
        la = r+l
        r = 1.0/r
        l2=l*l
        l3=l2*l
        l4=l3*l
        l5=l4*l
        ak[1] = 1.0
        ak[2] = (2 - l)*0.5
        ak[3] = (-9*l+6+2*l2)/6.0
        ak[4] = -(3*l3+36*l-22*l2-12)/12.0
        ak[5] = (60 + 350*l2 - 300*l -125*l3 + 12*l4)/60.0
        ak[6] = -(-120-274*l4+900*l-1700*l2+1125*l3+20*l5)/120.0
        la = la + l*@horner(r,0.0,ak[1],ak[2],ak[3],ak[4],ak[5],ak[6])
    end
    r = 1
    if (eta>-3.5 && eta<-0.03) || (eta>0.03 && eta<40)
        r=1
        q=la
        while r>1.0e-8
            la = q*(s+log(q))/(q-1.0)
            r = abs(q/la-1)
            q=la
        end
    end
    return la
end
function ratfun(x::Float64, ak::Array{Float64,1}, bk::Array{Float64,1})
  
    p = @horner(x, ak[1], ak[2], ak[3], ak[4], ak[5])
    q = @horner(x, bk[1], bk[2], bk[3], bk[4], bk[5])
    return p/q
end

function eps1(eta::Float64)
    ak=zeros(5)
    bk=zeros(5)
    if abs(eta) < 1.0
        ak[1] = -3.333333333438e-1
        ak[2] = -2.070740359969e-1
        ak[3] = -5.041806657154e-2
        ak[4] = -4.923635739372e-3
        ak[5] = -4.293658292782e-5
        bk[1] = 1.000000000000e+0
        bk[2] = 7.045554412463e-1
        bk[3] = 2.118190062224e-1
        bk[4] = 3.048648397436e-2
        bk[5] = 1.605037988091e-3
        eps1 = ratfun(eta,ak,bk)
    else
        la = lambdaeta(eta)
        eps1 = log(eta/(la-1.0))/eta
    end
    return eps1
end

function eps2(eta::Float64)
    ak=zeros(5)
    bk=zeros(5)

    if eta < -5.0
        x=eta*eta
        lnmeta = log(-eta)
        eps2 = (12.0 - x - 6.0*lnmeta*lnmeta)/(12.0*x*eta)
    elseif eta < -2.0
        ak[1] = -1.72847633523e-2
        ak[2] = -1.59372646475e-2
        ak[3] = -4.64910887221e-3
        ak[4] = -6.06834887760e-4
        ak[5] = -6.14830384279e-6
        bk[1] = 1.00000000000e+0
        bk[2] = 7.64050615669e-1
        bk[3] = 2.97143406325e-1
        bk[4] = 5.79490176079e-2
        bk[5] = 5.74558524851e-3
        eps2 = ratfun(eta,ak,bk)
    elseif eta < 2.0
        ak[1] = -1.72839517431e-2
        ak[2] = -1.46362417966e-2
        ak[3] = -3.57406772616e-3
        ak[4] = -3.91032032692e-4
        ak[5] = 2.49634036069e-6
        bk[1] = 1.00000000000e+0
        bk[2] = 6.90560400696e-1
        bk[3] = 2.49962384741e-1
        bk[4] = 4.43843438769e-2
        bk[5] = 4.24073217211e-3
        eps2 = ratfun(eta,ak,bk)
    elseif eta < 1000.0
        ak[1] = 9.99944669480e-1
        ak[2] = 1.04649839762e+2
        ak[3] = 8.57204033806e+2
        ak[4] = 7.31901559577e+2
        ak[5] = 4.55174411671e+1
        bk[1] = 1.00000000000e+0
        bk[2] = 1.04526456943e+2
        bk[3] = 8.23313447808e+2
        bk[4] = 3.11993802124e+3
        bk[5] = 3.97003311219e+3
        x = 1.0/eta
        eps2 = ratfun(x,ak,bk)/(-12.0*eta)
    else
        eps2 = -1.0/(12.0*eta)
    end
    return eps2
end
function eps3(eta::Float64)
    ak=zeros(5)
    bk=zeros(5)
    if eta < -8.0
        x=eta*eta
        y = log(-eta)/eta
        eps3=(-30.0+eta*y*(6.0*x*y*y-12.0+x))/(12.0*eta*x*x)
    elseif eta < -4.0
        ak[1] = 4.95346498136e-2
        ak[2] = 2.99521337141e-2
        ak[3] = 6.88296911516e-3
        ak[4] = 5.12634846317e-4
        ak[5] = -2.01411722031e-5
        bk[1] = 1.00000000000e+0
        bk[2] = 7.59803615283e-1
        bk[3] = 2.61547111595e-1
        bk[4] = 4.64854522477e-2
        bk[5] = 4.03751193496e-3
        eps3 = ratfun(eta,ak,bk)/(eta*eta)
    elseif eta < -2.0
        ak[1] = 4.52313583942e-3
        ak[2] = 1.20744920113e-3
        ak[3] = -7.89724156582e-5
        ak[4] = -5.04476066942e-5
        ak[5] = -5.35770949796e-6
        bk[1] = 1.00000000000e+0
        bk[2] = 9.12203410349e-1
        bk[3] = 4.05368773071e-1
        bk[4] = 9.01638932349e-2
        bk[5] = 9.48935714996e-3
        eps3 = ratfun(eta,ak,bk)
    elseif eta < 2.0
        ak[1] = 4.39937562904e-3
        ak[2] = 4.87225670639e-4
        ak[3] = -1.28470657374e-4
        ak[4] = 5.29110969589e-6
        ak[5] = 1.57166771750e-7
        bk[1] = 1.00000000000e+0
        bk[2] = 7.94435257415e-1
        bk[3] = 3.33094721709e-1
        bk[4] = 7.03527806143e-2
        bk[5] = 8.06110846078e-3
        eps3 = ratfun(eta,ak,bk)
    elseif eta < 10.0
        ak[1] = -1.14811912320e-3
        ak[2] = -1.12850923276e-1
        ak[3] = 1.51623048511e+0
        ak[4] = -2.18472031183e-1
        ak[5] = 7.30002451555e-2
        bk[1] = 1.00000000000e+0
        bk[2] = 1.42482206905e+1
        bk[3] = 6.97360396285e+1
        bk[4] = 2.18938950816e+2
        bk[5] = 2.77067027185e+2
        x=1.0/eta
        eps3 = ratfun(x,ak,bk)/(eta*eta)
    elseif eta < 100.0
        ak[1] = -1.45727889667e-4
        ak[2] = -2.90806748131e-1
        ak[3] = -1.33085045450e+1
        ak[4] = 1.99722374056e+2
        ak[5] = -1.14311378756e+1
        bk[1] = 1.00000000000e+0
        bk[2] = 1.39612587808e+2
        bk[3] = 2.18901116348e+3
        bk[4] = 7.11524019009e+3
        bk[5] = 4.55746081453e+4
        x=1.0/eta
        eps3 = ratfun(x,ak,bk)/(eta*eta)
    else
        eta3 = eta*eta*eta
        eps3 = -log(eta)/(12.0*eta3)
    end
    return eps3
end


"""
    gamma_inc_cf(a, x, ind)

Computes P(a,x) by continued fraction expansion given by : 
```math
R(a,x) * \\frac{1}{1-\\frac{z}{a+1+\\frac{z}{a+2-\\frac{(a+1)z}{a+3+\\frac{2z}{a+4-\\frac{(a+2)z}{a+5+\\frac{3z}{a+6-\\dots}}}}}}}
```
Used when 1 <= a <= BIG and x < x0.
DLMF : https://dlmf.nist.gov/8.9#E2
"""
function gamma_inc_cf(a::Float64, x::Float64, ind::Integer)
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
    q = rgammax(a,x)*a2n
    return (1.0 - q, q)
end
"""
    gamma_inc_taylor(a, x, ind)

Compute P(a,x) using Taylor Series for P/R given by : 
```math
R(a,x)/a * (1 + \\sum_{1}^{\\infty} x^{n}/((a+1)(a+2)...(a+n)))
```
Used when 1 <= a <= BIG and x <= max{a, ln 10}.
DLMF : https://dlmf.nist.gov/8.11#E2
"""
function gamma_inc_taylor(a::Float64, x::Float64, ind::Integer)
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
    p = (rgammax(a,x)/a)*(1.0 + sm)
    return (p, 1.0-p)
end
"""
    gamma_inc_asym(a, x, ind)

Compute P(a,x) using asymptotic expansion given by : 
```math
R(a,x)/a * (1 + \\sum_{1}^{N-1}(a_{n}/x^{n} + \\Theta _{n}a_{n}/x^{n}))
```
where R(a,x) = rgammax(a,x). Used when 1 <= a <= BIG and x >= x0.
DLMF : https://dlmf.nist.gov/8.11#E2
"""
function gamma_inc_asym(a::Float64, x::Float64, ind::Integer)
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
    q = (rgammax(a,x)/x)*(1.0 + sm)
    return (1.0 - q, q) 
end
"""
    gamma_inc_taylor_x(a,x,ind)

Computes P(a,x) based on Taylor expansion of P(a,x)/x**a given by:
```math
J = -a * \\sum_{1}^{\\infty} (-x)^{n}/((a+n)n!)
``` and P(a,x)/x**a is given by :
```math
(1 - J)/ \\Gamma(a+1)
``` resulting from term-by-term integration of gamma_inc(a,x,ind).
This is used when a < 1 and x < 1.1 - Refer Eqn (9) in the paper.
"""
function gamma_inc_taylor_x(a::Float64, x::Float64, ind::Integer)
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
       q = max((w*temp - l)*g - h, 0.0)
       return (1.0 - q, q)
    else
       w = exp(z)
       p = w*g*(1.0 - temp)
       return (p, 1.0 - p)
    end
end
"""
    gamma_inc_minimax(a,x,z)

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
function gamma_inc_minimax(a::Float64, x::Float64, z::Float64)
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
            p = c*(w - rt2pin*t/sqrt(a))
            return (p , 1.0 - p)
        else
            q = c*(w + rt2pin*t/sqrt(a))
            return (1.0 - q , q)
        end
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
        p = c*(w - rt2pin*t/sqrt(a))
        return (p , 1.0 - p)
    else
        q = c*(w + rt2pin*t/sqrt(a))
        return (1.0 - q , q)
    end
end
"""
    gamma_inc_temme(a, x, z)

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
function gamma_inc_temme(a::Float64, x::Float64, z::Float64)
    l = x/a
    y = -a*logmxp1(x/a)
    c = exp(-y)
    w = 0.5*erfcx(sqrt(y))
    c0 = @horner(z , d00 , d0[1] , d0[2] , d0[3] , d0[4] , d0[5] , d0[6]) 
    c1 = @horner(z , d10 , d1[1] , d1[2] , d1[3] , d1[4]) 
    c2 = @horner(z , d20 , d2[1])
    t = @horner(1.0/a , c0 , c1 , c2)
    if l < 1.0
        p = c*(w - rt2pin*t/sqrt(a))
        return (p , 1.0 - p)
    else
        q = c*(w + rt2pin*t/sqrt(a))
        return (1.0 - q , q)
    end
end
"""
    gamma_inc_temme_1(a, x, z, ind)

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
function gamma_inc_temme_1(a::Float64, x::Float64, z::Float64, ind::Integer)
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
        p = c*(w - rt2pin*t/sqrt(a))
        return (p , 1.0 - p)
    else
        q = c*(w + rt2pin*t/sqrt(a))
        return (1.0 - q , q)
    end
end
"""
    gamma_inc_fsum(a,x)

Compute using Finite Sums for Q(a,x) when a >= 1 && 2a is integer
Used when a <= x <= x0 and a = n/2.
Refer Eqn (14) in the paper.
"""
function gamma_inc_fsum(a::Float64, x::Float64)
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
    q = sm
    return (1.0 - q, q)

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
gamma_inc(a,x) or P(a,x) is the Incomplete gamma function ratio given by : ``1/\\Gamma (a) \\int_{0}^{x} e^{-t}t^{a-1} dt``
gamma_q(a,x) or Q(a,x) is the Incomplete gamma function ratio given by : 1 - P(a,x) ->  ``1/\\Gamma (a) \\int_{x}^{\\infty} e^{-t}t^{a-1} dt``
Returns a tuple (gamma_inc, gamma_q) where gamma_inc + gamma_q = 1.0
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
                end
                return gamma_inc_temme_1(a, x, z, ind)
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
                    return gamma_inc_minimax(a,x,z)
                elseif iop == 2
                    return gamma_inc_temme(a,x,z)
                else
                    t = @horner(z , d00 , d0[1] , d0[2] , d0[3])
                    return gamma_inc_temme_1(a, x, z, ind)
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
                return gamma_inc_taylor(a, x, ind)
            elseif x < x0[iop]
                return gamma_inc_cf(a, x, ind)
            else
                return gamma_inc_asym(a, x, ind)
            end
        else
            return gamma_inc_fsum(a,x)
            
        end
    elseif a == 0.5
        if x >= 0.25
            q = erfc(sqrt(x))
            return ( 1.0 - q , q )
        end
        p = erf(sqrt(x))
        return ( p , 1.0 - p )
    elseif x < 1.1
        return gamma_inc_taylor_x(a, x, ind)  
    end
    r = rgammax(a,x)
    if r == 0.0
        return (1.0, 0.0)
    else
        return gamma_inc_cf(a, x, ind)    
    end

    
    
end

function gamma_inc(a::BigFloat,x::BigFloat,ind::Integer) #BigFloat version from GNU MPFR wrapped via ccall
    z = BigFloat()
    ccall((:mpfr_gamma_inc, :libmpfr), Int32 , (Ref{BigFloat} , Ref{BigFloat} , Ref{BigFloat} , Int32) , z , a , x , ROUNDING_MODE[])
    q = z/gamma(a)
    return (1.0 - q, q)
end
gamma_inc(a::Float32,x::Float32,ind::Integer) = ( Float32(gamma_inc(Float64(a),Float64(x),ind)[1]) , Float32(gamma_inc(Float64(a),Float64(x),ind)[2]) )
gamma_inc(a::Float16,x::Float16,ind::Integer) = ( Float16(gamma_inc(Float64(a),Float64(x),ind)[1]) , Float16(gamma_inc(Float64(a),Float64(x),ind)[2]) )
gamma_inc(a::Real,x::Real,ind::Integer) = (gamma_inc(float(a),float(x),ind))
gamma_inc(a::Integer,x::Integer,ind::Integer) = gamma_inc(Float64(a),Float64(x),ind)
gamma_inc(a::AbstractFloat,x::AbstractFloat,ind::Integer) = throw(MethodError(gamma_inc,(a,x,ind,"")))

#EFFICIENT AND ACCURATE ALGORITHMS FOR THECOMPUTATION AND INVERSION OF THE INCOMPLETEGAMMA FUNCTION RATIOS by Amparo Gil, Javier Segura, Nico M. Temme
#SIAM Journal on Scientific Computing 34(6) (2012), A2965-A2981
# arXiv:1306.1754

"""
    gamma_inc_inv(a,x,IND)
    
DLMF: https://dlmf.nist.gov/8.2#E4 , https://dlmf.nist.gov/8.2#E5
Wiki: https://en.wikipedia.org/wiki/Incomplete_gamma_function

gamma_inc(a,x) or P(a,x) is the Incomplete gamma function ratio given by : ``1/\\Gamma (a) \\int_{0}^{x} e^{-t}t^{a-1} dt``
"""
function gamma_inc_inv(a::Float64, p::Float64, q::Float64)
    ck = zeros(20)
    if p < 0.5
        pcase = true
        porq = p
        s=-1
    else
        pcase = false
        porq = q
        s=1
    end

    logr = (1.0/a)*(log(p) + logabsgamma(a + 1.0)[1])
    if logr < log(0.2*(1+a))
        r = exp(logr)
        m=0
        a2=a*a
        a3=a2*a
        a4=a3*a
        ap1=a + 1.0
        ap12=(a+1.0)*ap1
        ap13=(a+1.0)*ap12
        ap14=ap12*ap12
        ap2 = a+2.0
        ap22 = ap2*ap2
        ck[1]= 1.0
        ck[2]= 1.0/(1.0+a)
        ck[3]=0.5*(3*a+5)/(ap12*(a+2))
        ck[4]= (1.0/3.0)*(31+8*a2+33*a)/(ap13*ap2*(a+3))
        ck[5]= (1.0/24.0)*(2888+1179*a3+125*a4+3971*a2+5661*a)/(ap14*ap22*(a+3)*(a+4))
        x0 = @horner(r, 0.0, ck[1], ck[2], ck[3], ck[4], ck[5])
    elseif ((q < min(0.02,exp(-1.5*a)/gamma(a))) && (a<10))
        m =0
        b=1.0-a
        b2=b*b
        b3=b2*b
        eta=sqrt(-2/a*log(q*gamstar(a)*sqrt(2*pi)/sqrt(a)))
        x0 = a*lambdaeta(eta)
        l = log(x0)

        if a > 0.12 || x0 > 5
            l2 = l*l
            l3 = l2*l
            l4 = l3*l
            r = 1.0/x0
            ck[1] = l - 1.0
            ck[2] = (3*b-2*b*l+l2-2*l+2)/2.0
            ck[3] = (24*b*l-11*b2-24*b-6*l2+12*l-12-9*b*l2+6*b2*l+2*l3)/6.0
            ck[4] = (-12*b3*l+84*b*l2-114*b2*l+72+36*l2+3*l4-72*l+162*b-168*b*l-12*l3+25*b3-22*b*l3+36*b2*l2+120*b2)/12.0
            x0 = x0 - l + b * r * @horner(r,ck[1],ck[2],ck[3],ck[4])
        else
            r = 1.0/x0
            l2 = l*l
            ck[1] = l - 1.0
            x0 = x0 - l + b * r * ck[1]
        end
    elseif abs(porq - 0.5) < 1.0e-05
        m=0
        x0=a-1.0/3.0+(8.0/405.0+184.0/25515.0/a)/a
    elseif abs(a-1.0) < 1.0e-4
        m=0
        if pcase
            x0 = -log(1.0-p)
        else
            x0 = -log(q)
        end
    elseif a < 1.0
        m=0
        if pcase
            x0=exp((1.0/a)*(log(porq)+logabsgamma(a+1.0)[1]))
        else
            x0=exp((1.0/a)*(log(1.0-porq)+logabsgamma(a+1.0)[1]))
        end
    else
        m=1
        r = erfcinv(2*porq)
        eta = s*r/sqrt(a*0.5)
        eta += (eps1(eta)+(eps2(eta)+eps3(eta)/a)/a)/a
        x0 = a*lambdaeta(eta)
    end

    t=1
    x=x0
    n=1
    a2=a*a
    a3=a2*a
    #println(x)
   
    #Newton iteration
    while t > 1.0e-15 && n < 15
        x=x0
        x2=x*x
        if m==0
            dlnr = (1.0-a)*log(x)+x+logabsgamma(a)[1]
            if dlnr > log(xmax/1000.0)
                n=20
            else
                r = exp(dlnr)
                if pcase
                    (px,qx) = gamma_inc(a,x,0)
                    ck[1] = -r*(px-p)
                else
                    (px,qx) = gamma_inc(a,x,0)
                    ck[1] = r*(qx-q)
                end
                ck[2] = (x-a+1.0)/(2.0*x)
                ck[3] = (2*x2-4*x*a+4*x+2*a2-3*a+1)/(6*x2)
                r = ck[1]
                if a > 0.1
                    x0 = x + @horner(r,0.0,1.0,ck[2],ck[3])
                else
                    if a > 0.05
                        x0 = x + @horner(r,0.0,1.0,ck[2])
                    else
                        x0=x+r
                    end
                end
            end
        else
            y=eta
            fp = -sqrt(a/(2*pi))*exp(-0.5*a*y*y)/gamstar(a)
            r = -(1/fp)*x
            if pcase
                (px,qx) = gamma_inc(a,x,0)
                ck[1] = -r*(px-p)
            else
                (px,qx) = gamma_inc(a,x,0)
                ck[1] = r*(qx-q)
            end
            ck[2] = (x-a+1.0)/(2.0*x)
            ck[3] = (2*x2-4*x*a+4*x+2*a2-3*a+1.0)/(6.0*x2)
            r = ck[1]

            if a > 0.1
                x0 = x + @horner(r,0.0,1.0,ck[2],ck[3])
            else
                if a > 0.05
                    x0 = x + @horner(r,0.0,1.0,ck[2])
                else
                    x0=x+r
                end
            end
        end
        t=abs(x/x0 - 1.0)
        n+=1
        x=x0
       
    end
    xr=x
    return xr
end



# function gamma_inc_inv(a::Float64, p::Float64, q::Float64)
#     x = 0.0
#     if a <= 0.0
#         @goto l500
#     end
#     t = p + q - 1.0
#     if abs(t) > eps()
#         @goto l520
#     end
#     if p == 0.0
#         return 0.0
#     end
#     if q == 0.0
#         @goto l400
#     end
#     if a == 1.0
#         @goto l410
#     end
    
#     e2 = 2*eps()
#     EPS = max(100.0*eps(), 1.0E-10)
#     amax = (0.4E-10)/(eps()*eps())
#     iop = 1
#     if eps() > 1.0E-10
#         iop = 2
#     end
#     epps = eps0[iop]
#     x=0.0
#     xn = x
# #-----SELECTION OF INIT. APPROX XN OF X WHEN A < 1---------

#     if a > 1.0
#         @goto l30
#     end
#     g = gamma(a + 1.0)
#     qg = q*g
#     if qg == 0.0
#         @goto l560
#     end
#     b = qg/a
#     if qg > 0.6*a
#         @goto l20
#     end
#     if a >= 0.3 || b < 0.35
#         @goto l10
#     else
#         t = exp(-(b+eu))
#         u = t*exp(t)
#         xn = t*exp(u)
#         @goto l100
#     end

#     @label l10
#      if b >= 0.45
#         @goto l20
#      end
#      if b == 0.0
#         @goto l560
#      end
#      y = -log(b)
#      s = 1.0 - a
#      z = log(y)
#      t = y - s*z
#      if b < 0.15
#         @goto l11
#      end
#      xn = y - s*log(t) - log(1.0 + s/(t + 1.0))
#      @goto l200

#     @label l11
#      if b <= 0.01
#         @goto l12
#      end
#      u = ((t + 2.0*(3.0 - a))*t + (2.0 - a)*(3.0 - a))/((t + (5.0 - a))*t + 2.0)
#      xn = y - s*log(t) - log(u)
#      @goto l200
     
#     @label l12
#      c1 = -s*z
#      c2 = -s*(1.0 + c1)
#      c3 =  s*((0.5*c1 + (2.0 - a))*c1 + (2.5 - 1.5*a))
#      c4 = -s*(((c1/3.0 + (2.5 - 1.5*a))*c1 + ((a - 6.0)*a + 7.0))*c1+ ((11.0*a - 46.0)*a + 47.0)/6.0)
#      c5 = -s*((((-c1/4.0 + (11.0*a - 17.0)/6.0)*c1 + ((-3.0*a + 13.0)*a - 13.0))*c1 + 0.5*(((2.0*a - 25.0)*a + 72.0)*a - 61.0))*c1 + (((25.0*a - 195.0)*a + 477.0)*a - 379.0)/12.0)
#      xn = ((((c5/y + c4)/y + c3)/y + c2)/y + c1) + y
#      if a > 1.0
#         @goto l200
#      end
#      if b > bmin[iop]
#         @goto l200
#      end
#      x = xn
#      return xn

#     @label l20
#      if b*q > 1.0E-8
#         @goto l21
#      end
#      xn = exp(-(q/a + eu))
#      @goto l23
#     @label l21
#      if p <= 0.9
#         @goto l22
#      end
#      xn = exp((loga1(-q) + log1gamma1(a))/a)
#      @goto l23
#     @label l22
#      xn = exp(log(p*g)/a)
#     @label l23
#      if xn == 0.0
#         @goto l510
#         t = 1.0 - xn/(a + 1.0)
#         xn /= t
#         @goto l100
#      end

# #SELECTION OF INIT APPROX XN OF X WHEN A>1
    
#     # @label l30
#     #  if q <= 0.5
#     #     @goto l31
#     #  else
#     #     w = log(p)
#     #     @goto l32
#     #  end
#     # @label l31
#     #  w = log(q)
#     # @label l32
#     #  t = sqrt(-2.0*w)
#     #  s = t - @horner(t, a0[1], a0[2], a0[3], a0[4]) / @horner(t, 1.0, b0[1], b0[2], b0[3], b0[4])
#     #  if q > 0.5
#     #     s = -s
#     #  end
#     @label l30
#      t = p - .5
#      if q < .5
#         t = .5 - q
#      end
#      s = normal_inv(p,q,t)
#      rta = sqrt(a)
#      s2 = s*s
#      xn = (((12.0*s2 - 243.0)*s2 - 923.0)*s2 + 1472.0)/204120.0
#      xn = (xn/a + s*((9.0*s2 + 256.0)*s2 - 433.0)/(38880.0*rta)) - ((3.0*s2 + 7.0)*s2 - 16.0)/810.0
#      xn = a + s*rta + (s2 - 1.0)/3.0 + s*(s2 - 7.0)/(36.0*rta) + xn/a
#      xn = max(xn, 0.0)
#      amin = 20.0
#      if eps() < 1.0E-8
#         amin = 250.0
#      end
#      if a < amin
#         @goto l40
#      end
#      x = xn
#      d = 1.0 - x/a
#      if abs(d) <= 0.1 #ch
#         return x
#      end
    
#     @label l40
#      if p <= 0.5
#         @goto l50
#      end
#      if xn < 3.0*a
#         @goto l200
#      end
#      w = log(q)
#      y = -(w + logabsgamma(a)[1])
#      d = max(2.0, a*(a - 1.0))
#      if y < d*log(10)
#         @goto l41
#      end
#      s = 1.0 - a
#      z = log(y)
#      @goto l12

#     @label l41
#      t = a - 1.0
#      println(xn)
#      xn = y + t*log(xn) - loga1(-t/(xn + 1.0))
#      println(xn)
#      xn = y + t*log(xn) - loga1(-t/(xn + 1.0))
#      @goto l200
    
#     @label l50
#      ap1 = a+1.0
#      if xn > 0.7*ap1
#         @goto l101
#      end
#      w = log(p) + logabsgamma(ap1)[1]
#      if xn > 0.15*ap1
#         @goto l60
#      end
#      ap2 = a + 2.0
#      ap3 = a + 3.0
#      x = exp((w+x)/a)
#      x = exp((w + x - log(1.0 + (x/ap1)*(1.0 + x/ap2)))/a)
#      x = exp((w + x - log(1.0 + (x/ap1)*(1.0 + x/ap2)))/a)
#      x = exp((w + x - log(1.0 + (x/ap1)*(1.0 + (x/ap2)*(1.0 + x/ap3))))/a)
#      xn = x
#      if xn > 1.0E-2 * ap1
#         @goto l60
#      end
#      if xn <= emin[iop]*ap1
#         return xn
#      end
#      @goto l101

#     @label l60
#      apn  = ap1
#      t = xn/apn
#      sm = 1.0 + t

#     @label l61
#      apn += 1.0
#      t = t*(xn/apn)
#      sm += t
#      if t > 1.0E-4
#         @goto l61
#      end
#      t = w - log(sm)
#      xn = exp((xn + t)/a)
#      xn = xn*(1.0 - (a*log(xn) - xn - t)/(a - xn))
#      @goto l101
    
#     #------SCHRODER ITERATION USING P--------

#     @label l100
#      if p > 0.5
#         @goto l200
#      end
#     @label l101 
#      if p <= xmin
#         @goto l550
#      end
#      am1 = a - 1.0
#     @label l102
#      if a <= amax
#         @goto l110
#      end
#      d = 1.0 - xn/a
#      if abs(d) <= 2.0*eps()
#         @goto l550
#      end

#     @label l110
#      (pn,qn) = gamma_inc(a,xn,0)
#      println(pn)
#      if pn == 0.0 || qn == 0.0
#         println("lol")
#         @goto l550
#      end
#      r = rgammax(a, xn)
#      if r == 0.0
#         println("hh")
#         @goto l550
#      end
#      t = (pn - p)/r
#      w = 0.5*(am1 - xn)
#      if abs(t) <= 0.1 && abs(w*t) <= 0.1
#         @goto l120
#      end
#      x = xn*(1.0 - t)
#      if x < 0.0
#         @goto l540
#      end
#      d = abs(t)
#      @goto l121

#     @label l120
#      h = t*(1.0 + w*t)
#      x = xn*(1.0 - h)
#      if x < 0.0
#         @goto l540
#      end
#      if abs(w) >= 1.0 && abs(w)*t*t <= epps
#         return x
#      end
#      d = abs(h)

#     @label l121
#      xn = x
#      if d > tol
#         @goto l102
#      end
#      if d <= epps
#         return xn
#      end
#      if abs(p - pn) <= tol*p
#         return xn
#      end
#      @goto l102

# #-----SCHRODER ITERATION USING Q----------
#     @label l200
#      if q <= 1.0E10 * xmin
#         @goto l550
#      end
#      am1 = a - 1.0
#     @label l201 
#      if a <= amax
#         @goto l210
#      end
#      d = 1.0 - xn/a
#      if abs(d) <= 2*eps()
#         @goto l550
#      end

#     @label l210
#      (pn,qn) = gamma_inc(a, xn, 0)
#      if pn == 0.0 || qn == 0.0
#         @goto l550
#      end
#      r = rgammax(a, xn)
#      if r == 0.0
#         @goto l550
#      end
#      t = (q - qn)/r
#      w = 0.5*(am1 - xn)
#      if abs(t) <= 0.1 && abs(w*t) <= 0.1
#         @goto l220
#      end
#      x = xn*(1.0 - t)
#      if x <= 0.0
#         @goto l540
#      end
#      d = abs(t)
#      @goto l221

#     @label l220
#      h = t*(1.0 + w*t)
#      x = xn*(1.0 - h)
#      if x < 0.0
#         @goto l540
#      end
#      if abs(w) >= 1.0 && abs(w)*t*t <= epps
#         return x
#      end
#      d = abs(h)

#     @label l221
#      xn = x
#      if d > tol
#         @goto l201
#      end
#      if d <= epps
#         return xn
#      end
#      if abs(q-qn) <= tol*q
#         return xn
#      end
#      @goto l201

#     @label l400
#      return xmax
#     @label l410
#      if q < 0.9
#         @goto l411
#      end
#      x = -loga1(-p)
#      return x
    
#     @label l411
#      x = -log(q)
#      return x

#     @label l500
#      print("-2 -- error")
#      return 0.0
#     @label l510
#      print("-3 -- error")
#      return 0.0
#     @label l520
#      print("-4 -- error")
#      return 0.0
#     @label l530
#      print("-6 -- error")
#      return 0.0
#     @label l540
#      print("-7 - error")
#      return 0.0
#     @label l550
#      print("-8 - error")
#      return xn
#     @label l560
#      x = xmax
#      print("-8 - error")
#      return x

# end

gamma_inc_inv(a::Float32, p::Float32, q::Float32) = Float32( gamma_inc_inv(Float64(a), Float64(p), Float64(q)) )
gamma_inc_inv(a::Float16, p::Float16, q::Float16) = Float16( gamma_inc_inv(Float64(a), Float64(p), Float64(q)) )
gamma_inc_inv(a::Real, p::Real, q::Real) = ( gamma_inc_inv(float(a), float(p), float(q)) )
gamma_inc_inv(a::Integer, p::Integer, q::Integer) = ( gamma_inc_inv(Float64(a), Float64(p), Float64(q)) )
gamma_inc_inv(a::Float32, p::Float32, q::Float32) = throw(MethodError(gamma_inc_inv,(a,p,q,"")))
