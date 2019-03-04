using Base.Math: @horner
#useful constants
const acc0 = [5.0e-15 , 5.0e-7 , 5.0e-4] #accuracy options
const big1 = [25.0 , 14.0, 10.0]
const e0 = [0.25e-3 , 0.25e-1 , 0.14]
const x0=[31.0 , 17.0 , 9.7]
const alog10 = log(10)
const rt2pin = 1.0/sqrt(2*pi)
const rtpi = sqrt(pi) 
const exparg = -745.1

wk = zeros(30)
#---------COEFFICIENTS FOR MINIMAX APPROX.-------------------

a0=[-.231272501940775E-02 , -.335378520024220E-01 , -.159840143443990E+00 , -.333333333333333E+00]
b0=[.633763414209504E-06 , -.939001940478355E-05 , .239521354917408E-02 , .376245718289389E-01 , .238549219145773E+00 , .729520430331981E+00]
a1=[-.398783924370770E-05 , -.587926036018402E-03 , -.491687131726920E-02 , -.185185185184291E-02 ]
b1=[.386325038602125E-02 , .506042559238939E-01 , .283344278023803E+00 , .780110511677243E+00]
a2=[.669564126155663E-03 , .413359788442192E-02]
b2=[-.421924263980656E-03 , .650837693041777E-02 , .682034997401259E-01 , .339173452092224E+00 , .810647620703045E+00]
a3=[.810586158563431E-03 , .649434157619770E-03]
b3=[-.632276587352120E-03 , .905375887385478E-02 , .906610359762969E-01 , .406288930253881E+00 , .894800593794972E+00]
a4=[-.105014537920131E-03 , -.861888301199388E-03]
b4=[.322609381345173E-01 , .178295773562970E+00 , .591353097931237E+00 , .103151890792185E+01]
a5=[-.435211415445014E-03 , -.336806989710598E-03]
b5=[.178716720452422E+00 , .600380376956324E+00 , .108515217314415E+01]
a6=[-.182503596367782E-03 , .531279816209452E-03]
b6=[.345608222411837E+00 , .770341682526774E+00]
a7=[.443219646726422E-03 , .344430064306926E-03]
b7=[.821824741357866E+00 , .115029088777769E+01]
a8=[.878371203603888E-03 , -.686013280418038E-03]

#----------------COEFFICIENTS FOR TEMME EXPANSION------------------

d00 = -.333333333333333E+00
d0=[.833333333333333E-01 , -.148148148148148E-01 , .115740740740741E-02 , .352733686067019E-03 , -.178755144032922E-03 , .391926317852244E-04]
d10 = -.185185185185185E-02
d1=[-.347222222222222E-02 , .264550264550265E-02 , -.990226337448560E-03 , .205761316872428E-03]
d20 = .413359788359788E-02
d2=[-.268132716049383E-02 , .771604938271605E-03]
d30 = .649434156378601E-03
d3=[.229472093621399E-03 , -.469189494395256E-03]
d40 = -.861888290916712E-03
d4=[.784039221720067E-03]
d50 = -.336798553366358E-03
d5=[-.697281375836586E-04]
d60 = .531307936463992E-03
d6=[-.592166437353694E-03]
d70 = .344367606892378E-03
d80 = -.652623918595309E-03

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
"""
function rgamma1pm1(a::Float64)
    t=a
    d = a - 0.5
    if d > 0.0
        t = d - 0.5
    end
    if t == 0.0
        return 0.0
    elseif t < 0.0
     top = @horner(t , -.422784335098468E+00 , -.771330383816272E+00 , -.244757765222226E+00 , .118378989872749E+00 , .930357293360349E-03 , -.118290993445146E-01 , .223047661158249E-02 , .266505979058923E-03 , -.132674909766242E-03)
     bot = @horner(t , 1.0 , .273076135303957E+00 , .559398236957378E-01)
     w = top/bot
     if d > 0.0
        return t*w/a
     else
        return a*(w + 1.0)
     end
    else
     top = @horner(t , .577215664901533E+00 , -.409078193005776E+00 , -.230975380857675E+00 , .597275330452234E-01 , .766968181649490E-02 , -.514889771323592E-02 , .589597428611429E-03)
     bot = @horner(t , 1.0 , .427569613095214E+00 , .158451672430138E+00 , .261132021441447E-01 , .423244297896961E-02) 
     w = top/bot
     if d > 0.0
        return (t/a)*(w - 1.0)
     else
        return a*w
     end
    end
end
"""
    exp1xgamma1(a,x)

Evaluation of exp(-x)*x^a/gamma(a) : ``1/\\Gamma(a) e^{-x} x^{a}``
"""
function exp1xgamma1(a::Float64,x::Float64)
    ans = 0.0
    if x == 0.0
        return ans
    elseif a >= 20.0
     u =x/a
     if u == 0.0
        return ans
     end
     t = (1.0/a)^2
     t1 = (((0.75*t - 1.0)*t + 3.5)*t - 105.0)/(a*1260.0)
     t1 = t1 + a*logmxp1(u)
     if t1 >= exparg
        ans = rt2pin*sqrt(a)*exp(t1)
     end
     return ans
    else
        t = a*log(x) - x
        if t < exparg
            return ans
        end
        if a >= 1.0
            ans = exp(t)/gamma(a)
            return ans
        else
            ans = (a*exp(t))*(1.0 + rgamma1pm1(a))
            return ans
        end
    end
end   

# Reference : 'Computation of the incomplete gamma function ratios and their inverse' by Armido R DiDonato , Alfred H Morris.
# Published in Journal: ACM Transactions on Mathematical Software (TOMS)
# Volume 12 Issue 4, Dec. 1986 Pages 377-393
# doi>10.1145/22721.23109

"""
    gamma_p(a,x,IND)
    
    DLMF: https://dlmf.nist.gov/8.2#E4 , https://dlmf.nist.gov/8.2#E5
    Wiki: https://en.wikipedia.org/wiki/Incomplete_gamma_function
    IND --> Accuracy desired ; IND=0 means 14 significant digits accuracy , IND=1 means 6 significant digit and IND=2 means only 3 digit accuracy suffices.
    gamma_p(a,x) or P(a,x) is the Incomplete gamma function ratio given by : ``1/\\Gamma (a) \\int_{0}^{x} e^{-t}t^{a-1} dt``
"""
function gamma_p(a::Float64,x::Float64,ind::Integer)
    iop = ind + 1
    acc = acc0[iop]
    if a<0.0 || x<0.0
        ans=2.0
        return ans
    elseif a==0.0 && x==0.0
        ans=2.0
        return ans
    elseif a*x==0.0 
        if x<=a
            ans=0.0
            return ans
        else
            ans=1.0
            return ans
        end               
    end
    
    if a >= 1.0
        @goto l10
    elseif a == 0.5
        @goto l320
    elseif x < 1.1
        @goto l110     
    end
    r = exp1xgamma1(a,x)
    if r == 0.0
        ans=1.0
        return ans
    else
        @goto l170    
    end

    @label l10
     if a >= big1[iop]
        @goto l20
     elseif a > x || x >= x0[iop]
        @goto l30 
     else
        twoa = a + a
        m = trunc(Int,twoa)
        if twoa != Float64(m)
            @goto l30 
        end
        i = a
        i = trunc(Int,i)
        if a == Float64(i)
            @goto l140
        else
            @goto l150
        end            

     end

    @label l20
     l = x/a
     if l == 0.0
        ans=0.0
        return ans
     end
     s = 1.0 - l
     z = -logmxp1(x)
     if z >= 700.0/a
        ans=0.0
        return ans
     end
     y = a*z
     rta = sqrt(a)
     if abs(s) <= e0[iop]/rta
        @goto l250
     end

     if abs(s) <= 0.4
        @goto l200
     end
    
    @label l30
      r = exp1xgamma1(a,x)
      if r == 0.0
        if x <= a
            ans=0.0
            return ans
        else
            ans=1.0
            return ans
        end
      end 
      if x <= max(a,alog10)
        @goto l50
      elseif x < x0[iop]
        @goto l170
      else
        @goto l80
      end

     #----TAYLOR SERIES FOR P/R---- 
    @label l50
     apn = a + 1.0
     t = x/apn
     wk[1] = t
     loop=2
     while true
        if loop > 20
            break
        end
        
        apn = apn + 1.0
        t = t*(x/apn)
        if t <= 1.0e-3
            @goto l60
        end
        wk[loop] = t
        loop = loop + 1
     end
    loop=20
    @label l60
     sum = t
     tol = 0.5*acc #tolerance
     while true
        apn = apn+1.0
        t = t*(x/apn)
        sum = sum + t
        if t <= tol
            break
        end
     end
     maxd = loop - 1
     for j = 1:maxd
        loop = loop - 1
        sum = sum + wk[loop]
     end
    ans = (r/a)*(1.0 + sum)
    return ans
    
    #----ASYMPTOTIC EXPANSION-----
    @label l80
     amn = a-1.0
     t=amn/x
     wk[1]=t
     loop=2
     while true
        if loop > 20
            break
        end
        amn = amn-1.0
        t=t*(amn/x)
        if abs(t) <= 1.0e-3
            @goto l90
        end
        wk[loop]=t
        loop=loop+1
     end
     loop=20
    @label l90 
     sum = t
     while true
        if abs(t) < acc
            @goto l100
        end
        amn=amn-1.0
        t=t*(amn/x)
        sum=sum+t
     end
    @label l100
     maxd=loop-1
     for j = 1:maxd
        loop=loop-1
        sum=sum+wk[loop]
     end
    qans = (r/x)*(1.0 + sum)
    ans = 1.0-qans
    return ans
    
    #---TAYLOR SERIES FOR P(A,X)/X**A---

    @label l110
     l=3.0
     c=x
     sum= x/(a + 3.0)
     tol = 3.0*acc/(a + 1.0)
     while true
        l=l+1.0
        c=-c*(x/l)
        t=c/(a+l)
        sum=sum+t
        if abs(t) <= tol 
            break
        end
     end
     temp = a*x*((sum/6.0 - 0.5/(a + 2.0))*x + 1.0/(a + 1.0))
     z = a*log(x)
     #GAM1 = 1/gamma(a+1) - 1
     h = rgamma1pm1(a)
     #H = 1.0/gamma(a+1.0) - 1.0
     g = 1.0 + h
     if x < 0.25
        @goto l120
     end
     if a < x/2.59
        @goto l135
     else
        @goto l130
     end

    @label l120
     if z > -.13394
        @goto l135
     end
    @label l130
     w = exp(z)
     ans = w*g*(1.0 - temp)
     return ans
    @label l135
     l = expm1(z)
     w = 1.0-l
     qans = (w*temp - l)*g - h
     if qans < 0.0
        ans=1.0
        return ans
     end
     ans = 1.0 - qans
     return ans
    
    #---FINITE SUMS FOR Q WHEN A>=1 && 2A IS INTEGER----
    @label l140
     sum = exp(-x)
     t = sum
     N = 1
     c=0.0
     @goto l160
    
    @label l150
     rtx = sqrt(x)
     sum = erfc(rtx)
     t = exp(-x)/(rtpi*rtx)
     N=0
     c=-0.5
     i = trunc(Int,a - 0.5)
     
    @label l160
     
     while true
        if N == i || i == 0
            break
        end
        N = N+1
        c = c+1.0
        t = (x*t)/c
        sum = sum + t
     end
    @label l161
     qans = sum
     ans = 1.0 - qans
     return ans
     
    #----CONTINUED FRACTION EXPANSION-----
    @label l170
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
     qans = r*a2n
     ans = 1.0 - qans
     return ans
    
    @label l200
     #Skipping invalid check
     c = exp(-y)
     w = 0.5*erfcx(sqrt(y))
     u = 1.0/a
     z = sqrt(z + z) 
     if l < 1.0 
        z=-z
     end
     if iop == 1
        @goto l210
     elseif iop == 2
        @goto l220
     else
        @goto l230
     end

    @label l210
     if abs(s) <= 1.0e-3
        @goto l260
     end
     #---USING THE MINIMAX APPROXIMATIONS---
     c0 = @horner(z , a0[4] , a0[3] , a0[2] , a0[1])/(@horner(z , 1.0 , b0[6] , b0[5], b0[4] , b0[3] , b0[2] , b0[1]) )
     c1 = @horner(z , a1[4] , a1[3] , a1[2] , a1[1])/(@horner(z , 1.0 , b1[4] , b1[3] , b1[2] , b1[1]))
     c2 = @horner(z , a2[2] , a2[1])/(@horner(z , 1.0 , b2[5] , b2[4] , b2[3] , b2[2] , b2[1]))
     c3 = @horner(z , A3[2] , A3[1])/(@horner(z , 1.0 , b3[5] , b3[4] , b3[3] , b3[2] , b3[1]))
     c4 = @horner(z , a4[2] , a4[1])/(@horner(z , 1.0 , b4[4] , b4[3] , b4[2] , b4[1]))
     c5 = @horner(z , a5[2] , a5[1])/(@horner(z , 1.0 , b5[3] , b5[2] , b5[1]))
     c6 = @horner(z , a6[2] , a6[1])/(@horner(Z , 1.0 , b6[2] , b6[1]))
     c7 = @horner(z , a7[2] , a7[1])/(@horner(Z , 1.0 , b7[2] , b7[1]))
     c8 = @horner(z , a8[2] , a8[1])

     t = @horner(u , c0 , c1 , c2 , c3 , c4 , c5 , c6 , c7 , c8)
     @goto l240
    
     #----TEMME EXPANSION----
    @label l220
     c0 = @horner(z , d00 , d0[1] , d0[2] , d0[3] , d0[4] , d0[5] , d0[6]) 
     c1 = @horner(Z , d10 , d1[1] , d1[2] , d1[3] , d1[4]) 

     t = @horner(u , c0 , c1 , c2)
     @goto l240
     
    @label l230
     t = @horner(z , d00 , d0[1] , d0[2] , d0[3])
    @label l240
     if l < 1.0
        @goto l241
     end
     qans = c*(w + rt2pin*t/rta)
     ans = 1.0 - qans
     return ans
    @label l241
     ans = c*(w - rt2pin*t/rta)
     qans = 1.0 - ans
     return ans
    
    #----TEMME EXPANSION FOR L = 1----
    @label l250
     #Skipping error return
     c = 1.0 - y
     w = (0.5 - sqrt(y)*(0.5 + (0.5 - y/3.0))/rtpi)/c
     u = 1.0/a
     z = sqrt(z + z)
     if l < 1.0
        z = -z
     end
     if iop == 1
        @goto l260
     elseif iop == 2
        @goto l270
     else
        @goto l280
     end
    
    @label l260
     c0 = @horner(z , d00 , d0[1] , d0[2] , d0[3])
     c1 = @horner(z , d10 , d1[1] , d1[2] , d1[3])
     c2 = @horner(z , d20 , d2[1] , d2[2])
     c3 = @horner(z , d30 , d3[1] , d3[2])
     c4 = @horner(z , d40 , d4[1])
     c5 = @horner(z , d50 , d5[1])
     c6 = @horner(z , d60 , d6[1])

     t = @horner(u , c0 , c1 , c2 , c3 , c4 , c5 , c6 , d70 , d80)
     @goto l240rep
    
    @label l270
     c0 = @horner(d00 , d0[1] , d0[2])
     c1 = @horner(d10 , d1[1])
     t = @horner(u , c0 , c1 , d20)
     @goto l240rep
    
    @label l280
     t = @horner(z , d00 , d0[1])
     @goto l240rep
    
    @label l240rep
     if l < 1.0
        @goto l241rep
     end
     qans = c*(w + rt2pin*t/rta)
     ans = 1.0 - qans
     return ans
    @label l241rep
     ans = c*(w - rt2pin*t/rta)
     qans = 1.0 - ans
     return ans
    
    @label l320
     if x >= 0.25
        @goto l321
     end
     ans = erf(sqrt(x))
     return ans
    @label l321
     qans = erfc(sqrt(x))
     ans = 1.0 - qans
     return ans
end

# Reference : 'Computation of the incomplete gamma function ratios and their inverse' by Armido R DiDonato , Alfred H Morris.
# Published in Journal: ACM Transactions on Mathematical Software (TOMS)
# Volume 12 Issue 4, Dec. 1986 Pages 377-393
# doi>10.1145/22721.23109

"""
    gamma_q(a,x,IND)
    
    DLMF: https://dlmf.nist.gov/8.2#E4 , https://dlmf.nist.gov/8.2#E5
    Wiki: https://en.wikipedia.org/wiki/Incomplete_gamma_function
    IND --> Accuracy desired ; IND=0 means 14 significant digits accuracy , IND=1 means 6 significant digit and IND=2 means only 3 digit accuracy suffices.
    gamma_q(a,x) or Q(a,x) is the Incomplete gamma function ratio given by : 1 - P(a,x) ->  ``1/\\Gamma (a) \\int_{x}^{\\infty} e^{-t}t^{a-1} dt``
"""
function gamma_q(a::Float64,x::Float64,ind::Integer)
    qans = 1.0 - gamma_p(a,x,ind)
    return qans
end

for f in (:gamma_p,:gamma_q)
    @eval begin
        $f(a::Float32,x::Float32,ind::Integer) = Float32($f(Float64(a),Float64(x),ind))
        $f(a::Float16,x::Float16,ind::Integer) = Float16($f(Float64(a),Float64(x),ind))
        $f(a::Real,x::Real,ind::Integer) = ($f(float(a),float(x),ind))
        $f(a::Integer,x::Integer,ind::Integer) = $f(Float64(a),Float64(x),ind)
        $f(a::AbstractFloat,x::AbstractFloat,ind::Integer) = throw(MethodError($f,(a,x,ind,"")))        
    end
end
