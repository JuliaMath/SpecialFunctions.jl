using Base.Math: @horner
using StatsFuns
#useful constants
ACC0=[5.0e-15 , 5.0e-7 , 5.0e-4]
BIG=[25.0 , 14.0, 10.0]
E0=[0.25e-3 , 0.25e-1 , 0.14]
X0=[31.0 , 17.0 , 9.7]
ALOG10 = log(10)
RT2PIN = 1.0/sqrt(2*pi)
RTPI = sqrt(pi)
EXPARG = -745.1
WK=zeros(30)
#---------COEFFICIENTS FOR MINIMAX APPROX.-------------------

A0=[-.231272501940775E-02 , -.335378520024220E-01 , -.159840143443990E+00 , -.333333333333333E+00]
B0=[.633763414209504E-06 , -.939001940478355E-05 , .239521354917408E-02 , .376245718289389E-01 , .238549219145773E+00 , .729520430331981E+00]
A1=[-.398783924370770E-05 , -.587926036018402E-03 , -.491687131726920E-02 , -.185185185184291E-02 ]
B1=[.386325038602125E-02 , .506042559238939E-01 , .283344278023803E+00 , .780110511677243E+00]
A2=[.669564126155663E-03 , .413359788442192E-02]
B2=[-.421924263980656E-03 , .650837693041777E-02 , .682034997401259E-01 , .339173452092224E+00 , .810647620703045E+00]
A3=[.810586158563431E-03 , .649434157619770E-03]
B3=[-.632276587352120E-03 , .905375887385478E-02 , .906610359762969E-01 , .406288930253881E+00 , .894800593794972E+00]
A4=[-.105014537920131E-03 , -.861888301199388E-03]
B4=[.322609381345173E-01 , .178295773562970E+00 , .591353097931237E+00 , .103151890792185E+01]
A5=[-.435211415445014E-03 , -.336806989710598E-03]
B5=[.178716720452422E+00 , .600380376956324E+00 , .108515217314415E+01]
A6=[-.182503596367782E-03 , .531279816209452E-03]
B6=[.345608222411837E+00 , .770341682526774E+00]
A7=[.443219646726422E-03 , .344430064306926E-03]
B7=[.821824741357866E+00 , .115029088777769E+01]
A8=[.878371203603888E-03 , -.686013280418038E-03]

#----------------COEFFICIENTS FOR TEMME EXPANSION------------------

D00 = -.333333333333333E+00
D0=[.833333333333333E-01 , -.148148148148148E-01 , .115740740740741E-02 , .352733686067019E-03 , -.178755144032922E-03 , .391926317852244E-04]
D10 = -.185185185185185E-02
D1=[-.347222222222222E-02 , .264550264550265E-02 , -.990226337448560E-03 , .205761316872428E-03]
D20 = .413359788359788E-02
D2=[-.268132716049383E-02 , .771604938271605E-03]
D30 = .649434156378601E-03
D3=[.229472093621399E-03 , -.469189494395256E-03]
D40 = -.861888290916712E-03
D4=[.784039221720067E-03]
D50 = -.336798553366358E-03
D5=[-.697281375836586E-04]
D60 = .531307936463992E-03
D6=[-.592166437353694E-03]
D70 = .344367606892378E-03
D80 = -.652623918595309E-03

#Computation of function 1/gamma(A+1) - 1 for -0.5<=A<=1.5
function rgamma1pm1(A::Float64)
    P=[.577215664901533E+00 , -.409078193005776E+00 , -.230975380857675E+00 , .597275330452234E-01 , .766968181649490E-02 , -.514889771323592E-02 , .589597428611429E-03]
    Q=[.100000000000000E+01 , .427569613095214E+00 , .158451672430138E+00 , .261132021441447E-01 , .423244297896961E-02]
    R=[-.422784335098468E+00 , -.771330383816272E+00 , -.244757765222226E+00 , .118378989872749E+00 , .930357293360349E-03 , -.118290993445146E-01 , .223047661158249E-02 , .266505979058923E-03 , -.132674909766242E-03]
    S1 = .273076135303957E+00
    S2 = .559398236957378E-01
    T=A
    D = A - 0.5
    if D > 0.0
        T = D - 0.5
    end
    if T == 0.0
        GAM1 = 0.0
        return GAM1
    elseif T < 0.0
        @goto ll30
    else
        @goto ll20
    end
    @label ll20
     TOP = @horner(T , P[1] , P[2] , P[3] , P[4] , P[5] , P[6] , P[7])
     BOT = @horner(T , 0.0 , Q[2] ,Q[3] ,Q[4] , Q[5]) + 1.0
     W = TOP/BOT
     if D > 0.0
        GAM1 = (T/A)*(W - 1.0)
        return GAM1
     else
        GAM1 = A*W
        return GAM1
     end
    @label ll30
     TOP = @horner(T , R[1] , R[2] , R[3] , R[4] , R[5] , R[6] , R[7] , R[8] , R[9])
     BOT = @horner(T , 0.0 , S1 , S2) + 1.0
     W = TOP/BOT
     if D > 0.0
        GAM1 = T*W/A
        return GAM1
     else
        GAM1 = A*(W + 1.0)
        return GAM1
     end
end

#Evaluation OF exp(-X)*X^A/gamma(A)
function RCOMP(A::Float64,X::Float64)
    ans = 0.0
    if X == 0.0
        return ans
    elseif A >= 20.0
        @goto twenty
    else
        T = A*log(X) - X
        if T < EXPARG
            return ans
        end
        if A >= 1.0
            ans = exp(T)/gamma(A)
            return ans
        else
            ans = (A*exp(T))*(1.0 + rgamma1pm1(A))
            return ans
        end
    end
    @label twenty
     U =X/A
     if U == 0.0
        return ans
     end
     T = (1.0/A)^2
     T1 = (((0.75*T - 1.0)*T + 3.5)*T - 105.0)/(A*1260.0)
     T1 = T1 + A*logmxp1(U)
     if T1 >= EXPARG
        ans = RT2PIN*sqrt(A)*exp(T1)
     end
     return ans
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
function gamma_p(a::Float64,x::Float64,IND::Integer)
    IOP = IND + 1
    ACC = ACC0[IOP]
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
    elseif a < 1.1
        @goto l110     
    end
    R = RCOMP(a,x)
    if R == 0.0
        ans=1.0
        return ans
    else
        @goto l170    
    end

    @label l10
     if a >= BIG[IOP]
        @goto l20
     elseif a > x || x >= X0[IOP]
        @goto l30 
     else
        TWOA = a + a
        M = trunc(Int,TWOA)
        if TWOA != Float64(M)
            @goto l30 
        end
        I = a
        I = trunc(Int,I)
        if a == Float64(I)
            @goto l140
        else
            @goto l150
        end            

     end

    @label l20
     L = x/a
     if L == 0.0
        ans=0.0
        return ans
     end
     S = 1.0 - L
     Z = -logmxp1(x)
     if Z >= 700.0/a
        ans=0.0
        return ans
     end
     Y = a*Z
     RTA = sqrt(a)
     if abs(S) <= E0[IOP]/RTA
        @goto l250
     end

     if abs(S) <= 0.4
        @goto l200
     end
    
    @label l30
      R = RCOMP(a,x)
      if R == 0.0
        if x <= a
            ans=0.0
            return ans
        else
            ans=1.0
            return ans
        end
      end 
      if x <= max(a,ALOG10)
        @goto l50
      elseif x < X0[IOP]
        @goto l170
      else
        @goto l80
      end

     #----TAYLOR SERIES FOR P/R---- 
    @label l50
     APN = a + 1.0
     T = x/APN
     WK[1] = T
     i=2
     while true
        if i > 20
            break
        end
        
        APN = APN + 1.0
        T = T*(x/APN)
        if T <= 1.0e-3
            @goto l60
        end
        WK[i] = T
        i = i + 1
     end
    i=20
    @label l60
     SUM = T
     TOL = 0.5*ACC
     while true
        APN = APN+1.0
        T = T*(x/APN)
        SUM = SUM + T
        if T < TOL
            break
        end
     end
     MAX = i - 1
     for j = 1:MAX
        i = i - 1
        SUM = SUM + WK[i]
     end
    ans = (R/a)*(1.0 + SUM)
    return ans
    
    #----ASYMPTOTIC EXPANSION-----
    @label l80
     AMN = a-1.0
     T=AMN/x
     WK[1]=T
     i=2
     while true
        if i > 20
            break
        end
        AMN = AMN-1.0
        T=T*(AMN/x)
        if abs(T) <= 1.0e-3
            @goto l90
        end
        WK[i]=T
        i=i+1
     end
     i=20
    @label l90 
     SUM = T
     while true
        if abs(T) < ACC
            @goto l100
        end
        AMN=AMN-1.0
        T=T*(AMN/x)
        SUM=SUM+T
     end
    @label l100
     MAX=i-1
     for j = 1:MAX
        i=i-1
        SUM=SUM+WK[i]
     end
    qans = (R/x)*(1.0 + SUM)
    ans = 1.0-qans
    return ans
    
    #---TAYLOR SERIES FOR P(A,X)/X**A---

    @label l110
     L=3.0
     C=x
     SUM= x/(a + 3.0)
     TOL = 3.0*ACC/(a + 1.0)
     while true
        L=L+1.0
        C=-C*(x/L)
        T=C/(a+L)
        SUM=SUM+T
        if abs(T) <= TOL 
            break
        end
     end
     J = a*x*((SUM/6.0 - 0.5/(a + 2.0))*x + 1.0/(a + 1.0))
     Z = a*log(x)
     #GAM1 = 1/gamma(a+1) - 1
     H = rgamma1pm1(a)
     #H = 1.0/gamma(a+1.0) - 1.0
     G = 1.0 + H
     if x < 0.25
        @goto l120
     end
     if a < x/2.59
        @goto l135
     else
        @goto l130
     end

    @label l120
     if Z > -.13394
        @goto l135
     end
    @label l130
     W = exp(Z)
     ans = W*G*(1.0 - J)
     return ans
    @label l135
     L = expm1(Z)
     W = 1.0-L
     qans = (W*J - L)*G - H
     if qans < 0.0
        ans=1.0
        return ans
     end
     ans = 1.0 - qans
     return ans
    
    #---FINITE SUMS FOR Q WHEN A>=1 && 2A IS INTEGER----
    @label l140
     SUM = exp(-x)
     T = SUM
     N = 1
     C=0.0
     @goto l160
    
    @label l150
     RTX = sqrt(x)
     SUM = erfc(RTX)
     T = exp(-x)/(RTPI*RTX)
     N=0
     C=-0.5
     I = trunc(Int,a - 0.5)
     
    @label l160
     
     while true
        if N == I || I == 0
            break
        end
        N = N+1
        C = C+1.0
        T = (x*T)/C
        SUM = SUM + T
     end
    @label l161
     qans = SUM
     ans = 1.0 - qans
     return ans
     
    #----CONTINUED FRACTION EXPANSION-----
    @label l170
     TOL = 4.0*ACC
     A2NM1 = 1.0
     A2N = 1.0
     B2NM1 = x
     B2N = x + (1.0 - a)
     C = 1.0
     while true
        A2NM1 = x*A2N + C*A2NM1
        B2NM1 = x*B2N + C*B2NM1
        C = C + 1.0
        T = C - a
        A2N = A2NM1 + T*A2N
        B2N = B2NM1 + T*B2N
        A2NM1 = A2NM1/B2N
        B2NM1 = B2NM1/B2N
        A2N = A2N/B2N
        B2N = 1.0
        if abs(A2N - A2NM1/B2NM1) < TOL*A2N
            break
        end
     end
     qans = R*A2N
     ans = 1.0 - qans
     return ans
    
    @label l200
     #Skipping invalid check
     C = exp(-Y)
     W = 0.5*erfcx(sqrt(Y))
     U = 1.0/a
     Z = sqrt(Z + Z) 
     if L < 1.0 
        Z=-Z
     end
     if IOP == 1
        @goto l210
     elseif IOP == 2
        @goto l220
     else
        @goto l230
     end

    @label l210
     if abs(S) <= 1.0e-3
        @goto l260
     end
     #---USING THE MINIMAX APPROXIMATIONS---
     C0 = @horner(Z , A0[4] , A0[3] , A0[2] , A0[1])/(@horner(Z , 0.0 , B0[6] , B0[5], B0[4] , B0[3] , B0[2] , B0[1]) + 1.0)
     C1 = @horner(Z , A1[4] , A1[3] , A1[2] , A1[1])/(1.0 + @horner(Z , 0.0 , B1[4] , B1[3] , B1[2] , B1[1]))
     C2 = @horner(Z , A2[2] , A2[1])/(1.0 + @horner(Z , 0.0 , B2[5] , B2[4] , B2[3] , B2[2] , B2[1]))
     C3 = @horner(Z , A3[2] , A3[1])/(1.0 + @horner(Z , 0.0 , B3[5] , B3[4] , B3[3] , B3[2] , B3[1]))
     C4 = @horner(Z , A4[2] , A4[1])/(1.0 + @horner(Z , 0.0 , B4[4] , B4[3] , B4[2] , B4[1]))
     C5 = @horner(Z , A5[2] , A5[1])/(1.0 + @horner(Z , 0.0 , B5[3] , B5[2] , B5[1]))
     C6 = @horner(Z , A6[2] , A6[1])/(1.0 + @horner(Z , 0.0 , B6[2] , B6[1]))
     C7 = @horner(Z , A7[2] , A7[1])/(1.0 + @horner(Z , 0.0 , B7[2] , B7[1]))
     C8 = @horner(Z , A8[2] , A8[1])

     T = @horner(U , C0 , C1 , C2 , C3 , C4 , C5 , C6 , C7 , C8)
     @goto l240
    
     #----TEMME EXPANSION----
    @label l220
     C0 = @horner(Z , 0.0 , D0[1] , D0[2] , D0[3] , D0[4] , D0[5] , D0[6]) + D00
     C1 = @horner(Z , 0.0 , D1[1] , D1[2] , D1[3] , D1[4]) + D10

     T = @horner(U , 0.0 , C1 , C2) + C0
     @goto l240
     
    @label l230
     T = @horner(Z , 0.0 , D0[1] , D0[2] , D0[3]) + D00
    @label l240
     if L < 1.0
        @goto l241
     end
     qans = C*(W + RT2PIN*T/RTA)
     ans = 1.0 - qans
     return ans
    @label l241
     ans = C*(W - RT2PIN*T/RTA)
     qans = 1.0 - ans
     return ans
    
    #----TEMME EXPANSION FOR L = 1----
    @label l250
     #Skipping error return
     C = 1.0 - Y
     W = (0.5 - sqrt(Y)*(0.5 + (0.5 - Y/3.0))/RTPI)/C
     U = 1.0/a
     Z = sqrt(Z + Z)
     if L < 1.0
        Z = -Z
     end
     if IOP == 1
        @goto l260
     elseif IOP == 2
        @goto l270
     else
        @goto l280
     end
    
    @label l260
     C0 = @horner(Z , D00 , D0[1] , D0[2] , D0[3])
     C1 = @horner(Z , D10 , D1[1] , D1[2] , D1[3])
     C2 = @horner(Z , D20 , D2[1] , D2[2])
     C3 = @horner(Z , D30 , D3[1] , D3[2])
     C4 = @horner(Z , D40 , D4[1])
     C5 = @horner(Z , D50 , D5[1])
     C6 = @horner(Z , D60 , D6[1])

     T = @horner(U , C0 , C1 , C2 , C3 , C4 , C5 , C6 , D70 , D80)
     @goto l240rep
    
    @label l270
     C0 = @horner(D00 , D0[1] , D0[2])
     C1 = @horner(D10 , D1[1])
     T = @horner(U , C0 , C1 , D20)
     @goto l240rep
    
    @label l280
     T = @horner(Z , D00 , D0[1])
     @goto l240rep
    
    @label l240rep
     if L < 1.0
        @goto l241rep
     end
     qans = C*(W + RT2PIN*T/RTA)
     ans = 1.0 - qans
     return ans
    @label l241rep
     ans = C*(W - RT2PIN*T/RTA)
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
function gamma_q(a::Float64,x::Float64,IND::Integer)
    qans = 1.0 - gamma_p(a,x,IND)
    return qans
end

for f in (:gamma_p,:gamma_q)
    @eval begin
        $f(a::Float32,x::Float32,IND::Integer) = Float32($f(Float64(a),Float64(x),IND))
        $f(a::Float16,x::Float16,IND::Integer) = Float16($f(Float64(a),Float64(x),IND))
        $f(a::Real,x::Real,IND::Integer) = ($f(float(a),float(x),IND))
        $f(a::Integer,x::Integer,IND::Integer) = $f(Float64(a),Float64(x),IND)
        $f(a::AbstractFloat,x::AbstractFloat,IND::Integer) = throw(MethodError($f,(a,x,IND,"")))        
    end
end
