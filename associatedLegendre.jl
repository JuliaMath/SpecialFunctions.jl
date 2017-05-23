
# solutions, w, of the differential equation are 
#    associated Legendre functions   
# (1-x*x)d/dx(dw/dx)⁢-2⁢x⁢(dw/dx)+(n(n+1)-m*m/(1-x*x))w=0
#for m =0 they are a.k.a. Legendre polynomials
#n is the degree and m is the order 
#
#There are more than a handful of programs. Some fast and inaccurate
#and some slow and more accurate; pick your poison (trade off) .
#
#       SSLegendreP2 and 3 : The finite sum is arranged so that 
# each term is a multiple of the preceeding term. SLegendreP2 and 3
#  are very slow and used for testing.
#
#       The recursion equation 8.5.3 Abramowitz & Stegun (AS) Handbook of
# Mathematical Functions  used in LegendreP2 is
#  (n-m+1)P(n+1,m,z)=(2n+1)zP(n,m,z)- (n+m)P(n-1,m,z)
#
#
#        ZZLegendreP2 and 3  uses a different recursion
#  P(n,m+2,x)+2(m+1)x(s(1-x*x))^(-1/2)P(n,m+1,x)+s(n-m)(n+m+1)P(n,m,x)=0
#from http://dlmf.nist.gov/14.10  equation  14.10.1  where s=sign(1-x*x)
#refer to Lebedev, Special Functions, Dover Pub. eqn. 7.12.8 page 194
#             See function BenchmarkZZLegemdreP2() for speed.
#
#These recursions are valid for Legendre functions of the first kind 
#(P(n,m,z))and second kind Q(n,m,z) , 
#(for types 1,2,3 ,definition from Mathematica documentation
# type 1 (?) is m=0, type 2 is -1<=x<=1 and type 3 is z> 1
# http://mpmath.org/doc/0.18/functions/orthogonal.html
#Thus we define functions  LegendreP2,....P3,....Q2,...Q3  
# for n and m positive integers 0<=m<=n
# refer to http://dlmf.nist.gov/14.6
#P(n,m,z) and Q(n,m,⁡z) (z> 1) are often referred to as the prolate
#spheroidal harmonics of the first and second kinds, respectively
#  We can define spherical harmonics (there are other definitions)
# Y(n,m,theta,phi)= LegendreP2(n,|m|,theta)*exp(im*m*phi)*((-)^(|m|-m)/2))
# *sqrt((2n+1)*factorial(n-|m|)/(4*pi*factorial(n+|m|))) 
# refer to https://en.wikipedia.org/wiki/Spherical_harmonics






"""
    LegendreP2(n::Integer,m::Integer,x::Number)
    associated Legendre function of first kind type 2 (-1.<=x<=1.).
    

"""

function LegendreP2(n::Integer,m::Integer,x::Number) 
  
     
    # 0 <= m  <= n
    # -1<=x <= 1
    M = (2*m -1) # M must be odd
    # (2m-1)!!= (2m)!/( m! 2^m)   double factorial 
    dblfac=1
    for j=1:M
        if  iseven(j)
            continue
        end    
        dblfac=j*dblfac  
    end
   
    pj2= ((-1)^m)*dblfac*(one(x)-x*x)^(m/2)   
    pj1=x*(2.*m+1.)*pj2  
    if n == m
        return pj2
    elseif n == m+1
        return pj1
    end
    

    for j = m+2 :n 
        pjj=(x*(2.*j-1.)*pj1 - pj2*(j +m-1.)) /(j-m)
        pj2=pj1
        pj1=pjj
    end
return pj1 
end




#checking LegendreP2 against scipy.special
#   for n< 18 the relative error is less than 2e-14
# the error depends on the n amd m values and x values to some extent

using PyCall

@pyimport scipy.special as s

"""
    function testLegendreP2()
    relative error LegendreP2 versus scipy special
"""


function testLegendreP2()
    x=.6
    rerr = 2e-14
    for n=0:17
        for m=0:n      #m=rand(0:n)
            a=LegendreP2(n,m,x)
            b=s.lpmv(m,n,x)
            df=norm(a-b)/norm(b) 
            if ( df > rerr)
            #if (norm(a-b)/norm(b)) > rerr)
            println(n," ", m," ", a," ",b," ",df)
            #@test_approx_eq_eps a b 1e-14
            end
        end
    
    end

end
 
testLegendreP2()

"""
    SLegendreP2(n,m,x)
    Legendre function 1st kind type 2 , -1.<=x<=1.
    series (to be optimized for in SSLegendreP2)
"""

function SLegendreP2(n,m,x)#n,m integers 

    v=(big(1.)-x*x)^(-m*big(1)/2)
    v=v*factorial(big(n+m))/(((big(2))^n)* factorial(big(n)));
    term =  big(0.);
    #if iseven(n+m); MXP=Int(trunc((n+m)/2));end;
    #if isodd(n+m);MXP=Int(trunc((n+m-1)/2));end;
    #MXP=Int(floor((n+m)/2));
    
    MXP=Int(trunc(((n+m)/2)))

    
    for p=0:MXP
                                             
        term= term + 
        ((-big(1.))^p)*binomial(big(2*(n-p)),big(n-m))*
        (x^(big(1)*(n + m -2*p)))*
        binomial(big(n),big(p)) 
        
    end
    return term * v
end 

# I derived this finite sum using  8.6.6,8.6.18,8.2.5 (AS),
#binomial expansion and differentiating





"""
    function testLegendreP2F()
    compare LegendreP2 against Float64(SLegendreP2)
"""
function testLegendreP2F()
#=
   checking LegendreP2 against Float64SLegendreP2
   for n < 18 the relative error is less than 7e-15
   the error depends on the n amd m values and 
        x values to some extent
 =#
    

    x=.6
    rerr = 7e-15
    for n=1:17
        for  m=0:n   #rand(1:n)
            a=LegendreP2(n,m,x)
            b=Float64(SLegendreP2(n,m,big(x)))
            df=norm(a-b)/norm(b) 
            if ( df > rerr)
            println(n," ", m," ", a," ",b," ",df)
            end
        end
end
end
testLegendreP2F()



    
"""
    SSLegendreP2(n::Integer,m::Integer,x::Number)
associated Legendre function 1st kind type 2 , -1.<=x<=1.
each term in sum is multiple of previous, from SLegendreP2 


"""

function SSLegendreP2(n::Integer,m::Integer,x::Number)#n,m integers 
   
    M = (2*n -1) # M must be odd
    # (2m-1)!!= (2m)!/( m! 2^m)   double factorial 
    dblfac=1
    for j=1:M
        if  iseven(j)
            continue
        end    
        dblfac=j*dblfac  
    end
    
    MXP=Int(trunc(((n+m)/2)))
    sum= (x^(n+m))*dblfac/(((one(x)-x*x)^(m/2))*prod(1:n-m))
    prev=sum
    for p=1:MXP
        term=(-prev)*(n-p+1)*(n+m-2*p+1)*(n+m-2*p+2
        )/(p*(2*n-2*p+1)*(2*n-2*p+2)*x*x)
        sum=sum+term
        prev=term
    end
    return sum
end



#checking SSLegendreP2 against scipy.special
#   for n< 17 the relative error is less than 3e-11
# the error depends on the n amd m values and x values to some extent

using PyCall

@pyimport scipy.special as s

"""
    function testSSLegendreP2()
    compare SSLegendreP2 with scipy special
"""
function testSSLegendreP2()
    x=.6
    rerr = 3e-11
    for n=1:17
        for m=0:n      #m=rand(1:n)
            a=SSLegendreP2(n,m,x)
            b=s.lpmv(m,n,x)
            df=norm(a-b)/norm(b) 
            if ( df > rerr)
            #if (norm(a-b)/norm(b)) > rerr)
            println(n," ", m," ", a," ",b," ",df)
            #@test_approx_eq_eps a b 1e-14
            end
        end
    
    end

end
testSSLegendreP2()

#checking SSLegendreP2 against Float64(SLegendreP2)
#   for n< 17 the relative error is less than 3e-11
# the error depends on the n amd m values and x values to some extent
"""
    function testSSLegendreP2F()
    compare SSLegendreP2 to float64SLegendreP2
"""
function testSSLegendreP2F()
    x=.6
    rerr = 3e-11
    for n=0:17
        for m=0:n    #m=rand(0:n)
            a=SSLegendreP2(n,m,x)
            b=Float64(SLegendreP2(n,m,big(x)))
            df=norm(a-b)/norm(b) 
            if ( df > rerr)
            println(n," ", m," ", a," ",b," ",df)
            end
        end
    end
end
testSSLegendreP2F()


"""
SLegendreP3(n::Integer,m::Integer,z::Number)
    associated Legendre function z>1 1st kind type 3
    series expansion
"""
function SLegendreP3(n::Integer,m::Integer,z::Number)
#function mylegenp3(n,m,z)#n,m integers 
    v=(z*z - big(1.))^(-m*big(1)/2)
  
    v=v*factorial(big(n+m))/(((big(2))^n)* factorial(big(n)));
    term =  big(0.);
         #if iseven(n+m); MXP=Int(trunc((n+m)/2));end;
         #if isodd(n+m);MXP=Int(trunc((n+m-1)/2));end;
         #MXP=Int(floor((n+m)/2));
    MXP=Int(trunc(((n+m)/2)))
    
    for p=0:MXP
                                                          
        term= term + 
        ((-big(1.))^p)*binomial(big(2*(n-p)),big(n-m))*
        (z^(big(1)*(n + m -2*p)))*
        binomial(big(n),big(p)) 
        
    end
    return term * v
end   

#print mylegenp3(2,1,2.)#agrees to 16 decimals
#print mylegenp3(2,1,mpf(2))#agrees
#print mylegenp3(2,1,mpf('2.'))#agrees 10.392304845413263761
#println(mylegenp3(2,1,big(2.)))  AGREES with mpmath

# I derived this finite sum using  8.6.6,8.6.18,8.2.5 (AS),
#binomial expansion and differentiating


function timeS2()
    for i=1:10
        SLegendreP2(100,50,big(.6))
    end
end   
TS2= @elapsed(timeS2())
println("time SLegendreP2(100,50,big(.6))  $TS2")

function timeP2()
    for i=1:10
        SSLegendreP2(100,50,big(.6))
    end
end   
TP2= @elapsed(timeP2())
println("time SSLegendreP2(100,50,big(.6))  $TP2")
#time SLegendreP2(100,50,big(.6))  0.0085234
#time SSLegendreP2(100,50,big(.6))  0.004890456
# the series is double speed when each term is a multiple of previous term



"""
    SSLegendreP3(n::Integer,m::Integer,x::Number)
    associated Legendre function 1st kind type 3 , x >1.
    each term in sum is multiple of previous. 

"""

function SSLegendreP3(n::Integer,m::Integer,x::Number)#n,m integers 
    
   
    M = (2*n -1) # M must be odd
    # (2m-1)!!= (2m)!/( m! 2^m)   double factorial 
    dblfac=1
    for j=1:M
        if  iseven(j)
            continue
        end    
        dblfac=j*dblfac  
    end
    
    MXP=Int(trunc(((n+m)/2)))
    sum= (x^(n+m))*dblfac/(((x*x - one(x))^(m/2))*prod(1:n-m))
    prev=sum
    for p=1:MXP
        term=(-prev)*(n-p+1)*(n+m-2*p+1)*(n+m-2*p+2
        )/(p*(2*n-2*p+1)*(2*n-2*p+2)*x*x)
        sum=sum+term
        prev=term
    end
    return sum
end


   

#checking SSLegendreP3 against Float64(SLegendreP3)
#   for n< 17 the relative error is less than 4e-14
# the error depends on the n amd m values and x values to some extent
"""
    function testSSLegendreP3F()
    compare SSLegendreP3 with Float64 SLegendreP3
"""
function testSSLegendreP3F()
    x=2.
    rerr = 4e-14
    for n=0:17
        for m=rand(0:n)
            a=SSLegendreP3(n,m,x)
            b=Float64(SLegendreP3(n,m,big(x)))
            df=norm(a-b)/norm(b) 
            if ( df > rerr)
            println(n," ", m," ", a," ",b," ",df)
            end
        end
end
end
testSSLegendreP3F()

# based on code published at
# i. a. selezneva, yu. l. ratis, e. hernández, j. pérez-quiles 
#and p. fernández de córdoba:
#a code to calculate high order legendre polynomials
#rev. acad. colomb. cienc.: volumen xxxvii, número 145 - diciembre 2013
#www.scielo.org.co/pdf/racefn/v37n145/v37n145a09.pdf   
    #Selezneva.....pdf    
       
"""
ZZLegendreP2(n::Integer,m::Integer,z::Number)
    associated Legendre function First kind type 2 |x|<=1.
    uses different three term recursion
"""

function ZZLegendreP2(n::Integer,m::Integer,z::Number)    
    nz=1
    pnm = zeros(typeof(z),2*n+1)  
    fac = prod(2.:n)
    sqz2 = sqrt((one(z)-z.*z))
    hsqz2 = 0.5*sqz2
    ihsqz2 = z./hsqz2
    if(n==0)
        pnm[1]=one(z)
        return ((-one(z))^m)*pnm[n+1+m] #pnm
    end
    if(n==1)
        pnm[1]=-.5*sqz2
        pnm[2]=z 
        pnm[3]=sqz2
        return ((-one(z))^m)*pnm[n+1+m]  #pnm
    end
    pnm[1] = (1-2*abs(n-2*floor(n/2)))*hsqz2.^n/fac
    pnm[2] = -pnm[1]*n.*ihsqz2
    for mr=1:2*n-1
        pnm[mr+2]=(mr-n).*ihsqz2.*pnm[mr+1]-(2*n-mr+1)*mr*pnm[mr]
    end
    return ((-one(z))^m)*pnm[n+1+m]
end
 

 

#Benchmarking ZZLegendreP2 against scipy.special

using PyCall

@pyimport scipy.special as s


function timePY()
    for i=1:100
        s.lpmv(15,30,.6) 
    end
end 

function timeZ2()
    for i=1:100
        ZZLegendreP2(30,15,.6)
    end
end 
"""
    function BenchmarkZZLegendreP2()
    timing for ZZLegendreP2 LegendreP2 
    more than five times as fast 
    compared to s.lpmv from scipy.special
    speed may depend on n,m,x values
"""

function BenchmarkZZLegendreP2()

TPY= @elapsed(timePY())
println("time s.lpmv(15,30,.6)   $TPY")


TZ2= @elapsed(timeZ2())
println("time ZZLegendreP2(30,15,.6)  $TZ2")

#time s.lpmv(15,30,.6)   0.000674941
#time ZZLegendreP2(30,15,.6)  0.000112901
    
end
BenchmarkZZLegendreP2()

#checking ZZLegendreP2 against scipy special
#   for n< 51 the relative error is less than 6e-7
# the error depends on the n amd m values and x values to some extent
#  
"""
    function testZZLegendreP2()
compare ZZLegendreP2 with scipy special
"""
function testZZLegendreP2()
    x=.6
    rerr = 6e-7
    for n=1:50
        for m= 0:n  #rand(1:n)
            a=ZZLegendreP2(n,m,x)
            b=s.lpmv(m,n,x)
            df=norm(a-b)/norm(b) 
            if ( df > rerr)
            println(n," ", m," ", a," ",b," ",df)
            end
        end
    end
end
testZZLegendreP2()


#checking ZZLegendreP2 against Float64(SLegendreP2)
#   for n< 51 the relative error is less than 6e-7
# the error depends on the n amd m values and x values to some extent
#  largest errors for n=m or n= m+1 (could be fixed?)
"""
    function testZZLegendreP2F()
    compare ZZLegendreP2 with Float64 SLegendreP2
"""
function testZZLegendreP2F()
    x=.6
    rerr = 6e-7
    for n=1:50
        for m= 0:n  #rand(1:n)
            a=ZZLegendreP2(n,m,x)
            b=Float64(SLegendreP2(n,m,big(x)))
            df=norm(a-b)/norm(b) 
            if ( df > rerr)
            println(n," ", m," ", a," ",b," ",df)
            end
        end
    end
end
testZZLegendreP2F()


# derived from work published at 
# i. a. selezneva, yu. l. ratis, e. hernández, j. pérez-quiles 
#and p. fernández de córdoba:
#a code to calculate high order legendre polynomials
#rev. acad. colomb. cienc.: volumen xxxvii, número 145 - diciembre 2013
#www.scielo.org.co/pdf/racefn/v37n145/v37n145a09.pdf   
    #Selezneva.....pdf    
       

"""
ZZLegendreP3(n::Integer,m::Integer,z::Number)
    associated Legendre function First kind type 3 ( z > 1.)
    uses different three term recursion 
"""

function ZZLegendreP3(n::Integer,m::Integer,z::Number)  
   
    nz=1
    pnm = zeros(typeof(z),2*n+1)#was nz,n+1
    fac = prod(2.:n)
    sqz2 = sqrt((z.*z - one(z)))
    hsqz2 = 0.5*sqz2
    ihsqz2 = z./hsqz2
    if(n==0)
        pnm[1]=one(z)
        return  ((-one(z))^m)*pnm[n+1+m]  
    end
    if(n==1)
        pnm[1]=-.5*sqz2
        pnm[2]=z
        pnm[3]=-sqz2 # sign
        return   ((-one(z))^m)*pnm[n+1+m] 
    end
    pnm[1] = (1-2*abs(n-2*floor(n/2)))*hsqz2.^n/fac
    pnm[2] = -pnm[1]*n.*ihsqz2
    for mr=1:2*n-1
        pnm[mr+2]=(mr-n).*ihsqz2.*pnm[mr+1]+(2*n-mr+1)*mr*pnm[mr] 
        
    end
    return ((-one(z))^m)*pnm[n+1+m]
end
 


#checking ZZLegendreP3 against Float64(SLegendreP3)
#   for n< 18 the relative error is less than 4e-6
# the error depends on the n amd m values and x values to some extent
#n=1 may have error
"""
    function testZZLegendreP3F()
    compare ZZLegendreP@ with float64 SLegendreP3
"""
function testZZLegendreP3F()
    x=2.
    rerr = 4e-6
    for n=0:17
        for m=0:n  #rand(1:n)
            a=ZZLegendreP3(n,m,x)
            b=Float64(SLegendreP3(n,m,big(x)))
            df=norm(a-b)/norm(b) 
            if ( df > rerr)
            println(n," ", m," ", a," ",b," ",df)
            end
        end
    end
end
testZZLegendreP3F()


