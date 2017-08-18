
# solutions, w, of the differential equation are 
#    associated Legendre functions   
# (1-x*x)d/dx(dw/dx)â¢-2â¢xâ¢(dw/dx)+(n(n+1)-m*m/(1-x*x))w=0
#n is the degree and m is the order.
#for m =0 they are a.k.a. Legendre polynomials
# see https://github.com/pjabardo/Jacobi.jl/blob/master/src/legendre.jl
# for m=0 also.
#
#       SSLegendreP2 and 3 : The finite sum I derived is arranged so that 
# each term is a multiple of the preceeding term. 
#
# LegendreP2 uses the recursion equation 8.5.3 Abramowitz & Stegun (AS) 
# https://www.math.hkbu.edu.hk/support/aands/frameindex.htm
#Handbook of Mathematical Functions  is
#  (n-m+1)P(n+1,m,z)=(2n+1)zP(n,m,z)- (n+m)P(n-1,m,z)
# By normalizing so that that P*P integrated is unity, the new P denoted 
# by PP satisfies an altered 3 term recursion. This is discussed in 
# Numerical Recipes 3rd ed. see https://github.com/milthorpe/
#   SphericalHarmonics.jl/blob/master/src/SphericalHarmonics.jl
# T.Fukushima has Fortran codes (https://static-content.springer.
# com/esm/art%3A10.1007%2Fs00190-011-0519-2/MediaObjects
# /190_2011_519_MOESM_ESM.txt);JamesBremerJr(github) and arXiv:1707.03287
#
#  ZZLegendreP2 and 3  uses a different recursion
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
#P(n,m,z) and Q(n,m,â¡z) (z> 1) are often referred to as the prolate
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
    dblfac=one(x)
    for j=1:2:M    
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




using BenchmarkTools

using PyCall
pyimport_conda("scipy.special","scipy") ####
@pyimport scipy.special as s


function timelpmv()
    for k=1:10^6
        k2=k-1
        x=.1 +k2*.9e-6
         s.lpmv(15,30,x)
    end
end

@btime timelpmv()     

function timeLegendreP2()
    for k=1:10^6
        k2=k-1
        x2= .1 + k2*.9e-6
    LegendreP2(30,15,x2)
    end
end

@btime timeLegendreP2()


"""
    SSLegendreP2(n::Integer,m::Integer,x::Number)
    associated Legendre function 1st kind type 2 , -1.<=x<=1.
    each term in sum is multiple of previous, from SLegendreP2 


"""

function SSLegendreP2(n::Integer,m::Integer,x::Number)#n,m integers 
   
    M = (2*n -1) # M must be odd
    # (2m-1)!!= (2m)!/( m! 2^m)   double factorial 

    DT=one(x) #1.0
    if n > 0
        for j=1:M
            if  isodd(j)
                DT=DT*j
            end
            if  j <= (n-m)
                DT=DT/j
            end
        end  
    end
    MXP=div(m+n,2)
    
    sum= (x^(n+m))*DT/((one(x)-x*x)^(m/2))
    prev=sum
    for p=1:MXP
        term=(-prev)*(n-p+1)*(n+m-2*p+1)*(n+m-2*p+2
        )/(p*(2*n-2*p+1)*(2*n-2*p+2)*x*x)
        sum=sum+term
        prev=term
    end
    return sum
end



"""
    ZZLegendreP2(n::Integer,m::Integer,z::Number)
    associated Legendre function First kind type 2 |x|<=1.
    uses different three term recursion
"""

function ZZLegendreP2(n::Integer,m::Integer,z::Number)    
    nz=1
    pnm = zeros(typeof(z),2*n+1)  
    #fac = prod(2.:n)
    fac= one(z) #1.0
    for k=1:n
        fac=fac*k
    end
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
 

#checking LegendreP2 against scipy.special
#   for n< 18 the relative error is less than 2e-14
#     for n < 45 rel error  less than 2e-13
#     for n < 62 the rel.error is less than 5e-12
#      for n < 75  the rel.error is less than 5e-11
# the error depends on the n amd m values and x values to s(ome extent

using PyCall
pyimport_conda("scipy.special","scipy") ####

@pyimport scipy.special as s

"""
    function testLegendreP2()
    relative error LegendreP2 versus scipy special
"""


function testLegendreP2()
    x=.6
    rerr = 5e-11
    for n=0:75 #17
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
    function testLegendreP2F()
    compare LegendreP2 against Float64(SSLegendreP2)
"""
function testLegendreP2F()
#=
   checking LegendreP2 against Float64SSLegendreP2
   for n < 18 the relative error is less than 7e-15
    error depends on the n amd m values and 
        x values to some extent
 =#
    

    x=.6
    rerr = 7e-15
    for n=1:75
        if n > 18
            rerr=2e-13
        end
        
            
        if n > 44 
            rerr= 5e-12
        end
        if n > 62 
            rerr= 5e-11
        end
        for  m=0:n   #rand(1:n)
            a=LegendreP2(n,m,x)
            b=Float64(SSLegendreP2(n,m,big(x)))
            df=norm(a-b)/norm(b) 
            if ( df > rerr)
            println(n," ", m," ", a," ",b," ",df)
            end
        end
    end        
    
end
testLegendreP2F()


#checking SSLegendreP2 against Float64(SSLegendreP2)
#   for n< 17 the relative error is less than 3e-11
# the error depends on the n amd m values and x values to some extent
"""
    function testSSLegendreP2F()
    compare SSLegendreP2 to float64SSLegendreP2
"""
function testSSLegendreP2F()
    x=.6
    rerr = 3e-11
    for n=0:17
        for m=0:n    #m=rand(0:n)
            a=SSLegendreP2(n,m,x)
            b=Float64(SSLegendreP2(n,m,big(x)))
            df=norm(a-b)/norm(b) 
            if ( df > rerr)
            println(n," ", m," ", a," ",b," ",df)
            end
        end
    end
end
testSSLegendreP2F()



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

#checking SSLegendreP2 against Float64(SSLegendreP2)
#   for n< 17 the relative error is less than 3e-11
# the error depends on the n amd m values and x values to some extent
"""
    function testSSLegendreP2F()
    compare SSLegendreP2 to float64SSLegendreP2
"""
function testSSLegendreP2F()
    x=.6
    rerr = 3e-11
    for n=0:17
        for m=0:n    #m=rand(0:n)
            a=SSLegendreP2(n,m,x)
            b=Float64(SSLegendreP2(n,m,big(x)))
            df=norm(a-b)/norm(b) 
            if ( df > rerr)
            println(n," ", m," ", a," ",b," ",df)
            end
        end
    end
end
testSSLegendreP2F()

# derived from work published at 
# i. a. selezneva, yu. l. ratis, e. hernÃ¡ndez, j. pÃ©rez-quiles 
#and p. fernÃ¡ndez de cÃ³rdoba:
#a code to calculate high order legendre polynomials
#rev. acad. colomb. cienc.: volumen xxxvii, nÃºmero 145 - diciembre 2013
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
    #fac = prod(2.:n)
    fac= one(z) #1.0
    for k=1:n
        fac=fac*k
    end
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
 


   
 """
    SSLegendreP3(n::Integer,m::Integer,x::Number)
    associated Legendre function 1st kind type 3 , x >1.
    each term in sum is multiple of previous. 

"""

function SSLegendreP3(n::Integer,m::Integer,x::Number)#n,m integers 
    
   
    M = (2*n -1) # M must be odd
    # (2m-1)!!= (2m)!/( m! 2^m)   double factorial 
   # dblfac=1
    DT=one(x)   #1.0
    if n > 0
        for j=1:M
            if  isodd(j)
                DT=DT*j
            end    
            if j <= (n-m)
                DT=DT/j
            end
        end
    end
    
    MXP=div(m+n,2)
    sum= (x^(n+m))*DT/((x*x - one(x))^(m/2))
    prev=sum
    for p=1:MXP
        term=(-prev)*(n-p+1)*(n+m-2*p+1)*(n+m-2*p+2
        )/(p*(2*n-2*p+1)*(2*n-2*p+2)*x*x)
        sum=sum+term
        prev=term
    end
    return sum
end


         



 

 


#A representation of associated Legendre function was derived from 8.6.6,
#    8.6.18,8.2.5,3.1.1 in AS (https://www.math.hkbu.edu.hk
#        /support/aands/frameindex.htm). 
#    14.7.8,14.7.10,14.7.11,14.7.14.1.2.2 in DLMF(dlmf.nist.gov),
#similar formula for type 2.
#(8.6.18 AS) Rodrigues' formula for integer n
#  P(n,z)=(1/((2^n)n!)) (d/dz)^n (z^2 - 1)^n
#(8.6.6 AS) P(n,m,z)= ((z^2 -1)^(m/2))(d/dz)^m P(n,z) for z > 1
#   P(n,m,x) = ((-)^m)((1-x^2)^(m/2))(d/dz)^m P(n,x)  for |x| <= 1   
# expand using binomial theorem where a=z^2 and b = -1
#(3.1.1 AS) define binomial coefficient = C(n,k) = n!/((n-k)! k!)
# defining  n!=n*(n-1)*(n-2)*...3*2*1 = factorial(n)
# (a + b)^n = a^n  + C(n,1)*(a^(n-1))*b  + C(n,2)*(a^(n-1))*b^2 
#      +C(n,3)*(a^(n-3))*b^3 +.......+b^n  where n is a positive integer
#differentiate to obtain the finite sum expression



#Benchmarking ZZLegendreP2 against scipy.special

using PyCall
using BenchmarkTools
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

function timeZZ2()
    for k=1:10^6
        k2=k-1
        x2=.1 +k2*.9e-6
        ZZLegendreP2(30,15,x2)
    end
end

function LPMV2()
    for k=1:10^6
        k2=k-1
        x2=.1 +k2*.9e-6
        s.lpmv(15,30,x2)
    end
end

println("lpmv time")
@btime LPMV2()
println("ZZLegendreP2 time")
@btime timeZZ2()

#=
    for k=1:10^6
        k2=k-1
        x2=.1 +k2*.9e-6
        s.lpmv(15,30,x2)
    end


println("lpmv time")
@btime LPMV2()
println("ZZLegendreP2 time")
@btime timeZZ2()

=#

        



        



#checking ZZLegendreP2 against scipy special
#   for n< 41 the relative error is less than 2.4e-9
#    for n < 61 the relative error is less than 6e-5
# the error depends on the n amd m values and x values to some extent
#
using PyCall
@pyimport scipy.special as s

"""
    function testZZLegendreP2(
    compare ZZLegendreP2 with scipy special
"""
function testZZLegendreP2()
    x=.6
    rerr = 2.4e-9# 5e-7
    #for n =51: 80
    for n=1:40
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


#checking ZZLegendreP2 against Float64(SSLegendreP2)
#   for n< 41 the relative error is less than 2.4e-9
#    for n <61 the relative error is less than 6e-5
# the error depends on the n amd m values and x values to some extent
#  largest errors for n=m or n= m+1 (could be fixed?)
"""
    function testZZLegendreP2F()
    compare ZZLegendreP2 with Float64 SSLegendreP2
"""
function testZZLegendreP2F()
    x=.6
    rerr = 2.4e-9 
    for n=1:40
        for m= 0:n  #rand(1:n)
            a=ZZLegendreP2(n,m,x)
            b=Float64(SSLegendreP2(n,m,big(x)))
            df=norm(a-b)/norm(b) 
            if ( df > rerr)
            println(n," ", m," ", a," ",b," ",df)
            end
        end
    end
end
testZZLegendreP2F()


# derived from work published at 
# i. a. selezneva, yu. l. ratis, e. hernÃ¡ndez, j. pÃ©rez-quiles 
#and p. fernÃ¡ndez de cÃ³rdoba:
#a code to calculate high order legendre polynomials
#rev. acad. colomb. cienc.: volumen xxxvii, nÃºmero 145 - diciembre 2013
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
    #fac = prod(2.:n)
    fac= one(z) #1.0
    for k=1:n
        fac=fac*k
    end
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
 




#checking ZZLegendreP3 against Float64(SSLegendreP3)
#   for n< 17 the relative error is less than 6.4e-7 for z=2.
#    for n<24 the relative error is less than .0054
# for z=3. n<17 the relative eror is less than 2e-7
# the error depends on the n amd m values and x values to some extent
#n=1 may have error
"""
    function testZZLegendreP3F()
    compare ZZLegendreP3 with float64 SSLegendreP3
"""
function testZZLegendreP3F()
    x=3.#2.
    rerr = 2e-7 #6.4e-7
    for n=1:16
        for m=0:n  #rand(1:n)
            a=ZZLegendreP3(n,m,x)
            b=Float64(SSLegendreP3(n,m,big(x)))
            df=norm(a-b)/norm(b) 
            if ( df > rerr)
            println(n," ", m," ", a," ",b," ",df)
            end
        end
    end
end
testZZLegendreP3F()






using PyCall
@pyimport scipy.special as s

"""
    function ZZSS2(n,m,x,ESTOP)
    Values near .99999 < x <= 1.000000 and -1.000000 <x<-.99999
    causing problems
    Identify n,m, x ,when zz is inaccurate

"""
function ZZSS2(n,m,x,ESTOP)
    valueZZ=Float64(ZZLegendreP2(n,m,big(x)))
    #valueS=Float64(SLegendreP2(n,m,big(x)))
    valueSS=Float64(SSLegendreP2(n,m,big(x)))
    valuezz=ZZLegendreP2(n,m,x)
    #valuess=SSLegendreP2(n,m,x)
    valuepy=s.lpmv(m,n,x)
    prerr=abs((valuezz-valuepy)/valuepy) 
    zzrerr=abs((valuezz-valueZZ)/valueZZ)
    if zzrerr > ESTOP
        println()
        println("n=$n, m=$m , x=$x , zz=$valuezz")
        println("SS=$valueSS,ZZ=$valueZZ, py=$valuepy ")
        println(" zzrerr= $zzrerr, prerr= $prerr")
    end 
    
end


ESTOP=.000000000001

n=10
for ix=1:10
    x=1. - (.1)^ix
        for m=0:n
            ZZSS2(n,m,x,ESTOP)
        end
end     
 
n=10
for ix=1:10
    x=-1. + (.1)^ix
        for m=0:n
            ZZSS2(n,m,x,ESTOP)
        end
end










 

"""
    function SSZZ1(n,m,x,ESTOP)
    Values near 1.< x < 1.00001 causing problems
    Identify n,m, x ,when zz goes negative, its underflow
    The cutoff does change with n values
    The user has the option to accept reduced accuracy
    for certain values of n,m,x or to increase the accuracy
    by using higher precision 
"""
function SSZZ1(n,m,x,ESTOP)
    valueZZ=Float64(ZZLegendreP3(n,m,big(x)))#high precision
    #valueS=Float64(SLegendreP3(n,m,big(x)))
    valueSS=Float64(SSLegendreP3(n,m,big(x)))
    valuezz=ZZLegendreP3(n,m,x)
    #valuess=SSLegendreP3(n,m,x)
    zzrerr=abs((valuezz-valueZZ)/valueZZ)
    #ssrerr=abs((valuess-valueSS)/valueSS)
    
    if zzrerr > ESTOP
    #if sign(valuezz) < 1 
        println()
        println("n= $n, m= $m, x= $x")
        println("SS=$valueSS,ZZ=$valueZZ, zz= $valuezz, zzrerr= $zzrerr")
        
        
        #println("S=$valueS,SS=$valueSS,ZZ=$valueZZ, n= $n, m= $m, x= $x")
        #println("  ssrerr= $ssrerr,ss= $valuess,zz= $valuezz, zzrerr= $zzrerr")
        
    end
end 





ESTOP = .000000000001 #e-12

n=10

for ix =1:7   
    for m=0:n
        x=1.+.1^(ix)
       
        SSZZ1(n,m,x,ESTOP)
    end
end

println("  ")

 n=80
ESTOP=.000000000001 # e-12
for ix=1:7
    for m=0:n 
        x= 1.+  (.1)^(ix)   
    
        SSZZ1(n,m,x,ESTOP)
    end
end





"""
    function TableSSZZ(n,m,x,T) 
    compare  SSLegendreP3,ZZLegendreP3 and Table1 value is T
    with and without arbitrary precision to 
    Table 1 from Computer Physics Comm.  181 (2010)2091-7 
    User has the option to accept reduced accuracy for certain
    values of n,m,x or to increase the accuracy by using
    higher values of precision in arbitrary precision arithmetic
"""
function TableSSZZ(n,m,x,T)
    valueZZ=Float64(ZZLegendreP3(n,m,big(x))) #higher precision
    #valueS=Float64(SLegendreP3(n,m,big(x)))
    valueSS=Float64(SSLegendreP3(n,m,big(x)))
    valuezz=ZZLegendreP3(n,m,x)
    #valuess=SSLegendreP3(n,m,x)
    zzrerr=abs((valuezz-T)/T)
    ZZrerr=abs((valueZZ-T)/T)
    println()
    println("n=$n, m=$m, x=$x")
    println("SS=$valueSS,ZZ=$valueZZ, T=$T ,zz=$valuezz " )
   # println("SS=$valueSS,ZZ=$valueZZ, T=$T , n= $n, m= $m, x= $x")
   # println(" ZZrerr=$ZZrerr, ss=$valuess,zz=$valuezz, zzrerr=$zzrerr")
     println(" ZZrerr=$ZZrerr, zzrerr=$zzrerr")
  
end
x=5.0 
m=0
n=0
T=1.0
x=5.0 
m=0
n=0
T=1.0
TableSSZZ(n,m,x,T)
n=1
T=5.0
TableSSZZ(n,m,x,T)
n=2
T=37.0
TableSSZZ(n,m,x,T)
n=3
T=305.
TableSSZZ(n,m,x,T)
n=4
T=2641.
TableSSZZ(n,m,x,T)
m=1
n=1
T=4.898979485566
TableSSZZ(n,m,x,T)
n=2
T=73.4846922835
TableSSZZ(n,m,x,T)
n=3
T=9.112101843153e2
TableSSZZ(n,m,x,T)
n=4
T=1.053280589397e4
TableSSZZ(n,m,x,T)
m=2
n=2
T=72.
TableSSZZ(n,m,x,T)
n=3
T=1800.
TableSSZZ(n,m,x,T)
n=4
T=3.132e4
TableSSZZ(n,m,x,T)
m=3
n=3
T=1.763632614804e3
TableSSZZ(n,m,x,T)
n=4
T=6.172714151814e4
TableSSZZ(n,m,x,T)
m=4
n=4
T=6.048e4
TableSSZZ(n,m,x,T)
x=2.0
m=0
n=0
T=1.0
TableSSZZ(n,m,x,T)
n=1
T=2.0
TableSSZZ(n,m,x,T)
n=2
T=5.5
TableSSZZ(n,m,x,T)
n=3
T=17.
TableSSZZ(n,m,x,T)
n=4
T=55.375
TableSSZZ(n,m,x,T)
m=1
n=1
T=1.732050807569
TableSSZZ(n,m,x,T)
n=2
T=1.039230484541e1
TableSSZZ(n,m,x,T)
n=3
T=4.936344801571e1
TableSSZZ(n,m,x,T)
n=4
T=2.165063509461e2
TableSSZZ(n,m,x,T)
m=2
n=2
T=9.0
TableSSZZ(n,m,x,T)
n=3
T=90.
TableSSZZ(n,m,x,T)
n=4
T=607.5
TableSSZZ(n,m,x,T)
m=3
n=3
T=7.794228634060e1
TableSSZZ(n,m,x,T)
n=4
T=1.091192008768e3
TableSSZZ(n,m,x,T)
m=4
n=4
T=945.0
TableSSZZ(n,m,x,T)
x=1.000001
m=0
n=0
T=1.0
TableSSZZ(n,m,x,T)
n=1
T=1.0000010
TableSSZZ(n,m,x,T)
n=2
T=1.000003000001
TableSSZZ(n,m,x,T)
n=3
T=1.000006000007
TableSSZZ(n,m,x,T)
n=4
T=1.000010000022
TableSSZZ(n,m,x,T)
m=1
n=1
T=1.414213915900e-3
TableSSZZ(n,m,x,T)
n=2
T=4.242645990341e-3
TableSSZZ(n,m,x,T)
n=3
T=8.485304708618e-3
TableSSZZ(n,m,x,T)
n=4
T=1.414220279870e-2
TableSSZZ(n,m,x,T)
m=2
n=2
T=6.000002999773e-6
TableSSZZ(n,m,x,T)
n=3
T=3.000004499888e-5
TableSSZZ(n,m,x,T)
n=4
T=9.000025499681e-5
TableSSZZ(n,m,x,T)
m=3
n=3
T=4.242643868860e-8
TableSSZZ(n,m,x,T)
n=4  #underflow
T=2.969853678052e-7
TableSSZZ(n,m,x,T)
m=4 #underflow
n=4
T=4.200004199683e-10
TableSSZZ(n,m,x,T)








