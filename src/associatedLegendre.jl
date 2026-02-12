
# solutions, w, of the differential equation are 
#    associated Legendre functions   

# (1-x*x)d/dx(dw/dx)⁢-2⁢x⁢(dw/dx)+(n(n+1)-m*m/(1-x*x))w=0

    #for m =0 they are a.k.a. Legendre polynomials
#n is the degree and m is the order 
#The recursion equation 8.5.3 Abramowitz & Stegun (AS) Handbook of
#   Mathematical Functions  used is
#  (n-m+1)P(n+1,m,z)=(2n+1)zP(n,m,z)- (n+m)P(n-1,m,z)
# (for a recent comparison of transforms and recursions
#   see  https:doi.org/10.1515/jag-2016-0032)
#this recursion is valid for legendre functions of the first kind 
#(P(n,m,z))and second kind Q(n,m,z) , 
#   (for types 1,2,3 ,definition from Mathematica documentation
#   type 1 (?) is m=0, type 2 is -1<=x<=1 and type 3 is z> 1
#   http://mpmath.org/doc/0.18/functions/orthogonal.html
#Thus we define functions  LegendreP2,....P3,....Q2,...Q3  
# for n and m positive integers 0<=m<=n
# refer to http://dlmf.nist.gov/14.6
#P(n,m,z) and Q(n,m,⁡z) (z> 1) are often referred to as the prolate
#spheroidal harmonics of the first and second kinds, respectively
#  We can define spherical harmonics (there are other definitions)
# Y(n,m,theta,phi)= LegendreP2(n,|m|,theta)*exp(im*m*phi)*((-)^(|m|-m)/2))
# *sqrt((2n+1)*factorial(n-|m|)/(4*pi*factorial(n+|m|))) 
# refer to https://en.wikipedia.org/wiki/Spherical_harmonics



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






function LegendreP3(n::Integer,m::Integer,z::Number) 
  
    #reliable for values of n <= 17 for z=2. 
    # 0 <= m  <= n
    # z > 1.
    M = (2*m -1) # M must be odd
    # (2m-1)!!= (2m)!/( m! 2^m)   double factorial 
          dblfac=1
          for j=1:M
            if  iseven(j)
            continue
           end
          dblfac=j*dblfac  
          end
    
    
    pj2=dblfac*(z*z - one(z))^(m/2)
        pj1=z*(2.*m+1.)*pj2  
    
    if n==m 
        return (pj2) 
        elseif n== m+1 
        return (pj1)
     end
    
   
    for j = m+2 :n 
     pjj=(z*(2.*j-1.)*pj1 - pj2*(j +m-1.)) /(j-m)
                pj2=pj1
                pj1=pjj
    end
return pj1 
end





