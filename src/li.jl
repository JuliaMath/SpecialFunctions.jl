"""
    Li(s, z)

 Calculates the Polylogarithm function ``Li_s(z)`` defined by

```math
    L_s = \sum_{n=1}^{\infty} \frac{z^n}{n^s} 
```

 For ideas going into this see
       
 + Crandall, "Note on fast polylogarithm computation", 2006, 
           which focusses on the case where s=n (integer and real)    
           http://www.wolfgang-ehrhardt.de/Polylog.pdf

 + Vepstas, "AN EFFICIENT ALGORITHM FOR ACCELERATING THE CONVERGENCE
                 OF OSCILLATORY SERIES, USEFUL FOR COMPUTING THE
                 POLYLOGARITHM AND HURWITZ ZETA FUNCTIONS", 2007
           which treats the general case, but presumes arbitrary precision arithmetic
           https://arxiv.org/abs/math/0702243
             
 + Wood, "The computation of Plylogarithms", 1992
           which focusses on s=n, integer and real, and which has formatting issues making
           it hard to read correctly.
           https://www.cs.kent.ac.uk/pubs/1992/110/

 + Maximon, "The dilogarithm function for complex argument", 2003
           which provides useful test cases for s=2

 + Zagier,  "The dilogarithm function in geometry and number theory", 1989 
           similar to Maximon
           
  Of these the only one that actually specifies a full algorithm is
  Crandall, and it also treats special cases more carefully, so this
  is the one that I have paid most attention to. However, extending it
  for s on the whole complex plane requires some additions, and many
  of these are actually most nicely documented on the wikipedia page
    
  + https://en.wikipedia.org/wiki/Polylogarithm

  With further details at

  + http://mathworld.wolfram.com/Polylogarithm.html
  + http://dlmf.nist.gov/25.12#ii
  + http://mathworld.wolfram.com/Trilogarithm.html
  + http://functions.wolfram.com/ZetaFunctionsandPolylogarithms/PolyLog/

  The wiki page points out some errors in earlier works, but not all
  parts on the page have references, and not all statements seem to
  come from any of the listed references?

  The code draws heavily on existing functions, in particular the
  Hurwitz-zeta function, which is aliased to zeta(s,q) in Julia.

  Accuracy has been tested using many of the identities known for Li
  and relations to known functions as special cases, and by comparison
  to `polylog(s, z)` in the `mpmath` arbitrary-precision package in Python. 

      http://mpmath.org/doc/current/functions/zeta.html

  The latter shows deviations of the order of 

  + 10^{Im(s) - 20} for Im(s) < 0
  + 10^{Im(s) - 20} for Im(s) > 0
                                                  
  It isn't clear whether we can do better than this with
  double-precision arithmetic.

## Arguments
* `s::Complex`: the 'fractional' coefficient
* `z::Complex`: the point at which to calculate it
* `accuracy::Real=1.0e-18`: nominal accuracy of calculation, but mainly useful for testing

## Examples
```jldoctest
julia> Li(-1.0, 0.0) 
(0.0,1)
```
"""
function Li(s::Number, z::Number, accuracy::Real=1.0e-18)
    T = 0.5 # threshold at which we change algorithms
    if z ≈ 1.0
        if real(s) > 1
            return zeta(s)
        else
            return Inf
        end
    elseif z ≈ -1.0
        return -eta(s)
    elseif s ≈ 0.0
        return z ./ (1-z)
    elseif abs(z) <= T
        return Li_direct(s, z, accuracy)
    elseif abs(z) >= 1/T && isinteger(s) && real(s) < 0
        # use inversion formula to calculate in terms of Li_n(1/z)
        # but note for negative integer s, it collapses to something small
        return -(-1.0)^s .*Li_direct(s, 1/z, accuracy)
    elseif  abs(z) >= 1/T
        # use inversion formula to calculate in terms of Li_s(1/z)
        G = (2*pi*im).^s .* zeta( 1-s, 0.5 - log(-1/complex(z))/(2*pi*im) ) ./  gamma(s)
        F = complex(-1.0)^s .* Li_direct(s, 1/z, accuracy)
        A = 2*pi*im*log(complex(z))^(s-1)/(gamma(s))
        if ( isreal(z) && real(z)>=1 )
            Θ = 1
        else
            Θ = 0
        end
        # println("G = $G, F=$F, Θ=$Θ, A=$A")
        return ( G - F - Θ*A )
    else 
        # power series around mu=0, for z = e^mu
        Li_series_mu(s, z, accuracy)
    end
end
     

"""
    Polylog(s, z)

    This is just an alias of Li(s, z)

"""
polylog(s::Number, z::Number ) = Li(s, z)

    
####################################################################
#### these are component functions and aren't exported at this point
     
# Dirichlet beta function, for testing results
#   https://en.wikipedia.org/wiki/Dirichlet_beta_function
function Dbeta(s::Number)
    β = 4.0^(-s) * ( zeta(s,0.25) - zeta(s,0.75) )
end
    
function Li_direct(s::Number, z::Number, accuracy=1.0e-18)
    # calculate using direct definition
    if abs(z) > 1 || ( abs(z) ≈ 1  && real(s) <= 2)
        warn("Should have |z| < 1 or (|z|==1 and Re(s)>2)")
        throw(DomainError)
    end
    if abs(z) > 1/2
        warn("Slow convergence for  |z| > 1/2")
    end
    total = 0.0
    L = ceil(-log10(accuracy)*log2(10)) # summation limit from Crandall, which is conservative
    for n=1:L
        a = z^n / n^s
        total += a
        # println("   total = $total")
    end
    return total
end

function Li_series_mu(s::Number, z::Number, accuracy=1.0e-18)
    # calculate using power series around μ = log(z) = 0
    μ = log(complex(z))
    # println("μ = $μ") 
    if abs(μ) > 2*pi
        warn("Should have |μ| < 2*pi, but |μ| = $(abs(μ))")
        throw(DomainError)
    end
    L = Int(ceil(-log10(accuracy)*log2(10))) # summation limit from Crandall, which is conservative
    if isinteger(s)
        n = Int(round(s))
        if n>0
            # Crandall's 1.4 for s integer
            total = μ^(n-1)*(harmonic(n-1) - log(-μ))/gamma(n)
            # println("   μ=$μ, total = $total")
            for m=0:L
                if n - m != 1
                    total += μ^m * zeta(n - m) / gamma(m+1)
                end
                # println("   m=$m, total = $total")
            end
            # println("   μ=$μ, total = $total")
            A = 2*pi*im*log(complex(z))^(s-1)/gamma(n)
            if  isreal(z) && real(z)>=1 
                total -= A
            end
            # println("   μ=$μ, total = $total")
        elseif n==0
            total = z ./ (1-z)
        elseif n<0
            # Crandall's 1.5 for s integer 
            total = factorial(-n) * (-μ)^(n-1)
            for k=0:L
                total -= μ^k * bernoulli(k-n+1, 0.0) / ( gamma(k+1)*(k-n+1) )
            end
        else
            error("Should not get this case")
        end
    else
        # equivalent of Crandalls 1.4 for s non-integer 
        total = gamma(1-s) * (-μ)^(s-1)
        for k=0:L
            total += μ^k * zeta(s-k)/factorial(Float64(k))
        end
        A = 2*pi*im*log(complex(z))^(s-1)/(gamma(s))
        if isreal(z) && real(z)>=1 
            total -= A
        end
    end
    return total
end
