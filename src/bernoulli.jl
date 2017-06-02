#
# Bernoulli numbers (as rationals up to B_{34}) and Bernoulli polynomials

"""
    bernoulli(n)

 Calculates the first 34 Bernoulli numbers B_n  (of the first-kind or NIST type) 
 e.g., see

 + http://mathworld.wolfram.com/BernoulliNumber.html
 + https://en.wikipedia.org/wiki/Bernoulli_number
 + http://dlmf.nist.gov/24

 N.B. Bernoulli numbers of second kind only seem to differ in that B_1 = + 1/2 (instead of -1/2)

## Arguments
* `n::Integer`: the index into the series, n=0,1,2,3,...,34 (for larger n use ``bernoulli(n,0.0)`` )

 We only provide the 1st 34 as beyond this, we can't return Int64 rationals, so best to compute
 the real approximation using ``bernoulli(n,0.0)``

## Examples
```jldoctest
julia> bernoulli(6)
1 // 42
```
"""
function bernoulli(n::Int)
    # this just does a lookup -- seemed like it would be easier to code and faster
    # for the size of numbers I am working with
    if n<0
        throw(DomainError)
    end
    if n>34
        # warn("If n > 34, then the numerator needs Int128 at least, and worse, so this code is not the code you want. Try using bernoulli(n, 0.0) to get a floating point approximation to the result.")
        throw(DomainError)
    end

    # Denominator of Bernoulli number B_n
    #   http://oeis.org/A027642
    D = [2, 6, 1, 30, 1, 42, 1, 30, 1, 66, 1, 2730, 1, 6, 1, 510, 1, 798, 1, 330, 1, 138, 1, 2730, 1, 6, 1, 870, 1, 14322, 1, 510, 1, 6, 1, 1919190, 1, 6, 1, 13530, 1, 1806, 1, 690, 1, 282, 1, 46410, 1, 66, 1, 1590, 1, 798, 1, 870, 1, 354, 1, 56786730]

    # Numerator of Bernoulli number B_n (storing 62 of these because they are easy)
    #   http://oeis.org/A027641
    N = [-1, 1, 0, -1, 0, 1, 0, -1, 0, 5, 0, -691, 0, 7, 0, -3617, 0, 43867, 0, -174611, 0, 854513, 0, -236364091, 0, 8553103, 0, -23749461029, 0, 8615841276005, 0, -7709321041217, 0, 2577687858367, 1]
    
    if n==0
        return 1 
    else
        return N[n] // D[n]
    end
end

# get the Bernoulli polynomials from the Hurwitz-zeta function (which is already defined)
#  
# 
"""
    bernoulli(n, x)

 created: 	Tue May 23 2017 
 email:   	matthew.roughan@adelaide.edu.au
 (c) M Roughan, 2017

 Calculates Bernoulli polynomials from the Hurwitz-zeta function using

```ζ(-n,x) = -B_{n+1}(x)/(n+1), for Re(x)>0,```

 which is faster than direct calculation of the polynomial except for n<=4.

 e.g., see
 
 + https://en.wikipedia.org/wiki/Bernoulli_polynomials
 + http://dlmf.nist.gov/24

## Arguments
* `n::Integer`: the index into the series, n=0,1,2,3,...
* `x::Real`: the point at which to calculate the polynomial

## Examples
```jldoctest
julia> bernoulli(6, 1.2)
0.008833523809524069
```
"""
function bernoulli(n::Int, x::Real)
    if n<0
        throw(DomainError)
    end
    if n == 0
        return 1 # zeta formula doesn't hold for n=0, so return explicit value
    elseif n == 1
        return x-0.5 # get some easy cases out of the way quickly
    elseif n == 2
        return x^2 - x + 1.0/6.0
    elseif n == 3
        return x^3 - 1.5*x^2 + 0.5*x
    elseif n == 4
        return x^4 - 2.0*x^3 +     x^2 - 1/30.0
    # elseif n == 5
    #     return x^5 - 2.5*x^4 +(5.0/3.0)*x^3 - x/6.0 # this isn't faster than just doing the zeta directly
    end

    # direct summation is slower than the zeta function below, even for small n
    # if n <= 34
    #     # direct summation for reasonably small values of coefficients
    #     total = 0.0
    #     for k=0:n
    #         total +=  binomial.(n,k) .* bernoulli.(k) .* x.^(n-k)
    #     end
    #     return total
    # else

    if x > 0
        return -n*zeta(1-n, x)
    else
        # comments in the gamma.jl note that identity with zeta(s,z) only works for Re(z)>0
        # so exploit symmetries in B_n(x) to compute recursively for x<=0
        return bernoulli(n, x+1) - n*x^(n-1)
    end
end
