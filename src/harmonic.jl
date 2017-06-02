"""
    harmonic(n)

 Calculates harmonic numbers
   e.g., see http://mathworld.wolfram.com/HarmonicNumber.html

## Arguments
* `n::Integer`: index of the Harmonic number to calculate

## Examples
```jldoctest
julia> harmonic(2)
1.5
```
"""
function harmonic(n::Integer)
    # γ = euler_mascheroni_const = 0.577215664901532860606512090082402431042 # http://oeis.org/A001620
    if n <= 0
        throw(DomainError())
    end
    if n <= 10
        # get exact values for small n
        total = 0.0
        for k=1:n
            total +=  1.0 / k
        end
        return total
    else
        # numerical approximation for larger n
        return γ + digamma(n+1)
    end
end


"""
    harmonic(n,r)

 Calculates generalized harmonic numbers
   e.g., see http://mathworld.wolfram.com/HarmonicNumber.html

## Arguments
* `n::Integer`: index 1 of the Harmonic number to calculate
* `r::Real`: index 2 of the Harmonic number to calculate

It should be possible to extend this to complex r, but that requires more testing.

## Examples
```jldoctest
julia> harmonic(2,1)
1.5
```
"""
function harmonic(n::Integer, r::Real)
    total = 0.0
    if n <= 0
        throw(DomainError())
    end
    for k=1:n
        total +=  1.0 / k^r
    end
    return total
end
