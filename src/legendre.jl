@doc raw"""
    legendreP(n, x)

Evaluate the Legendre Polynomial ``P_n(x)`` of order ``n`` at position ``x``, defined by the
Bonnet recursion

```math
P_n(x)
=
\begin{cases}
    1
    & \text{for} \quad n = 0
    \\
    x
    & \text{for} \quad n = 1
    \\
    \frac{2n-1}{n} P_{n-1}(x) - \frac{n-1}{n} P_{n-2}(x)
    & \text{for} \quad n \geq 2
    \,.
\end{cases}
```

External links: [DLMF](https://dlmf.nist.gov/14.7.E1),
[Wikipedia](https://en.wikipedia.org/wiki/Legendre_polynomials).
"""
function legendreP(n::Integer, x::Real)
    if n < 0
        throw(DomainError(n, "must be non-negative"))
    end
    if x < -1 || x > 1
        throw(DomainError(x, "must be in the range [-1,1]"))
    end

    legendreP_bonnet(n, x)
end

function legendreP_bonnet(n::Integer, x::Real)

    legP_prev_prev  = 1
    legP_prev       = x

    if      n == 0
        return legP_prev_prev
    elseif  n == 1
        return legP_prev
    end

    legP = missing
    for m = 2:n
        legP = ((2*m-1)//m) * x * legP_prev - ((m-1)/m) * legP_prev_prev

        legP_prev_prev  = legP_prev
        legP_prev       = legP
    end

    legP
end
