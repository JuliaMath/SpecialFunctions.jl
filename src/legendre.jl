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

    ABC_recurrence(n, x,
        m->(2m-1)//m, 0, m->(m-1)//m,       # A_m, B_m, C_m
        1, 0, 1)                            # P_{0,0}, P_{1,0}, P_{1,1}
end

@doc raw"""
    ABC_recurrence(n, x, A, B, C, c0_0, c1_0, c1_1)

Evaluate polynomials by the recurrence
```math
Q_n(x)
=
\begin{cases}
    c^{(0)}_0
    & \text{for} \quad n = 0
    \\
    c^{(1)}_0 + c^{(0)}_1 x
    & \text{for} \quad n = 1
    \\
    (A_n x + B_n) p_{n-1}(x) - C_n p_{n-2}(x)
    & \text{for} \quad n \geq 2
    \,.
\end{cases}
```

This will be used to compute Legendre ``P_n``, Laguerre ``L_n``, Hermite ``H_n``,
and Chebyshev polynomials of 1st kind ``T_n`` and 2nd kind ``U_n``.

# Arguments
- `n::Integer`:     Polynomial order
- `x::Real`:        evaluation point
- `A::Union{Integer,Function}`:    recurrence coefficient ``A_m``
- `B::Union{Integer,Function}`:    recurrence coefficient ``B_m``
- `C::Union{Integer,Function}`:    recurrence coefficient ``C_m``
- `c0_0::Real`:     polynomial coefficient ``c^{(0)}_0``
- `c1_0::Real`:     polynomial coefficient ``c^{(1)}_0``
- `c1_1::Real`:     polynomial coefficient ``c^{(1)}_1``

#Implementation
For some polynomials the coefficients ``A_m, B_m, C_m`` are independent of ``m`` and for
others they aren't. Which is why we use `Union{Integer,Function}`.
"""
function ABC_recurrence(n::Integer, x::Real,
        A::Union{Integer,Function}, B::Union{Integer,Function}, C::Union{Integer,Function},
        c0_0::Integer, c1_0::Integer, c1_1::Integer)

    p_prev_prev = c0_0
    p_prev      = c1_0 + c1_1 * x

    if      n == 0
        return p_prev_prev
    elseif  n == 1
        return p_prev
    end

    p = missing
    for m = 2:n
        p = ABC_recurrence_step(m, x, A, B, C, p_prev, p_prev_prev)

        p_prev_prev  = p_prev
        p_prev       = p
    end

    p
end

function ABC_recurrence_step(m::Integer, x::Real,
    A::Function, B::Integer, C::Function,
    p_prev::Real, p_prev_prev::Real)

    (A(m)*x + B) * p_prev - C(m) * p_prev_prev
end
