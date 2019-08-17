@doc raw"""
    chebyshevT(n, x)

Evaluate the Chebyshev polynomials of the first kind ``T_n(x)`` of order ``n`` at position
``x``, defined by

```math
T_n(x)
=
\begin{cases}
    1
    & \text{for} \quad n = 0\,, \; x \in [-1,1]
    \\
    x
    & \text{for} \quad n = 1\,, \; x \in [-1,1]
    \\
    2 x T_{n-1}(x) - T_{n-2}(x)
    & \text{for} \quad n \geq 2\,, \; x \in [-1,1]
    \,.
\end{cases}
```

External links: [DLMF](https://dlmf.nist.gov/18.3.T1),
[Wikipedia](https://en.wikipedia.org/wiki/Chebyshev_polynomials).
"""
function chebyshevT(n::Integer, x::Real)
    if n < 0
        throw(DomainError(n, "must be non-negative"))
    end
    if x < -1 || x > 1
        throw(DomainError(x, "must be in the range [-1,1]"))
    end

    ABC_recurrence(n, x,
        2, 0, 1,                # A_m, B_m, C_m
        1, 0, 1)                # P_{0,0}, P_{1,0}, P_{1,1}
end

@doc raw"""
    hermiteH(n, x)

Evaluate the (physicists') Hermite polynomials ``H_n(x)`` of order ``n`` at position ``x``,
defined by

```math
H_n(x)
= (-1)^n e^{x^2} \frac{\mathrm{d}^n}{\mathrm{d}x^n} e^{-x^2}
\quad \text{for} \quad
x \in \mathbb{R}, \; n = 0, 1, 2, \dots
```

External links: [DLMF](https://dlmf.nist.gov/18.3.T1),
[Wikipedia](https://en.wikipedia.org/wiki/Hermite_polynomials).
"""
function hermiteH(n::Integer, x::Real)
    if n < 0
        throw(DomainError(n, "must be non-negative"))
    end

    ABC_recurrence(n, x,
        2, 0, m->2(m-1),        # A_m, B_m, C_m
        1, 0, 2)                # P_{0,0}, P_{1,0}, P_{1,1}
end

@doc raw"""
    laguerreL(n, x)

Evaluate the Laguerre polynomial ``L_n(x)`` of order ``n`` at position ``x``, defined by

```math
L_n(x)
=\frac{e^x}{n!} \frac{\mathrm{d}^n}{\mathrm{d}x^n} \left(e^{-x} x^n\right)
\quad \text{for} \quad
x \geq 0, \; n = 0, 1, 2, \dots
```

External links: [DLMF](https://dlmf.nist.gov/18.3.T1),
[Wikipedia](https://en.wikipedia.org/wiki/Laguerre_polynomials).
"""
function laguerreL(n::Integer, x::Real)
    if n < 0
        throw(DomainError(n, "must be non-negative"))
    end
    if x < 0
        throw(DomainError(x, "must be nonnegative"))
    end

    ABC_recurrence(n, x,
        m->-1//m, m->2-1//m, m->1-1//m,     # A_m, B_m, C_m
        1, 1, -1)                           # P_{0,0}, P_{1,0}, P_{1,1}
end

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
others they aren't. Which is why we use `Union{Integer,Function}` and functions can be
easily created by anonymous functions.

# Reference
> JIN, J. M.;
> JJIE, Zhang Shan.
> "Computation of special functions".
> Wiley, 1996.
> (p. 23)
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

# differentiate for combinations of coefficients: Real vs Function
function ABC_recurrence_step(m::Integer, x::Real,
    A::Function, B::Function, C::Function,
    p_prev::Real, p_prev_prev::Real)

    (A(m)*x + B(m)) * p_prev - C(m) * p_prev_prev
end
function ABC_recurrence_step(m::Integer, x::Real,
    A::Function, B::Integer, C::Function,
    p_prev::Real, p_prev_prev::Real)

    (A(m)*x + B) * p_prev - C(m) * p_prev_prev
end
function ABC_recurrence_step(m::Integer, x::Real,
    A::Integer, B::Integer, C::Function,
    p_prev::Real, p_prev_prev::Real)

    (A*x + B) * p_prev - C(m) * p_prev_prev
end
