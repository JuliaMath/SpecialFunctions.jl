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
        2, 0, 1,                            # A_m, B_m, C_m
        y->1,                               # P_0(y)
        y->y)                               # P_1(y)
end

@doc raw"""
    chebyshevU(n, x)

Evaluate the Chebyshev polynomials of the second kind ``U_n(x)`` of order ``n`` at position
``x``, defined by

```math
U_n(x)
=
\begin{cases}
    1
    & \text{for} \quad n = 0\,, \; x \in [-1,1]
    \\
    2x
    & \text{for} \quad n = 1\,, \; x \in [-1,1]
    \\
    2x U_{n-1}(x) - U_{n-2}(x)
    & \text{for} \quad n \geq 2\,, \; x \in [-1,1]
    \,.
\end{cases}
```

External links: [DLMF](https://dlmf.nist.gov/18.3.T1),
[Wikipedia](https://en.wikipedia.org/wiki/Chebyshev_polynomials).
"""
function chebyshevU(n::Integer, x::Real)
    if n < 0
        throw(DomainError(n, "must be non-negative"))
    end
    if x < -1 || x > 1
        throw(DomainError(x, "must be in the range [-1,1]"))
    end

    ABC_recurrence(n, x,
        2, 0, 1,                            # A_m, B_m, C_m
        y->1,                               # P_0(y)
        y->2y)                              # P_1(y)
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
        2, 0, m->2(m-1),                    # A_m, B_m, C_m
        y->1,                               # P_0(y)
        y->2y)                              # P_1(y)
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
        y->1,                               # P_0(y)
        y->1-y)                             # P_1(y)
end

@doc raw"""
    legendreP(n[, m], x)

Evaluate the Legendre polynomial ``P_n(x)`` of degree ``n`` at position ``x``.
If the order ``m`` is supplied then the associated Legendre polynomials ``P_n^{(m)}(x)`` is
computed.
Their definitions are given by
```math
\begin{align}
    P_n^{(m)}(x)
    &=
    (-1)^m (1-x^2)^\frac{m}{2} \frac{\mathrm{d}^m}{\mathrm{d}x^m} \left( P_n(x) \right)
    \,, \quad m = 0, 1, \dots, n \,, \; x \in [-1,1]
    \\
    P_n(x)
    &=
    \begin{cases}
        1
        & \text{for} \quad n = 0
        \\
        x
        & \text{for} \quad n = 1
        \\
        \frac{2n-1}{n} x P_{n-1}(x) - \frac{n-1}{n} P_{n-2}(x)
        & \text{for} \quad n \geq 2
        \,.
    \end{cases}
\end{align}
```

External links: [DLMF - Legendre polynomial](https://dlmf.nist.gov/14.7.E1),
[DLMF - associated Legendre polynomial](https://dlmf.nist.gov/14.7.E8),
[Wikipedia - Legendre polynomial](https://en.wikipedia.org/wiki/Legendre_polynomials),
[Wikipedia - associated Legendre polynomial](https://en.wikipedia.org/wiki/Associated_Legendre_polynomials).
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
        y->1,                               # P_0(y)
        y->y)                               # P_1(y)
end
function legendreP(n::Integer, m::Integer, x::Real)
    if n < 0
        throw(DomainError(n, "must be non-negative"))
    end
    if m < 0 || m > n
        throw(DomainError(m, "must be in the range m = 0,..,n"))
    end
    if x < -1 || x > 1
        throw(DomainError(x, "must be in the range [-1,1]"))
    end

    if m == 0
        return legendreP(n, x)
    end

    # step 1: compute P_m^m(x)
    p = 1                                           # = P_0^0(x)
    for k = 1:m
        p_prev  = p
        p       = -(2k-1) * sqrt(1-x^2) * p_prev    # = P_k^k(x)
    end                                             # on output: p = P_m^m(x)
    if n == m
        return p
    end

    # step 2: compute P_{m+1}^m(x)
    p_prev  = p                                     # = P_m^m(x)
    p       = (2m+1) * x * p_prev                   # = P_{m+1}^m(x)

    # step 3: compute P_n^m(x)
    for k = m+2:n
        p_prev, p_prev_prev = p, p_prev
        p = (2k-1)/(k-m)*x*p_prev - (k+m-1)/(k-m)*p_prev_prev       # = P_k^m(x)
    end                                             # on output: p = P_n^m(x)

    p
end

@doc raw"""
    legendreQ(n, x)

Evaluate the Legendre function of second kind ``Q_n(x)`` of integer degree ``n`` at position
``x``, defined by the Bonnet recursion

```math
Q_n(x)
=
\begin{cases}
    \frac{1}{2} \log \frac{1+x}{1-x}
    & \text{for} \quad n = 0\,, \; x \in (-1,1)
    \\
    P_1(x) Q_0(x) - 1
    & \text{for} \quad n = 1\,, \; x \in (-1,1)
    \\
    \frac{2n-1}{n} x Q_{n-1}(x) - \frac{n-1}{n} Q_{n-2}(x)
    & \text{for} \quad n \geq 2\,, \; x \in (-1,1)
\end{cases}
```
where ``P_1`` is the Legendre polynomial of first kind and degree one.
The function diverges at the interval boundaries
``\lim_{x\searrow -1} Q_n(x) = (-1)^{n+1} \infty\,, \lim_{x\nearrow 1} Q_n(x) = \infty``.

External links: [DLMF](https://dlmf.nist.gov/14.7.E2),
[Wikipedia](https://en.wikipedia.org/wiki/Legendre_function#Legendre_functions_of_the_second_kind_(Qn)).
"""
function legendreQ(n::Integer, x::Real)
    if n < 0
        throw(DomainError(n, "must be non-negative"))
    end
    if x < -1 || x > 1
        throw(DomainError(x, "must be in the range [-1,1]"))
    end

    if      x == 1
        return Inf
    elseif  x == -1
        return (-1)^(n+1) * Inf
    end

    Q0(y) = 0.5 * log((1+y)/(1-y))
    Q1(y) = legendreP(1, x) * Q0(y) - 1

    ABC_recurrence(n, x,
        m->(2m-1)//m, 0, m->(m-1)//m,       # A_m, B_m, C_m
        Q0, Q1)
end

@doc raw"""
    ABC_recurrence(n, x, A, B, C, p0, p1)

Evaluate polynomials by the recurrence
```math
p_n(x)
=
\begin{cases}
    p_0(x)
    & \text{for} \quad n = 0
    \\
    p_1(x)
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
- `n::Integer`:                     Polynomial order
- `x::Real`:                        evaluation point
- `A::Union{Integer,Function}`:     recurrence coefficient ``A_m``
- `B::Union{Integer,Function}`:     recurrence coefficient ``B_m``
- `C::Union{Integer,Function}`:     recurrence coefficient ``C_m``
- `p0::Function`:                   function of degree 0, i.e. ``p_0``
- `p1::Function`:                   function of degree 1, i.e. ``p_1``

#Implementation
For some polynomials the coefficients ``A_m, B_m, C_m`` are independent of ``m`` and for
others they aren't. Which is why we use `Union{Integer,Function}` and functions can be
easily created by anonymous functions.

Degree 0, and 1 are supplied as a function. Efficiency is not an issue as it is only invoked
once.

# Reference
> JIN, J. M.;
> JJIE, Zhang Shan.
> "Computation of special functions".
> Wiley, 1996.
> (p. 23)
"""
function ABC_recurrence(n::Integer, x::Real,
        A::Union{Integer,Function}, B::Union{Integer,Function}, C::Union{Integer,Function},
        p0::Function, p1::Function)

    # evaluate initial functions
    p_prev_prev = p0(x)
    p_prev      = p1(x)

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
function ABC_recurrence_step(m::Integer, x::Real,
    A::Integer, B::Integer, C::Integer,
    p_prev::Real, p_prev_prev::Real)

    (A*x + B) * p_prev - C * p_prev_prev
end
