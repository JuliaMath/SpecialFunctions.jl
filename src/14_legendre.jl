@doc raw"""
    chebyshevt(n, x)

Evaluate the Chebyshev polynomials of the first kind ``T_n(x)`` of order ``n`` at position
``x``, defined by

```math
T_n(x)
=
\begin{cases}
    1
    & \text{for} \quad n = 0
    \\
    x
    & \text{for} \quad n = 1
    \\
    2 x T_{n-1}(x) - T_{n-2}(x)
    & \text{for} \quad n \geq 2
    \,.
\end{cases}
```

External links: [DLMF](https://dlmf.nist.gov/18.3.T1),
[Wikipedia](https://en.wikipedia.org/wiki/Chebyshev_polynomials).
"""
function chebyshevt(n::Integer, x)
    (n < 0) && throw(DomainError(n, "must be non-negative"))

    ABC_recurrence(n, x,
        ( (2one(x), zero(x), one(x)) for m=2:n),            # (A_m, B_m, C_m)
        one(x),                                             # P_0(x)
        x)                                                  # P_1(x)
end

@doc raw"""
    chebyshevu(n, x)

Evaluate the Chebyshev polynomials of the second kind ``U_n(x)`` of order ``n`` at position
``x``, defined by

```math
U_n(x)
=
\begin{cases}
    1
    & \text{for} \quad n = 0
    \\
    2x
    & \text{for} \quad n = 1
    \\
    2x U_{n-1}(x) - U_{n-2}(x)
    & \text{for} \quad n \geq 2
    \,.
\end{cases}
```

External links: [DLMF](https://dlmf.nist.gov/18.3.T1),
[Wikipedia](https://en.wikipedia.org/wiki/Chebyshev_polynomials).
"""
function chebyshevu(n::Integer, x)
    (n < 0) && throw(DomainError(n, "must be non-negative"))

    ABC_recurrence(n, x,
        ( (2one(x), zero(x), one(x)) for m=2:n),            # (A_m, B_m, C_m)
        one(x),                                             # P_0(x)
        2x)                                                 # P_1(x)
end

@doc raw"""
    hermiteh(n, x)

Evaluate the (physicists') Hermite polynomials ``H_n(x)`` of order ``n`` at position ``x``,
defined by

```math
H_n(x)
= (-1)^n e^{x^2} \frac{\mathrm{d}^n}{\mathrm{d}x^n} e^{-x^2}
\quad \text{for} \quad
n = 0, 1, 2, \dots
```

External links: [DLMF](https://dlmf.nist.gov/18.3.T1),
[Wikipedia](https://en.wikipedia.org/wiki/Hermite_polynomials).
"""
function hermiteh(n::Integer, x)
    (n < 0) && throw(DomainError(n, "must be non-negative"))

    ABC_recurrence(n, x,
        ( (2one(x), zero(x), 2(m-one(x))) for m=2:n),               # (A_m, B_m, C_m)
        one(x),                                                     # P_0(x)
        2x)                                                         # P_1(x)
end

@doc raw"""
    laguerrel(n, x)

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
function laguerrel(n::Integer, x)
    (n < 0) && throw(DomainError(n, "must be non-negative"))

    ABC_recurrence(n, x,
        ( (-one(x)/m, 2-one(x)/m, 1-one(x)/m) for m=2:n),               # (A_m, B_m, C_m)
        one(x),                                                         # P_0(x)
        one(x)-x)                                                       # P_1(x)
end

@doc raw"""
    legendrep(n[, m], x)

Evaluate the Legendre polynomial ``P_n(x)`` of degree ``n`` at position ``x``.
If the order ``m`` is supplied then the associated Legendre polynomial ``P_n^{(m)}(x)`` is
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
The Legendre polynomial ``P_n(x)`` are also called Legendre polynomial of first kind,
respectively ``P_n^{(m)}(x)`` associated Legendre functions of first kind.

External links: [DLMF - Legendre polynomial](https://dlmf.nist.gov/14.7.E1),
[DLMF - associated Legendre polynomial](https://dlmf.nist.gov/14.7.E8),
[Wikipedia - Legendre polynomial](https://en.wikipedia.org/wiki/Legendre_polynomials),
[Wikipedia - associated Legendre polynomial](https://en.wikipedia.org/wiki/Associated_Legendre_polynomials).
"""
function legendrep(n::Integer, x)
    (n < 0) && throw(DomainError(n, "must be non-negative"))

    ABC_recurrence(n, x,
        ( ((2m-one(x))/m, zero(x), (m-one(x))/m) for m=2:n),            # (A_m, B_m, C_m)
        one(x),                                                         # P_0(x)
        x)                                                              # P_1(x)
end
function legendrep(n::Integer, m::Integer, x)

    # special case: m=0  =>  P_n^m = P_n
    # Legendre polynomial P_n(x) can be evaluated outside of [-1,1].
    # Thus, we use that routine before checking x interval for associated Legendre polynomial P_n^m(x).
    (0 == m) && return legendrep(n, x)

    # general argument checks
    (n < 0)          && throw(DomainError(n, "must be non-negative"))
    (m < 0 || m > n) && throw(DomainError(m, "must be in the range m = 0,..,n"))

    # step 1: compute P_m^m(x)
    p = one(x)                                                      # = P_0^0(x)
    for k = 1:m
        p_prev = p
        p      = -(2k-1) * sqrt(1-x^2) * p_prev                     # = P_k^k(x)
    end                                                             # on output: p = P_m^m(x)
    (n == m) && return p

    # step 2: compute P_{m+1}^m(x)
    p_prev = p                                                      # = P_m^m(x)
    p      = (2m+1) * x * p_prev                                    # = P_{m+1}^m(x)

    # step 3: compute P_n^m(x)
    for k = m+2:n
        p_prev, p_prev_prev = p, p_prev
        p = (2k-1)/(k-m)*x*p_prev - (k+m-1)/(k-m)*p_prev_prev       # = P_k^m(x)
    end                                                             # on output: p = P_n^m(x)

    p
end

@doc raw"""
    legendreq(n[, m], x)

Evaluate the Legendre function of second kind ``Q_n(x)`` of integer degree ``n`` at position
``x``.
If the order ``m`` is supplied then the associated Legendre function of second kind
``Q_n^{(m)}(x)`` is computed.
Their definitions are given by

```math
\begin{align}
    Q_n^{(m)}(x)
    &=
    (-1)^m (1-x^2)^\frac{m}{2} \frac{\mathrm{d}^m}{\mathrm{d}x^m} \left( Q_n(x) \right)
    \,, \quad m \in \mathbb{N}_0 \,, \; x \in (-1,1)
    \\
    Q_n(x)
    &=
    \begin{cases}
        \frac{1}{2} \log \frac{1+x}{1-x}
        & \text{for} \quad n = 0
        \\
        P_1(x) Q_0(x) - 1
        & \text{for} \quad n = 1
        \\
        \frac{2n-1}{n} x Q_{n-1}(x) - \frac{n-1}{n} Q_{n-2}(x)
        & \text{for} \quad n \geq 2
        \,.
    \end{cases}
\end{align}
```
where ``P_1`` is the Legendre polynomial of first kind and degree one.
The functions diverge at the interval boundaries
``\lim_{x\rightarrow \pm 1} |Q_n(x)| = \infty``.

External links: [DLMF - Legendre function](https://dlmf.nist.gov/14.7.E2),
[DLMF - associated Legendre function](https://dlmf.nist.gov/14.7.E9),
[Wikipedia - Legendre function](https://en.wikipedia.org/wiki/Legendre_functions#Legendre_functions_of_the_second_kind_(Qn)),
[Wikipedia - associated Legendre function](https://en.wikipedia.org/wiki/Legendre_function#Associated_Legendre_functions_of_the_second_kind).
"""
function legendreq(n::Integer, x)
    (0 > n)      && throw(DomainError(n, "must be non-negative"))
    (1 < abs(x)) && throw(DomainError(x, "must be in the range [-1,1]"))

    if x == 1
        return Inf
    elseif x == -1
        return iseven(n) ? -Inf : +Inf
    end

    Q0 = 0.5 * log((1+x)/(1-x))
    Q1 = legendrep(1, x) * Q0 - 1

    ABC_recurrence(n, x,
        ( ((2m-one(x))/m, zero(x), (m-one(x))/m) for m=2:n),            # (A_m, B_m, C_m)
        Q0, Q1)
end
function legendreq(n::Integer, m::Integer, x)
    (n <  0) && throw(DomainError(n, "must be non-negative"))
    (m <  0) && throw(DomainError(m, "must be non-negative"))
    (m == 0) && return legendreq(n, x)

    # x check after: Q^0_n is also implemented for x = +-1
    (1 <= abs(x)) && throw(DomainError(x, "must be in the range (-1,1)"))

    # step 1: compute Q_n^0(x)
    q_n_0 = legendreq(n, x)

    # step 2: compute Q_n^1(x)
    Q01 = -(1-x^2)^(-0.5)                                                   # = Q_0^1(x)
    Q11 = -sqrt(1-x^2) * (0.5*log((1+x)/(1-x)) + x/(1-x^2))                 # = Q_1^1(x)
    if n == 0
        q = Q01
    elseif n == 1
        q = Q11
    else
        q, q_prev = Q11, Q01
        for k = 2:n
            q_prev, q_prev_prev = q, q_prev
            q = (2k-1) / (k-1) * x * q_prev - k / (k-1) * q_prev_prev       # = Q_k^1(x)
        end                                                                 # on output: q = Q_n^1(x)
    end
    q_n_1 = q

    # step 3: compute Q_n^m(x)
    q_prev  = q_n_0
    q       = q_n_1
    for k = 2:m
        q_prev, q_prev_prev = q, q_prev
        q = -2*(k-1)*x*q_prev / sqrt(1-x^2) - (n+k-1)*(n-k+2)*q_prev_prev   # = Q_n^k(x)
    end                                                                     # on output: q = Q_n^m(x)

    q
end

@doc raw"""
    ABC_recurrence(n, x, ABC, p0, p1)

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
- `x`:                              evaluation point
- `ABC::Base.Generator`:            recurrence coefficients ``(A_m,B_m,C_m)`` for ``m=2,\dots,n``
- `p0`:                             function of degree 0, i.e. ``p_0``, evaluated at ``x``
- `p1`:                             function of degree 1, i.e. ``p_1``, evaluated at ``x``

#Implementation
The coefficients ``A_m, B_m, C_m`` are supplied via a generator which yields them as a tuple ``(A_m,B_m,C_m)``.

Type stability:
- The generator `ABC` yields tuples `(Am,Bm,Cm)`. The tuple's elements are of the same type as `x`.
- `p0,p1` are functions which have the same return type as the type of `x`.
- Thus, all variables: `p_prev_prev,p_prev,p,Am,Bm,Cm` are always of type `x`.

# Reference
> JIN, J. M.;
> JJIE, Zhang Shan.
> "Computation of special functions".
> Wiley, 1996.
> (p. 23)
"""
function ABC_recurrence(n::Integer, x::T,
        ABC::Base.Generator,
        p0::T, p1::T) where {T}

    # handle low orders
    if n == 0
        return p0
    elseif n == 1
        return p1
    end

    # use initial functions for recurrence relation
    p_prev_prev = p0
    p_prev      = p1

    p = zero(x)
    for AmBmCm in ABC
        Am, Bm, Cm = AmBmCm

        # recurrence step
        p = (Am*x + Bm) * p_prev - Cm * p_prev_prev

        p_prev_prev = p_prev
        p_prev      = p
    end

    p
end
