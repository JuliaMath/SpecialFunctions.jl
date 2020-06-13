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
    2 T_m^2(x) - 1
    & \text{for} \quad n = 2m\,, \ m>0
    \\
    2 T_m(x) T_{m+1}(x) - x
    & \text{for} \quad n = 2m+1\,, \ m>0
    \,.
\end{cases}
```

External links: [DLMF](https://dlmf.nist.gov/18.3.T1),
[Wikipedia](https://en.wikipedia.org/wiki/Chebyshev_polynomials).
"""
function chebyshevt(n::Integer, x)
    (n < 0) && throw(DomainError(n, "must be non-negative"))

    if     n == 0
        return one(x)
    elseif n == 1
        return x
    elseif iseven(n)
        Tn2 = chebyshevt(n÷2, x)
        return 2*Tn2*Tn2-one(x)
    else
        Tn2   = chebyshevt(n÷2,   x)
        Tn2p1 = chebyshevt(n÷2+1, x)
        return 2*Tn2*Tn2p1-x
    end
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

    U0 = one(x)
    U1 = 2x

    if n == 0
        return U0
    elseif n == 1
        return U1
    else
        # Um = U_m(x), Um1 = U_{m-1}(x), Um2 = U_{m-2}(x)
        Um, Um1 = U1, U0

        for m = 2:n
            Um1, Um2 = Um, Um1
            Um = 2x*Um1-Um2
        end

        return Um
    end
end

@doc raw"""
    hermiteh(n, x)

Evaluate the (physicists') Hermite polynomials ``H_n(x)`` of order ``n`` at position ``x``,
defined by

```math
H_n(x)
= (-1)^n e^{x^2} \frac{\mathrm{d}^n}{\mathrm{d}x^n} e^{-x^2}
\quad \text{for} \quad
n = 0, 1, 2, \dots \,.
```

As a recurrence relation it is defined via
```math
H_n(x)
=
\begin{cases}
    1
    & \text{for} \quad n = 0
    \\
    2x
    & \text{for} \quad n = 1
    \\
    2x H_{n-1}(x) - 2 (n-1) H_{n-2}(x)
    & \text{for} \quad n \geq 2
    \,.
\end{cases}
```

External links: [DLMF](https://dlmf.nist.gov/18.3.T1),
[Wikipedia](https://en.wikipedia.org/wiki/Hermite_polynomials).
"""
function hermiteh(n::Integer, x)
    (n < 0) && throw(DomainError(n, "must be non-negative"))

    H0 = one(x)
    H1 = 2x

    if n == 0
        return H0
    elseif n == 1
        return H1
    else
        # Hm = H_m(x), Hm1 = H_{m-1}(x), Hm2 = H_{m-2}(x)
        Hm, Hm1 = H1, H0

        for m = 2:n
            Hm1, Hm2 = Hm, Hm1
            Hm = 2x*Hm1 - 2(m-1)Hm2
        end

        return Hm
    end
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

As a recurrence relation it is defined via
```math
L_n(x)
=
\begin{cases}
    1
    & \text{for} \quad n = 0
    \\
    1-x
    & \text{for} \quad n = 1
    \\
    \frac{2n-1-x}{n} L_{n-1}(x) - \frac{n-1}{n} L_{n-2}(x)
    & \text{for} \quad n \geq 2
    \,.
\end{cases}
```

External links: [DLMF](https://dlmf.nist.gov/18.3.T1),
[Wikipedia](https://en.wikipedia.org/wiki/Laguerre_polynomials).
"""
function laguerrel(n::Integer, x)
    (n < 0) && throw(DomainError(n, "must be non-negative"))

    L0 = one(x)
    L1 = 1-x

    if n == 0
        return L0
    elseif n == 1
        return L1
    else
        # Lm = L_m(x), Lm1 = L_{m-1}(x), Lm2 = L_{m-2}(x)
        Lm, Lm1 = L1, L0

        for m = 2:n
            Lm1, Lm2 = Lm, Lm1
            Lm = (2m-1-x)/m*Lm1 - (m-1)/m*Lm2
        end

        return Lm
    end
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

    P0 = one(x)
    P1 = x

    if n == 0
        return P0
    elseif n == 1
        return P1
    else
        # Pm = P_m(x), Pm1 = P_{m-1}(x), Pm2 = P_{m-2}(x)
        Pm, Pm1 = P1, P0

        for m = 2:n
            Pm1, Pm2 = Pm, Pm1
            Pm = (2m-1)/m*x*Pm1 - (m-1)/m*Pm2
        end

        return Pm
    end                                                # P_1(x)
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
    p = float(one(x))                                               # = P_0^0(x)
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

    if n == 0
        return Q0
    elseif n == 1
        return Q1
    else
        # Qm = Q_m(x), Qm1 = Q_{m-1}(x), Qm2 = Q_{m-2}(x)
        Qm, Qm1 = Q1, Q0

        for m = 2:n
            Qm1, Qm2 = Qm, Qm1
            Qm = (2m-1)/m*x*Qm1 - (m-1)/m*Qm2
        end

        return Qm
    end
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