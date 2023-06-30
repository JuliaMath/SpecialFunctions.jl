using Base.Math: @horner
using Base.MPFR: ROUNDING_MODE
#useful constants
const acc0 = [5.0e-15, 5.0e-7, 5.0e-4] #accuracy options
const big1 = [25.0, 14.0, 10.0]
const e0 = [0.25e-3, 0.25e-1, 0.14]
const x0 = [31.0, 17.0, 9.7]
const alog10 = log(10)
const rt2pin = Float64(invsqrt2π)
const rtpi = Float64(sqrtπ)
const stirling_coef = [1.996379051590076518221, -0.17971032528832887213e-2, 0.131292857963846713e-4, -0.2340875228178749e-6, 0.72291210671127e-8, -0.3280997607821e-9, 0.198750709010e-10, -0.15092141830e-11, 0.1375340084e-12, -0.145728923e-13, 0.17532367e-14, -0.2351465e-15, 0.346551e-16, -0.55471e-17, 0.9548e-18, -0.1748e-18, 0.332e-19, -0.58e-20]
const auxgam_coef = [-1.013609258009865776949, 0.784903531024782283535e-1, 0.67588668743258315530e-2, -0.12790434869623468120e-2, 0.462939838642739585e-4, 0.43381681744740352e-5, -0.5326872422618006e-6, 0.172233457410539e-7, 0.8300542107118e-9, -0.10553994239968e-9, 0.39415842851e-11, 0.362068537e-13, -0.107440229e-13, 0.5000413e-15, -0.62452e-17, -0.5185e-18, 0.347e-19, -0.9e-21]

#----------------COEFFICIENTS FOR TEMME EXPANSION------------------

const d00 = -.333333333333333E+00
const d0  = [.833333333333333E-01, -.148148148148148E-01, .115740740740741E-02, .352733686067019E-03, -.178755144032922E-03, .391926317852244E-04]
const d10 = -.185185185185185E-02
const d1  = [-.347222222222222E-02, .264550264550265E-02, -.990226337448560E-03, .205761316872428E-03]
const d20 = .413359788359788E-02
const d2  = [-.268132716049383E-02, .771604938271605E-03]
const d30 = .649434156378601E-03
const d3  =[.229472093621399E-03, -.469189494395256E-03]
const d40 = -.861888290916712E-03
const d4  = [.784039221720067E-03]
const d50 = -.336798553366358E-03
const d5  = [-.697281375836586E-04]
const d60 = .531307936463992E-03
const d6  = [-.592166437353694E-03]
const d70 = .344367606892378E-03
const d80 = -.652623918595309E-03

"""
   rgamma1pm1(a)

   Computation of ``1/Gamma(a+1) - 1`` for `-0.5<=a<=1.5` : ``1/\\Gamma (a+1) - 1``
   Uses the relation `gamma(a+1) = a*gamma(a)`.
"""
function rgamma1pm1(a::Float64)
    @assert -0.5 <= a <= 1.5
    t = a
    rangereduce = a > 0.5
    t = rangereduce ? a-1 : a #-0.5<= t <= 0.5
    if t == 0.0
        return 0.0
    elseif t < 0.0
        top = @horner(t, -.422784335098468E+00, -.771330383816272E+00, -.244757765222226E+00, .118378989872749E+00, .930357293360349E-03, -.118290993445146E-01, .223047661158249E-02, .266505979058923E-03, -.132674909766242E-03)
        bot = @horner(t, 1.0, .273076135303957E+00, .559398236957378E-01)
        w = top/bot
        return rangereduce ? t*w/a : a*(w + 1)
    else
        top = @horner(t, .577215664901533E+00, -.409078193005776E+00, -.230975380857675E+00, .597275330452234E-01, .766968181649490E-02, -.514889771323592E-02, .589597428611429E-03)
        bot = @horner(t, 1.0, .427569613095214E+00, .158451672430138E+00, .261132021441447E-01, .423244297896961E-02)
        w = top/bot
        return rangereduce ? (t/a)*(w - 1.0) : a*w
    end
end

"""
    rgammax(a,x)

Evaluation of ``1/\\Gamma(a) e^{-x} x^{a}``. Based on `DRCOMP` from the `NSWC` library.
"""
function rgammax(a::Float64, x::Float64)
    if x == 0.0
        return 0.0
    elseif a > 20.0
        t = x/a
        if t == 0.0
            return 0.0
        end
        w = -(stirling_error(a) - a*LogExpFunctions.logmxp1(t))
        return 1/sqrt(2*Float64(π))*sqrt(a)*exp(w)
    else
        t = a*log(x) - x
        if a >= 1.0
            exp(t)/gamma(a)
        else
            return a*exp(t)*(1.0 + rgamma1pm1(a))
        end
    end
end

"""
    auxgam(x)

Compute function `g` in ``1/\\Gamma(x+1) = 1+x*(x-1)*g(x)``, `-1 <= x <= 1`.
"""
function auxgam(x::Float64)
    @assert -1 <= x <= 1
    if x < 0
        return -(1.0 + (1.0 + x)*(1.0 + x)*auxgam(1.0 + x))/(1.0 - x)
    else
        t = 2*x - 1.0
        return chepolsum(t, auxgam_coef)
    end
end

"""
    loggamma1p(x)

Compute ``log(\\Gamma(1+x))`` for `-1 < x <= 1`.
"""
function loggamma1p(x::Float64)
    @assert -1 < x <= 1
    return -log1p(x*(x - 1.0)*auxgam(x))
end

"""
    chepolsum(n,x,a)

Computes a series of Chebyshev Polynomials given by: `a[1]/2 + a[2]T1(x) + .... + a[n]T{n-1}(X)`.
"""
function chepolsum(x::Float64, a::Array{Float64,1})
    n = length(a)
    if n == 1
        return a[1]/2.0
    elseif n == 2
        return a[1]/2.0 + a[2]*x
    else
        tx = 2*x
        r = a[n]
        h = a[n - 1] + r*tx
        for k = n-2:-1:2
            s = r
            r = h
            h = a[k] + r*tx - s
        end
        return a[1]/2.0 - r + h*x
    end
end

"""
    stirling_error(x)

Compute ``\\ln{\\Gamma(x)} - (x-0.5)*\\ln{x} + x - \\ln{(2\\pi)}/2``. Adapted from `stirling` in `IncgamFI`.
"""
function stirling_error(x::Float64)
    if x < floatmin(Float64)*1000.0
        return floatmax(Float64)/1000.0
    elseif x < 1
        return loggamma1p(x) - (x + 0.5)*log(x) + x - log2π/2.0
    elseif x < 2
        return loggamma1p(x - 1) - (x - 0.5)*log(x) + x - log2π/2.0
    elseif x < 3
        return loggamma1p(x - 2) - (x - 0.5)*log(x) + x  - log2π/2.0 + log(x - 1)
    elseif x < 12
        z=18.0/(x*x)-1.0
        return chepolsum(z, stirling_coef)/(12.0*x)
    else
        z = 1.0/(x*x)
        if x < 1000
            return @horner(z, 0.25721014990011306473e-1, 0.82475966166999631057e-1, -0.25328157302663562668e-2, 0.60992926669463371e-3, -0.33543297638406e-3, 0.250505279903e-3)/((0.30865217988013567769 + z)*x)
        else
            return (((-z/1680.0 + 1.0/1260.0)*z - 1.0/360.0)*z + 1.0/12.0)/x
        end
    end
end

@doc raw"""
    gammax(x)

```math
\operatorname{gammax}(x) = \begin{cases}e^{\operatorname{stirling}(x)}\quad\quad\quad \text{if} \quad x>0,\\
\frac{\Gamma(x)}{\sqrt{2 \pi}e^{-x + (x-0.5)\operatorname{log}(x)}},\quad \text{if} \quad x\leq 0.
\end{cases}
```
"""
function gammax(x::Float64)
    if x >= 3
        return exp(stirling_error(x))
    elseif x > 0
        return gamma(x)/(exp(-x + (x - 0.5)*log(x))*sqrt2π)
    else
        return floatmax(Float64)/1000.0
    end
end

"""
    lambdaeta(eta)

Compute the value of ``\\lambda`` satisfying ``\\eta^{2}/2 = \\lambda-1-\\log{\\lambda}``.
"""
function lambdaeta(eta::Float64)
    s = eta*eta*0.5
    if eta == 0.0
        la = 1
    elseif eta < -1.0
        r = exp(-1 - s)
        la = @horner(r, 0.0, 1.0, 1.0, 1.5, 8.0/3.0, 125.0/24.0, 15.0/5.0)
    elseif eta < 1.0
        r = eta
        la = @horner(r, 1.0, 1.0, 1.0/3.0, 1.0/36.0, -1.0/270.0, 1.0/4320.0, 1.0/17010.0)
    else
        r = 11 + s
        l = log(r)
        la = r + l
        r = 1.0/r
        ak1 = 1.0
        ak2 = (2 - l)*0.5
        ak3 = (@horner(l, 6, -9, 2))/6.0
        ak4 = -(@horner(l, -12, 36, -22, 3))/12.0
        ak5 = (@horner(l, 60, -300, 350, -125, 12))/60.0
        ak6 = -(@horner(l, -120, 900, -1700, 1125, -274, 20))/120.0
        la = la + l*@horner(r,0.0, ak1, ak2, ak3, ak4, ak5, ak6)
    end
    r = 1
    if (eta > -3.5 && eta < -0.03) || (eta > 0.03 && eta < 40)
        r = 1
        q = la
        while r > 1.0e-8
            la = q*(s + log(q))/(q - 1.0)
            r = abs(q/la - 1)
            q = la
        end
    end
    return la
end

"""
Computing the first coefficient for the expansion :
```math
\\epsilon (\\eta_{0},a) = \\epsilon_{1} (\\eta_{0},a)/a + \\epsilon_{2} (\\eta_{0},a)/a^{2} + \\epsilon_{3} (\\eta_{0},a)/a^{3}
```
Refer Eqn (3.12) in the paper
"""
function coeff1(eta::Float64)
    if abs(eta) < 1.0
        coeff1 = @horner(
            eta,
            -3.333333333438e-1, -2.070740359969e-1, -5.041806657154e-2,
            -4.923635739372e-3, -4.293658292782e-5
        ) / @horner(
            eta,
            1.000000000000e+0, 7.045554412463e-1, 2.118190062224e-1,
            3.048648397436e-2, 1.605037988091e-3
        )
    else
        la = lambdaeta(eta)
        coeff1 = log(eta/(la - 1.0))/eta
    end
    return coeff1
end
"""
Computing the second coefficient for the expansion :
```math
\\epsilon (\\eta_{0},a) = \\epsilon_{1} (\\eta_{0},a)/a + \\epsilon_{2} (\\eta_{0},a)/a^{2} + \\epsilon_{3} (\\eta_{0},a)/a^{3}
```
Refer Eqn (3.12) in the paper
"""
function coeff2(eta::Float64)

    if eta < -5.0
        x = eta*eta
        lnmeta = log(-eta)
        coeff2 = (12.0 - x - 6.0*lnmeta*lnmeta)/(12.0*x*eta)
    elseif eta < -2.0
        coeff2 = @horner(
            eta,
            -1.72847633523e-2, -1.59372646475e-2, -4.64910887221e-3,
            -6.06834887760e-4, -6.14830384279e-6
        ) / @horner(
            eta,
            1.00000000000e+0, 7.64050615669e-1, 2.97143406325e-1,
            5.79490176079e-2, 5.74558524851e-3)
    elseif eta < 2.0
        coeff2 = @horner(
            eta,
            -1.72839517431e-2, -1.46362417966e-2, -3.57406772616e-3,
            -3.91032032692e-4, 2.49634036069e-6
        ) / @horner(
            eta,
            1.00000000000e+0, 6.90560400696e-1, 2.49962384741e-1,
            4.43843438769e-2, 4.24073217211e-3
        )
    elseif eta < 1000.0
        x = 1.0/eta
        coeff2 = @horner(
            x,
            9.99944669480e-1, 1.04649839762e+2, 8.57204033806e+2,
            7.31901559577e+2, 4.55174411671e+1
        ) / @horner(
            x,
            1.00000000000e+0, 1.04526456943e+2, 8.23313447808e+2,
            3.11993802124e+3, 3.97003311219e+3
        )/(-12.0*eta)
    else
        coeff2 = -1.0/(12.0*eta)
    end
    return coeff2
end
"""
Computing the third and last coefficient for the expansion :
```math
\\epsilon (\\eta_{0},a) = \\epsilon_{1} (\\eta_{0},a)/a + \\epsilon_{2} (\\eta_{0},a)/a^{2} + \\epsilon_{3} (\\eta_{0},a)/a^{3}
```
Refer Eqn (3.12) in the paper
"""
function coeff3(eta::Float64)
    if eta < -8.0
        x=eta*eta
        y = log(-eta)/eta
        coeff3=(-30.0 + eta*y*(6.0*x*y*y - 12.0 + x))/(12.0*eta*x*x)
    elseif eta < -4.0
        coeff3 = (
            @horner(
                eta,
                4.95346498136e-2, 2.99521337141e-2, 6.88296911516e-3,
                5.12634846317e-4, -2.01411722031e-5
            ) / @horner(
                eta,
                1.00000000000e+0, 7.59803615283e-1, 2.61547111595e-1,
                4.64854522477e-2, 4.03751193496e-3
            )
        )/(eta*eta)
    elseif eta < -2.0
        coeff3 = @horner(
            eta,
            4.52313583942e-3, 1.20744920113e-3, -7.89724156582e-5,
            -5.04476066942e-5, -5.35770949796e-6
        ) / @horner(
            eta,
            1.00000000000e+0, 9.12203410349e-1, 4.05368773071e-1,
            9.01638932349e-2, 9.48935714996e-3
        )
    elseif eta < 2.0
        coeff3 = @horner(
            eta,
            4.39937562904e-3, 4.87225670639e-4, -1.28470657374e-4,
            5.29110969589e-6, 1.57166771750e-7
        ) / @horner(
            eta,
            1.00000000000e+0, 7.94435257415e-1, 3.33094721709e-1,
            7.03527806143e-2, 8.06110846078e-3
        )
    elseif eta < 10.0
        x = 1.0/eta
        coeff3 = (
            @horner(
                x,
                -1.14811912320e-3, -1.12850923276e-1, 1.51623048511e+0,
                -2.18472031183e-1, 7.30002451555e-2
            ) / @horner(
                x,
                1.00000000000e+0, 1.42482206905e+1, 6.97360396285e+1,
                2.18938950816e+2, 2.77067027185e+2
            )
        )/(eta*eta)
    elseif eta < 100.0
        x = 1.0/eta
        coeff3 = (
            @horner(
                x,
                -1.45727889667e-4, -2.90806748131e-1, -1.33085045450e+1,
                1.99722374056e+2, -1.14311378756e+1
            ) / @horner(
                x,
                1.00000000000e+0, 1.39612587808e+2, 2.18901116348e+3,
                7.11524019009e+3, 4.55746081453e+4
            )
        )/(eta*eta)
    else
        eta3 = eta*eta*eta
        coeff3 = -log(eta)/(12.0*eta3)
    end
    return coeff3
end


"""
    gamma_inc_cf(a, x, ind)

Computes ``P(a,x)`` by continued fraction expansion given by :
```math
R(a,x) * \\frac{1}{1-\\frac{z}{a+1+\\frac{z}{a+2-\\frac{(a+1)z}{a+3+\\frac{2z}{a+4-\\frac{(a+2)z}{a+5+\\frac{3z}{a+6-\\dots}}}}}}}
```
Used when `1 <= a <= BIG` and `x < x0`.

External links: [DLMF](https://dlmf.nist.gov/8.9.2)

See also: [`gamma_inc(a,x,ind)`](@ref SpecialFunctions.gamma_inc)
"""
function gamma_inc_cf(a::Float64, x::Float64, ind::Integer)
    acc = acc0[ind + 1]
    tol = 4.0*acc
    a2nm1 = 1.0
    a2n = 1.0
    b2nm1 = x
    b2n = x + (1.0 - a)
    c = 1.0
    while true
       a2nm1 = x*a2n + c*a2nm1
       b2nm1 = x*b2n + c*b2nm1
       c = c + 1.0
       t = c - a
       a2n = a2nm1 + t*a2n
       b2n = b2nm1 + t*b2n
       a2nm1 = a2nm1/b2n
       b2nm1 = b2nm1/b2n
       a2n = a2n/b2n
       b2n = 1.0
       if abs(a2n - a2nm1/b2nm1) < tol*a2n
           break
       end
    end
    q = rgammax(a, x)*a2n
    return (1.0 - q, q)
end
"""
    gamma_inc_taylor(a, x, ind)

Compute ``P(a,x)`` using Taylor Series for P/R given by :
```math
R(a,x)/a * (1 + \\sum_{n=1}^{\\infty} x^{n}/((a+1)(a+2)...(a+n)))
```
Used when `1 <= a <= BIG` and `x <= max{a, ln 10}`.

External links: [DLMF](https://dlmf.nist.gov/8.11.2)

See also: [`gamma_inc(a,x,ind)`](@ref SpecialFunctions.gamma_inc)
"""
function gamma_inc_taylor(a::Float64, x::Float64, ind::Integer)
    acc = acc0[ind + 1]
    wk = zeros(30)
    flag = false
    apn = a + 1.0
    t = x/apn
    wk[1] = t
    loop = 2
    for indx = 2:20
       apn += 1.0
       t *= x/apn
       if t <= 1.0e-3
           loop = indx
           flag = true
           break
       end
       wk[indx] = t
    end
    if !flag
        loop = 20
    end
    sm = t
    tol = 0.5*acc #tolerance
    while true
       apn += 1.0
       t *= x/apn
       sm += t
       if t <= tol
           break
       end
    end
    for j = loop-1:-1:1
       sm += wk[j]
    end
    p = (rgammax(a, x)/a)*(1.0 + sm)
    return (p, 1.0 - p)
end
"""
    gamma_inc_asym(a, x, ind)

Compute ``P(a,x)`` using asymptotic expansion given by:
```math
R(a,x)/a * (1 + \\sum_{n=1}^{N-1}(a_{n}/x^{n} + \\Theta _{n}a_{n}/x^{n}))
```
where `R(a,x) = rgammax(a,x)`. Used when `1 <= a <= BIG` and `x >= x0`.

External links: [DLMF](https://dlmf.nist.gov/8.11.2)

See also: [`gamma_inc(a,x,ind)`](@ref SpecialFunctions.gamma_inc)
"""
function gamma_inc_asym(a::Float64, x::Float64, ind::Integer)
    wk = zeros(30)
    flag = false
    acc = acc0[ind + 1]
    amn = a - 1.0
    t = amn/x
    wk[1] = t
    loop = 2
    for indx = 2:20
       amn -= 1.0
       t *= amn/x
       if abs(t) <= 1.0e-3
           loop = indx
           flag = true
           break
       end
       wk[indx] = t
    end
    if !flag
        loop = 20
    end
    sm = t
    while true
       if abs(t) < acc
           break
       end
       amn -= 1.0
       t *= amn/x
       sm += t
    end
    for j = loop-1:-1:1
       sm += wk[j]
    end
    q = (rgammax(a, x)/x)*(1.0 + sm)
    return (1.0 - q, q)
end
"""
    gamma_inc_taylor_x(a,x,ind)

Computes ``P(a,x)`` based on Taylor expansion of ``P(a,x)/x**a`` given by:
```math
J = -a * \\sum_{1}^{\\infty} (-x)^{n}/((a+n)n!)
```
and ``P(a,x)/x**a`` is given by :
```math
(1 - J)/ \\Gamma(a+1)
```
resulting from term-by-term integration of `gamma_inc(a,x,ind)`.
This is used when `a < 1` and `x < 1.1` - Refer Eqn (9) in the paper.

See also: [`gamma_inc(a,x,ind)`](@ref SpecialFunctions.gamma_inc)
"""
function gamma_inc_taylor_x(a::Float64, x::Float64, ind::Integer)
    acc = acc0[ind + 1]
    l = 3.0
    c = x
    sm = x/(a + 3.0)
    tol = 3.0*acc/(a + 1.0)
    while true
       l += 1.0
       c *= -x/l
       t = c/(a + l)
       sm += t
       if abs(t) <= tol
           break
       end
    end
    temp = a*x*((sm/6.0 - 0.5/(a + 2.0))*x + 1.0/(a + 1.0))
    z = a*log(x)
    h = rgamma1pm1(a)
    g = 1.0 + h
    if (x < 0.25 && z > -.13394) || a < x/2.59
       l = expm1(z)
       w = 1.0 + l
       q = max((w*temp - l)*g - h, 0.0)
       return (1.0 - q, q)
    else
       w = exp(z)
       p = w*g*(1.0 - temp)
       return (p, 1.0 - p)
    end
end

"""
    gamma_inc_minimax(a,x,z)

Compute ``P(a,x)`` using minimax approximations given by :
```math
1/2 * erfc(\\sqrt{y}) - e^{-y}/\\sqrt{2\\pi*a}* T(a,\\lambda)
``` where
```math
T(a,\\lambda) = \\sum_{0}^{N} c_{k}(z)a^{-k}
```

This is a higher accuracy approximation of Temme expansion, which deals with the region near `a ≈ x` with `a` large.
Refer Appendix F in the paper for the extensive set of coefficients calculated using Brent's multiple precision arithmetic(set at 50 digits) in BRENT, R. P. A FORTRAN multiple-precision arithmetic package, ACM Trans. Math. Softw. 4(1978), 57-70 .

External links: [DLMF](https://dlmf.nist.gov/8.12.8)

See also: [`gamma_inc(a,x,ind)`](@ref SpecialFunctions.gamma_inc)
"""
function gamma_inc_minimax(a::Float64, x::Float64, z::Float64)
    l = x/a
    s = 1.0 - l
    y = -a*LogExpFunctions.logmxp1(l)
    c = exp(-y)
    w = 0.5*erfcx(sqrt(y))

    if abs(s) <= 1.0e-3
        c0 = @horner(z, d00, d0[1], d0[2], d0[3])
        c1 = @horner(z, d10, d1[1], d1[2], d1[3])
        c2 = @horner(z, d20, d2[1], d2[2])
        c3 = @horner(z, d30, d3[1], d3[2])
        c4 = @horner(z, d40, d4[1])
        c5 = @horner(z, d50, d5[1])
        c6 = @horner(z, d60, d6[1])

        t = @horner(z, c0, c1, c2, c3, c4, c5, c6, d70, d80)
        if l < 1.0
            p = c*(w - rt2pin*t/sqrt(a))
            return (p, 1.0 - p)
        else
            q = c*(w + rt2pin*t/sqrt(a))
            return (1.0 - q, q)
        end
    end
    #---USING THE MINIMAX APPROXIMATIONS---
    c0 = @horner(z, -.333333333333333E+00, -.159840143443990E+00, -.335378520024220E-01, -.231272501940775E-02)/(@horner(z, 1.0, .729520430331981E+00, .238549219145773E+00,  .376245718289389E-01, .239521354917408E-02, -.939001940478355E-05, .633763414209504E-06))
    c1 = @horner(z, -.185185185184291E-02, -.491687131726920E-02, -.587926036018402E-03, -.398783924370770E-05)/(@horner(z, 1.0, .780110511677243E+00, .283344278023803E+00,  .506042559238939E-01, .386325038602125E-02))
    c2 = @horner(z,  .413359788442192E-02,  .669564126155663E-03)/(@horner(z, 1.0, .810647620703045E+00, .339173452092224E+00, .682034997401259E-01, .650837693041777E-02, -.421924263980656E-03))
    c3 = @horner(z,  .649434157619770E-03,  .810586158563431E-03)/(@horner(z, 1.0, .894800593794972E+00, .406288930253881E+00, .906610359762969E-01, .905375887385478E-02, -.632276587352120E-03))
    c4 = @horner(z, -.861888301199388E-03, -.105014537920131E-03)/(@horner(z, 1.0, .103151890792185E+01, .591353097931237E+00, .178295773562970E+00, .322609381345173E-01))
    c5 = @horner(z, -.336806989710598E-03, -.435211415445014E-03)/(@horner(z, 1.0, .108515217314415E+01, .600380376956324E+00, .178716720452422E+00))
    c6 = @horner(z,  .531279816209452E-03, -.182503596367782E-03)/(@horner(z, 1.0, .770341682526774E+00, .345608222411837E+00))
    c7 = @horner(z,  .344430064306926E-03,  .443219646726422E-03)/(@horner(z, 1.0, .115029088777769E+01, .821824741357866E+00))
    c8 = @horner(z, -.686013280418038E-03,  .878371203603888E-03)

    t = @horner(1.0/a, c0, c1, c2, c3, c4, c5, c6, c7, c8)
    if l < 1.0
        p = c*(w - rt2pin*t/sqrt(a))
        return (p, 1.0 - p)
    else
        q = c*(w + rt2pin*t/sqrt(a))
        return (1.0 - q, q)
    end
end

"""
    gamma_inc_temme(a, x, z)

Compute ``P(a,x)`` using Temme's expansion given by :
```math
1/2 * erfc(\\sqrt{y}) - e^{-y}/\\sqrt{2\\pi*a}* T(a,\\lambda)
``` where
```math
T(a,\\lambda) = \\sum_{0}^{N} c_{k}(z)a^{-k}
```
This mainly solves the problem near the region when `a ≈ x` with a large, and is a lower accuracy version of the minimax approximation.

External links: [DLMF](https://dlmf.nist.gov/8.12.8)

See also: [`gamma_inc(a,x,ind)`](@ref SpecialFunctions.gamma_inc)
"""
function gamma_inc_temme(a::Float64, x::Float64, z::Float64)
    l = x/a
    y = -a*LogExpFunctions.logmxp1(x/a)
    c = exp(-y)
    w = 0.5*erfcx(sqrt(y))
    c0 = @horner(z, d00, d0[1], d0[2], d0[3], d0[4], d0[5], d0[6])
    c1 = @horner(z, d10, d1[1], d1[2], d1[3], d1[4])
    c2 = @horner(z, d20, d2[1])
    t  = @horner(1.0/a, c0, c1, c2)
    if l < 1.0
        p = c*(w - rt2pin*t/sqrt(a))
        return (p, 1.0 - p)
    else
        q = c*(w + rt2pin*t/sqrt(a))
        return (1.0 - q, q)
    end
end

"""
    gamma_inc_temme_1(a, x, z, ind)

Computes ``P(a,x)`` using simplified Temme expansion near ``y=0`` by :
```math
E(y) - (1 - y)/\\sqrt{2\\pi*a} * T(a,\\lambda)
```
where
```math
E(y) = 1/2 - (1 - y/3)*(\\sqrt(y/\\pi))
```
Used instead of it's previous function when ``\\sigma <= e_{0}/\\sqrt{a}``.

External links: [DLMF](https://dlmf.nist.gov/8.12.8)
"""
function gamma_inc_temme_1(a::Float64, x::Float64, z::Float64, ind::Integer)
    iop = ind + 1
    l = x/a
    y = -a * LogExpFunctions.logmxp1(l)
    if a*eps()*eps() > 3.28e-3
        throw(DomainError((a, x, ind), "P(a,x) or Q(a,x) is computationally indeterminant in this case."))
    end
    c = exp(-y)
    w = 0.5*erfcx(sqrt(y))
    u = 1.0/a
    if iop == 1
        c0 = @horner(z, d00, d0[1], d0[2], d0[3])
        c1 = @horner(z, d10, d1[1], d1[2], d1[3])
        c2 = @horner(z, d20, d2[1], d2[2])
        c3 = @horner(z, d30, d3[1], d3[2])
        c4 = @horner(z, d40, d4[1])
        c5 = @horner(z, d50, d5[1])
        c6 = @horner(z, d60, d6[1])

        t = @horner(u, c0, c1, c2, c3, c4, c5, c6, d70, d80)

    elseif iop == 2
        c0 = @horner(d00, d0[1], d0[2])
        c1 = @horner(d10, d1[1])
        t = @horner(u, c0, c1, d20)

    else
        t = @horner(z, d00, d0[1])

    end
    if l < 1.0
        p = c*(w - rt2pin*t/sqrt(a))
        return (p, 1.0 - p)
    else
        q = c*(w + rt2pin*t/sqrt(a))
        return (1.0 - q, q)
    end
end

"""
    gamma_inc_fsum(a,x)

Compute using Finite Sums for ``Q(a,x)`` when `a >= 1 && 2a` is integer.
Used when `a <= x <= x0` and `a = n/2`.
Refer Eqn (14) in the paper.

See also: [`gamma_inc(a,x,ind)`](@ref SpecialFunctions.gamma_inc)
"""
function gamma_inc_fsum(a::Float64, x::Float64)
    if isinteger(a)
        sm = exp(-x)
        t = sm
        N = 1
        c = 0.0
        i = trunc(Int, a)
    else
        rtx = sqrt(x)
        sm = erfc(rtx)
        t = exp(-x)/(rtpi*rtx)
        N = 0
        c = -0.5
        i = trunc(Int, a)
    end
    for lp = N:(i - 1)
        if i == 0
            break
        end
        c += 1.0
        t = (x*t)/c
        sm += t
    end
    q = sm
    return (1.0 - q, q)

end

"""
    gamma_inc_inv_psmall(a, logr)

Compute `x0` - initial approximation when `p` is small.
Here we invert the series in Eqn (2.20) in the paper and write the inversion problem as:
```math
x = r\\left[1 + a\\sum_{k=1}^{\\infty}\\frac{(-x)^{n}}{(a+n)n!}\\right]^{-1/a},
```
where ``r = (p\\Gamma(1+a))^{1/a}``
Inverting this relation we obtain ``x = r + \\sum_{k=2}^{\\infty}c_{k}r^{k}``.
"""
function gamma_inc_inv_psmall(a::Float64, logr::Float64)
    r    = exp(logr)
    ap1  = a + 1.0
    ap1² = ap1*ap1
    ap1³ = ap1*ap1²
    ap1⁴ = ap1²*ap1²
    ap2  = a + 2.0
    ap2² = ap2*ap2
    ck1  = 1.0
    ck2  = 1.0/(1.0 + a)
    ck3  = 0.5*(3*a + 5)/(ap1²*(a + 2))
    ck4  = (1.0/3.0)*(@horner(a, 31, 33, 8))/(ap1³*ap2*(a + 3))
    ck5  = (1.0/24.0)*(@horner(a, 2888, 5661, 3971, 1179, 125))/(ap1⁴*ap2²*(a + 3)*(a + 4))
    x0   = @horner(r, 0.0, ck1, ck2, ck3, ck4, ck5)
    return x0
end

"""
    gamma_inc_inv_qsmall(a, q, qgammaxa)

Compute `x0` - initial approximation when `q` is small from ``e^{-x_{0}} x_{0}^{a} = q \\Gamma(a)``.
Asymptotic expansions Eqn (2.29) in the paper is used here and higher approximations are obtained using
```math
x \\sim x_{0} - L + b \\sum_{k=1}^{\\infty} d_{k}/x_{0}^{k}
```
where ``b = 1-a``, ``L = \\ln{x_0}``.
"""
function gamma_inc_inv_qsmall(a::Float64, q::Float64, qgammaxa::Float64)
    b = 1.0 - a
    eta = sqrt(-2/a*log(qgammaxa))
    x0 = a*lambdaeta(eta)
    l = log(x0)

    if a > 0.12 || x0 > 5
        r = 1.0/x0
        ck1 = l - 1.0
        ck2 = (@horner(l, @horner(b, 2, 3), @horner(b, -2, -2), 1))/2.0
        ck3 = (@horner(l, @horner(b, -12, -24, -11), @horner(b, 12, 24, 6), @horner(b, -6, -9), 2))/6.0
        ck4 = (@horner(l, @horner(b, 72, 162, 120, 25), @horner(b, -72, -168, -114, -12), @horner(b, 36, 84, 36), @horner(b, -12, -22), 3))/12.0
        x0 = x0 - l + b*r*@horner(r, ck1, ck2, ck3, ck4)
    elseif x0 > 1
        # The x0 > 1 condition isn't in the original version but without it
        # the update in the branch can cause negative initial x0
        r = 1.0/x0
        l² = l*l
        ck1 = l - 1.0
        x0 = x0 - l + b*r*ck1
    end
    return x0
end

"""
    gamma_inc_inv_alarge(a, minpq, pcase)

Compute `x0` - initial approximation when `a` is large.
The inversion problem is rewritten as :
```math
0.5 \\operatorname{erfc}(\\eta \\sqrt{a/2}) + R_{a}(\\eta) = q
```
For large values of `a` we can write: ``\\eta(q,a) = \\eta_{0}(q,a) + \\epsilon(\\eta_{0},a)``
and it is possible to expand:
```math
\\epsilon(\\eta_{0},a) = \\epsilon_{1}(\\eta_{0},a)/a + \\epsilon_{2}(\\eta_{0},a)/a^{2} + \\epsilon_{3}(\\eta_{0},a)/a^{3} + ...
```
which is calculated by coeff1, coeff2 and coeff3 functions below.
This returns a tuple `(x0,fp)`, where `fp` is computed since it's an approximation for the coefficient after inverting the original power series.
"""
function gamma_inc_inv_alarge(a::Float64, minpq::Float64, pcase::Bool)
    r = erfcinv(2*minpq)
    s = r/sqrt(a*0.5)
    eta = pcase ? -s : s
    eta += (coeff1(eta) + (coeff2(eta) + coeff3(eta)/a)/a)/a
    x0 = a*lambdaeta(eta)
    fp = -sqrt(inv2π*a)*exp(-0.5*a*eta*eta)/gammax(a)
    return (x0, fp)
end
# Reference : 'Computation of the incomplete gamma function ratios and their inverse' by Armido R DiDonato, Alfred H Morris.
# Published in Journal: ACM Transactions on Mathematical Software (TOMS)
# Volume 12 Issue 4, Dec. 1986 Pages 377-393
# doi>10.1145/22721.23109

"""
    gamma_inc(a,x,IND=0)

Returns a tuple ``(p, q)`` where ``p + q = 1``, and
``p=P(a,x)`` is the Incomplete gamma function ratio given by:
```math
P(a,x)=\\frac{1}{\\Gamma (a)} \\int_{0}^{x} e^{-t}t^{a-1} dt.
```
and ``q=Q(a,x)`` is the Incomplete gamma function ratio given by:
```math
Q(a,x)=\\frac{1}{\\Gamma (a)} \\int_{x}^{\\infty} e^{-t}t^{a-1} dt.
```
In terms of these, the lower incomplete gamma function is
``\\gamma(a,x) = P(a,x) \\Gamma(a)`` and the upper incomplete
gamma function is ``\\Gamma(a,x) = Q(a,x) \\Gamma(a)``.

`IND ∈ [0,1,2]` sets accuracy: `IND=0` means 14 significant digits accuracy, `IND=1` means 6 significant digit, and `IND=2` means only 3 digit accuracy.

External links: [DLMF](https://dlmf.nist.gov/8.2.4), [Wikipedia](https://en.wikipedia.org/wiki/Incomplete_gamma_function)

See also [`gamma(z)`](@ref SpecialFunctions.gamma), [`gamma_inc_inv(a,p,q)`](@ref SpecialFunctions.gamma_inc_inv)
"""
gamma_inc(a::Real,x::Real,ind::Integer=0) = _gamma_inc(promote(float(a),float(x))...,ind)

function _gamma_inc(a::Float64, x::Float64, ind::Integer)
    iop = ind + 1
    acc = acc0[iop]
    if a < 0.0 || x < 0.0
        throw(DomainError((a, x, ind), "`a` and `x` must be greater than 0 ---- Domain : (0, Inf)"))
    elseif a == 0.0 && x == 0.0
        throw(DomainError((a, x, ind), "`a` and `x` must be greater than 0 ---- Domain : (0, Inf)"))
    elseif isnan(a) || isnan(x)
        ax = a*x
        return (ax, ax)
    elseif a == 0.0 || isinf(x)
        return (1.0, 0.0)
    elseif x == 0.0
        return (0.0, 1.0)
    end

    if a >= 1.0
        if a >= big1[iop]
            l = x/a
            if l == 0.0
                return (0.0, 1.0)
            end
            s = 1.0 - l
            z = -LogExpFunctions.logmxp1(l)
            if z >= 700.0/a
                if abs(s) <= 2*eps(Float64)
                    throw(DomainError((a, x, ind), "P(a,x) or Q(a,x) is computationally indeterminant in this case."))
                end
                if x <= a
                    return (0.0, 1.0)
                else
                    return (1.0, 0.0)
                end
            end

            y = a*z
            rta = sqrt(a)
            if abs(s) <= e0[iop]/rta
                z = sqrt(z + z)
                if l < 1.0
                    z = -z
                end
                return gamma_inc_temme_1(a, x, z, ind)
            end

            if abs(s) <= 0.4
                if abs(s) <= 2.0*eps() && a*eps()*eps() > 3.28e-3
                    throw(DomainError((a, x, ind), "P(a,x) or Q(a,x) is computationally indeterminant in this case."))
                end
                c = exp(-y)
                w = 0.5*erfcx(sqrt(y))
                u = 1.0/a
                z = sqrt(z + z)
                if l < 1.0
                    z = -z
                end
                if iop == 1
                    return gamma_inc_minimax(a, x, z)
                elseif iop == 2
                    return gamma_inc_temme(a, x, z)
                else
                    t = @horner(z, d00, d0[1], d0[2], d0[3])
                    return gamma_inc_temme_1(a, x, z, ind)
                end
            else
                return _gamma_inc_choose_algorithm(a, x, ind)
            end
        elseif a > x || x >= x0[iop] || !isinteger(2*a)
            return _gamma_inc_choose_algorithm(a, x, ind)
        else
            return gamma_inc_fsum(a,x)
        end
    elseif a == 0.5
        if x >= 0.25
            q = erfc(sqrt(x))
            return (1.0 - q, q)
        end
        p = erf(sqrt(x))
        return (p, 1.0 - p)
    elseif x < 1.1
        return gamma_inc_taylor_x(a, x, ind)
    end
    r = rgammax(a, x)
    if r == 0.0
        return (1.0, 0.0)
    else
        return gamma_inc_cf(a, x, ind)
    end
end

function _gamma_inc(a::BigFloat,x::BigFloat,ind::Integer) #BigFloat version from GNU MPFR wrapped via ccall
    z = BigFloat()
    ccall((:mpfr_gamma_inc, :libmpfr), Int32, (Ref{BigFloat}, Ref{BigFloat}, Ref{BigFloat}, Int32), z, a, x, ROUNDING_MODE[])
    q = z/gamma(a)
    return (1.0 - q, q)
end
_gamma_inc(a::Float32,x::Float32,ind::Integer) = Float32.(_gamma_inc(Float64(a),Float64(x),ind))
_gamma_inc(a::Float16,x::Float16,ind::Integer) = Float16.(_gamma_inc(Float64(a),Float64(x),ind))

function _gamma_inc_choose_algorithm(a::Float64, x::Float64, ind::Int)
    r = rgammax(a, x)
    if r == 0.0
        if x <= a
            return (0.0, 1.0)
        else
            return (1.0, 0.0)
        end
    end
    if x <= max(a, alog10)
        return gamma_inc_taylor(a, x, ind)
    elseif x < x0[ind + 1]
        return gamma_inc_cf(a, x, ind)
    else
        return gamma_inc_asym(a, x, ind)
    end
end

#EFFICIENT AND ACCURATE ALGORITHMS FOR THECOMPUTATION AND INVERSION OF THE INCOMPLETEGAMMA FUNCTION RATIOS by Amparo Gil, Javier Segura, Nico M. Temme
#SIAM Journal on Scientific Computing 34(6) (2012), A2965-A2981
# arXiv:1306.1754

"""
    gamma_inc_inv(a,p,q)

Inverts the `gamma_inc(a,x)` function, by computing `x` given `a`,`p`,`q` in ``P(a,x)=p`` and ``Q(a,x)=q``.

External links: [DLMF](https://dlmf.nist.gov/8.2.4), [Wikipedia](https://en.wikipedia.org/wiki/Incomplete_gamma_function)

See also: [`gamma_inc(a,x,ind)`](@ref SpecialFunctions.gamma_inc).
"""
function gamma_inc_inv(a::Real, p::Real, q::Real)
    return _gamma_inc_inv(map(float, promote(a, p, q))...)
end

# `gamma inc_inv` ensures that arguments of `_gamma_inc_inv` are
# floating point numbers of the same type
function _gamma_inc_inv(a::T, p::T, q::T) where {T<:Real}
    if p + q != 1
        throw(ArgumentError("p + q must equal one but is $(p + q)"))
    end

    if iszero(p)
        return zero(T)
    elseif iszero(q)
        return T(Inf)
    end

    pcase = p < 0.5
    minpq = pcase ? p : q
    return __gamma_inc_inv(a, minpq, pcase)
end

function __gamma_inc_inv(a::Float64, minpq::Float64, pcase::Bool)
    haseta = false

    logp = pcase ? log(minpq) : log1p(-minpq)
    loggamma1pa = a <= 1.0 ? loggamma1p(a) : loggamma(a + 1.0)
    logr = (logp + loggamma1pa) / a
    if logr < log(0.2*(1 + a)) #small value of p
        x0 = gamma_inc_inv_psmall(a, logr)
    elseif !pcase && a < 10 && minpq < 0.02 && (qgammaxa = minpq*gammax(a)*sqrt(twoπ/a)) < 1 #small q
        # This deviates from the original version. The qgammaxa variable
        # here ensures that the argument of sqrt in gamma_inc_inv_qsmall
        # is positive
        x0 = gamma_inc_inv_qsmall(a, minpq, qgammaxa)
    elseif abs(minpq - 0.5) < 1.0e-05
        x0 = a - 1.0/3.0 + (8.0/405.0 + 184.0/25515.0/a)/a
    elseif abs(a - 1.0) < 1.0e-4
        x0 = pcase ? -log1p(-minpq) : -log(minpq)
    elseif a < 1.0 # small value of a
        x0 = exp(logr)
    else    #large a
        haseta = true
        x0, fp = gamma_inc_inv_alarge(a, minpq, pcase)
    end

    t = 1
    x = x0
    n = 1
    logabsgam = logabsgamma(a)[1]
    # Newton-like higher order iteration with upper limit as 15.
    while t > 1.0e-15 && n < 15
        if !haseta
            dlnr = (1.0 - a)*log(x) + x + logabsgam
            if dlnr > log(floatmax(Float64)/1000.0)
                break
            else
                r = exp(dlnr)
            end
        else
            r = -(1/fp)*x
        end

        px, qx = gamma_inc(a, x, 0)

        ck1 = pcase ? -r*(px - minpq) : r*(qx - minpq)
        if a > 0.05
            ck2 = (x - a + 1.0)/(2.0*x)

            # This check is not in the invincgam subroutine from IncgamFI but it probably
            # should be since the routine fails to compute a finite value for very small p
            if !isfinite(ck2)
                break
            end

            if a > 0.1
                ck3 = (@horner(x, @horner(a, 1, -3, 2), @horner(a, 4, -4), 2))/(6*x^2)

                # This check is not in the invincgam subroutine from IncgamFI but it probably
                # should be since the routine fails to compute a finite value for very small p
                if !isfinite(ck3)
                    break
                end
                Δx = ck1 * @horner(ck1, 1.0, ck2, ck3)
            else
                Δx = ck1 * @horner(ck1, 1.0, ck2)
            end
        else
            Δx = ck1
        end

        x0 = x + Δx
        t = abs(Δx/x0)
        n += 1
        x = x0
    end
    return x
end

function __gamma_inc_inv(a::T, minpq::T, pcase::Bool) where {T<:Union{Float16,Float32}}
    return T(__gamma_inc_inv(Float64(a), Float64(minpq), pcase))
end

# like promote(x,y), but don't complexify real values
promotereal(x::Real, y::Real) = promote(x, y)
promotereal(x::Real, y::Number) = let (u, v) = promote(x, y); real(u), v; end
promotereal(x::Number, y::Real) = let (u, v) = promote(x, y); u, real(v); end
promotereal(x, y) = promote(x,y)

"""
    gamma(a,x)

Returns the upper incomplete gamma function
```math
\\Gamma(a,x) = \\int_x^\\infty t^{a-1} e^{-t} dt \\,
```
supporting arbitrary real or complex `a` and `x`.

(The ordinary gamma function [`gamma(x)`](@ref) corresponds to ``\\Gamma(a) = \\Gamma(a,0)``.
See also the [`gamma_inc`](@ref) function to compute both the upper and lower
(``\\gamma(a,x)``) incomplete gamma functions scaled by ``\\Gamma(a)``.

External links: [DLMF](https://dlmf.nist.gov/8.2.2), [Wikipedia](https://en.wikipedia.org/wiki/Incomplete_gamma_function)
"""
gamma(a::Number,x::Number) = _gamma(promotereal(float(a), float(x))...)
gamma(a::Integer,x::Number) = _gamma(a, float(x))

function _gamma(a::Number, x::Number)
    if a isa Real && x isa Real && !isfinite(a*x)
        if isinf(x) && isfinite(a)
            if x > 0 # == +Inf
                return one(a)*zero(x)
            elseif a > 0 && isinteger(a) # x == -Inf
                return one(a)*x # -Inf
            end
        elseif isinf(a) && isfinite(x)
            if a > 0
                if x ≥ 0
                    return a*one(x) # +Inf
                else
                    throw(DomainError((a, x), "gamma will only return a complex result if called with a complex argument"))
                end
            elseif a < 0
                return zero(a)*one(x)
            end
        end
    end
    return iszero(x) ? gamma(one(x)*a) : x^a * expint(1 - a, x)
end

_gamma(a::Integer, x::BigFloat) = _gamma_big(a, x)
_gamma(a::BigInt, x::Real) = _gamma_big(a, x)
_gamma(a::BigInt, x::BigFloat) = _gamma_big(a, x)
_gamma(a::BigFloat,x::BigFloat) = _gamma_big(a, x)
function _gamma_big(a::Real,x::Real)
    if x < 0
        # MPFR returns NaN in this case
        if isinteger(a) && a > 0
            return invoke(_gamma, Tuple{Number,Number}, a, BigFloat(x))
        else
            throw(DomainError((a, x), "gamma will only return a complex result if called with a complex argument"))
        end
    end
    z = BigFloat()
    ccall((:mpfr_gamma_inc, :libmpfr), Int32, (Ref{BigFloat}, Ref{BigFloat}, Ref{BigFloat}, Int32), z, a, x, ROUNDING_MODE[])
    return z
end

"""
    loggamma(a,x)

Returns the log of the upper incomplete gamma function [`gamma(a,x)`](@ref):
```math
\\log \\Gamma(a,x) = \\log \\int_x^\\infty t^{a-1} e^{-t} dt \\,
```
supporting arbitrary real or complex `a` and `x`.

If `a` and/or `x` is complex, then `exp(loggamma(a,x))` matches `gamma(a,x)` (up to floating-point error),
but `loggamma(a,x)` may differ from `log(gamma(a,x))` by an integer multiple of ``2\\pi i``
(i.e. it may employ a different branch cut).

See also [`loggamma(x)`](@ref).
"""
loggamma(a::Number, x::Number) = _loggamma(promotereal(float(a), float(x))...)
loggamma(a::Integer, x::Number) = _loggamma(a, float(x))
function _loggamma(a::Number, x::Number)
    if a isa Real && x isa Real && !isfinite(a*x)
        if isinf(x) && isfinite(a)
            if x > 0 # == +Inf
                return -one(a)*x
            elseif x < 0
                throw(DomainError((a, x), "loggamma will only return a complex result if called with a complex argument"))
            end
        elseif isinf(a) && isfinite(x)
            if a > 0 && x ≥ 0
                return a*one(x) # +Inf
            elseif a < 0
                return a*one(x) # -Inf
            end
        end
    end
    # from gamma(a,x) = x^a * expintx(1-a, x) * exp(-x):
    iszero(x) && return loggamma(one(x)*a)
    if x isa Real && x < 0 && a isa Integer && isodd(a)
        # minus signs in expintx and x^a may cancel
        return a*log(-x) + log(-expintx(1-a, x)) - x
    end
    return a*log(x) + log(expintx(1 - a, x)) - x
end
