# Functions
Here the *Special Functions* are listed according to the structure of [NIST Digital Library of Mathematical Functions](https://dlmf.nist.gov/).

## Error Functions, Dawson’s and Fresnel Integrals

| Function | Description |
|:-------- |:----------- |
| [`erf(x)`](@ref SpecialFunctions.erf)         | [error function](https://en.wikipedia.org/wiki/Error_function) at ``x`` |
| [`erfc(x)`](@ref SpecialFunctions.erfc)       | complementary error function, i.e. the accurate version of ``1-\operatorname{erf}(x)`` for large ``x`` |
| [`erfcinv(x)`](@ref SpecialFunctions.erfcinv) | inverse function to [`erfc()`](@ref SpecialFunctions.erfc) |
| [`erfcx(x)`](@ref SpecialFunctions.erfcx)     | scaled complementary error function, i.e. accurate ``e^{x^2} \operatorname{erfc}(x)`` for large ``x`` |
| [`erfi(x)`](@ref SpecialFunctions.erfi)       | imaginary error function defined as ``-i \operatorname{erf}(ix)`` |
| [`erfinv(x)`](@ref SpecialFunctions.erfinv)   | inverse function to [`erf()`](@ref SpecialFunctions.erf) |
| [`dawson(x)`](@ref SpecialFunctions.dawson)   | scaled imaginary error function, a.k.a. Dawson function, i.e. accurate ``\frac{\sqrt{\pi}}{2} e^{-x^2} \operatorname{erfi}(x)`` for large ``x`` |


## Elliptic Integrals

| Function | Description |
|:-------- |:----------- |
| [`ellipk(m)`](@ref SpecialFunctions.ellipk)    | [complete elliptic integral of 1st kind](https://en.wikipedia.org/wiki/Elliptic_integral#Notational_variants) ``K(m)``  |


# WIP not grouped
| Function                                                      | Description                                                                                                                                                     |
|:------------------------------------------------------------- |:--------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| [`sinint(x)`](@ref SpecialFunctions.sinint)                   | [sine integral](https://en.wikipedia.org/wiki/Trigonometric_integral) at `x` |
| [`cosint(x)`](@ref SpecialFunctions.cosint)                   | [cosine integral](https://en.wikipedia.org/wiki/Trigonometric_integral) at `x` |
| [`digamma(x)`](@ref SpecialFunctions.digamma)                 | [digamma function](https://en.wikipedia.org/wiki/Digamma_function) (i.e. the derivative of `lgamma` at `x`)                                                     |
| [`invdigamma(x)`](@ref SpecialFunctions.invdigamma)   | [invdigamma function](http://bariskurt.com/calculating-the-inverse-of-digamma-function/) (i.e. inverse of `digamma` function at `x` using fixed-point iteration algorithm) |
| [`trigamma(x)`](@ref SpecialFunctions.trigamma)     | [trigamma function](https://en.wikipedia.org/wiki/Trigamma_function) (i.e the logarithmic second derivative of `gamma` at `x`) |
| [`polygamma(m,x)`](@ref SpecialFunctions.polygamma)  | [polygamma function](https://en.wikipedia.org/wiki/Polygamma_function) (i.e the (m+1)-th derivative of the `lgamma` function at `x`) |
| [`gamma_inc(a,x,IND)`](@ref SpecialFunctions.gamma_inc)  | [incomplete gamma function ratio P(a,x) and Q(a,x)](https://en.wikipedia.org/wiki/Incomplete_gamma_function) (i.e evaluates P(a,x) and Q(a,x)for accuracy specified by IND and returns tuple (p,q)) |
| [`beta_inc(a,b,x,y)`](@ref SpecialFunctions.beta_inc)  | [incomplete beta function ratio Ix(a,b) and Iy(a,b)](https://en.wikipedia.org/wiki/Beta_function#Incomplete_beta_function) (i.e evaluates Ix(a,b) and Iy(a,b) and returns tuple (p,q)) |
| [`gamma_inc_inv(a,p,q)`](@ref SpecialFunctions.gamma_inc_inv)  | [inverse of incomplete gamma function ratio P(a,x) and Q(a,x)](https://en.wikipedia.org/wiki/Incomplete_gamma_function) (i.e evaluates x given P(a,x)=p and Q(a,x)=q  |
| [`ellipe(x)`](@ref SpecialFunctions.ellipe)    | [complete elliptic integral of 2nd kind](https://en.wikipedia.org/wiki/Elliptic_integral) (i.e evaluates complete elliptic integral of 2nd kind at `x`)  |
| [`eta(x)`](@ref SpecialFunctions.eta)                         | [Dirichlet eta function](https://en.wikipedia.org/wiki/Dirichlet_eta_function) at `x`                                                                           |
| [`zeta(x)`](@ref SpecialFunctions.zeta)                       | [Riemann zeta function](https://en.wikipedia.org/wiki/Riemann_zeta_function) at `x`                                                                             |
| [`airyai(z)`](@ref SpecialFunctions.airyai)                   | [Airy Ai function](https://en.wikipedia.org/wiki/Airy_function) at `z`                                                                                          |
| [`airyaiprime(z)`](@ref SpecialFunctions.airyaiprime)         | derivative of the Airy Ai function at `z`                                                                                                                       |
| [`airybi(z)`](@ref SpecialFunctions.airybi)                   | [Airy Bi function](https://en.wikipedia.org/wiki/Airy_function) at `z`                                                                                          |
| [`airybiprime(z)`](@ref SpecialFunctions.airybiprime)         | derivative of the Airy Bi function at `z`                                                                                                                       |
| [`airyaix(z)`](@ref SpecialFunctions.airyaix), [`airyaiprimex(z)`](@ref SpecialFunctions.airyaiprimex), [`airybix(z)`](@ref SpecialFunctions.airybix), [`airybiprimex(z)`](@ref SpecialFunctions.airybiprimex) | scaled Airy Ai function and `k`th derivatives at `z`                                                                                                           |
| [`besselj(nu,z)`](@ref SpecialFunctions.besselj)              | [Bessel function](https://en.wikipedia.org/wiki/Bessel_function) of the first kind of order `nu` at `z`                                                         |
| [`besselj0(z)`](@ref SpecialFunctions.besselj0)               | `besselj(0,z)`                                                                                                                                                  |
| [`besselj1(z)`](@ref SpecialFunctions.besselj1)               | `besselj(1,z)`                                                                                                                                                  |
| [`besseljx(nu,z)`](@ref SpecialFunctions.besseljx)            | scaled Bessel function of the first kind of order `nu` at `z`                                                                                                   |
| [`bessely(nu,z)`](@ref SpecialFunctions.bessely)              | [Bessel function](https://en.wikipedia.org/wiki/Bessel_function) of the second kind of order `nu` at `z`                                                        |
| [`bessely0(z)`](@ref SpecialFunctions.bessely0)               | `bessely(0,z)`                                                                                                                                                  |
| [`bessely1(z)`](@ref SpecialFunctions.bessely1)               | `bessely(1,z)`                                                                                                                                                  |
| [`besselyx(nu,z)`](@ref SpecialFunctions.besselyx)            | scaled Bessel function of the second kind of order `nu` at `z`                                                                                                  |
| [`besselh(nu,k,z)`](@ref SpecialFunctions.besselh)            | [Bessel function](https://en.wikipedia.org/wiki/Bessel_function) of the third kind (a.k.a. Hankel function) of order `nu` at `z`; `k` must be either `1` or `2` |
| [`hankelh1(nu,z)`](@ref SpecialFunctions.hankelh1)            | `besselh(nu, 1, z)`                                                                                                                                             |
| [`hankelh1x(nu,z)`](@ref SpecialFunctions.hankelh1x)          | scaled `besselh(nu, 1, z)`                                                                                                                                      |
| [`hankelh2(nu,z)`](@ref SpecialFunctions.hankelh2)            | `besselh(nu, 2, z)`                                                                                                                                             |
| [`hankelh2x(nu,z)`](@ref SpecialFunctions.hankelh2x)          | scaled `besselh(nu, 2, z)`                                                                                                                                      |
| [`besseli(nu,z)`](@ref SpecialFunctions.besseli)              | modified [Bessel function](https://en.wikipedia.org/wiki/Bessel_function) of the first kind of order `nu` at `z`                                                |
| [`besselix(nu,z)`](@ref SpecialFunctions.besselix)            | scaled modified Bessel function of the first kind of order `nu` at `z`                                                                                          |
| [`besselk(nu,z)`](@ref SpecialFunctions.besselk)              | modified [Bessel function](https://en.wikipedia.org/wiki/Bessel_function) of the second kind of order `nu` at `z`                                               |
| [`besselkx(nu,z)`](@ref SpecialFunctions.besselkx)            | scaled modified Bessel function of the second kind of order `nu` at `z`                                                                                         |
| [`gamma(x)`](@ref SpecialFunctions.gamma)                                            | [gamma function](https://en.wikipedia.org/wiki/Gamma_function) at `x`                                                                                           |
| [`loggamma(x)`](@ref SpecialFunctions.loggamma)                                           | accurate `log(gamma(x))` for large `x`                 |
| [`logabsgamma(x)`](@ref SpecialFunctions.logabsgamma)                                           | accurate `log(abs(gamma(x)))` for large `x`                                                                                  |
| [`lgamma(x)`](@ref SpecialFunctions.lgamma)                                           | accurate `log(gamma(x))` for large `x`                                                                                                                          |
| [`logfactorial(x)`](@ref SpecialFunctions.logfactorial)                                            | accurate `log(factorial(x))` for large `x`; same as `lgamma(x+1)` for `x > 1`, zero otherwise                                                                   |
| [`beta(x,y)`](@ref SpecialFunctions.beta)                                           | [beta function](https://en.wikipedia.org/wiki/Beta_function) at `x,y`                                                                                           |
| [`logbeta(x,y)`](@ref SpecialFunctions.logbeta)                                          | accurate `log(beta(x,y))` for large `x` or `y`      |
| [`logabsbeta(x,y)`](@ref SpecialFunctions.logabsbeta)                                          | accurate `log(abs(beta(x,y)))` for large `x` or `y`     |
| [`logabsbinomial(x,y)`](@ref SpecialFunctions.logabsbinomial)                                          | accurate `log(abs(beta(x,y)))` for large `x` or `y`     |
