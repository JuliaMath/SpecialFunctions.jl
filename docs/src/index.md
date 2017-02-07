# Special Functions

This package provides a comprehensive collection of special functions based on the
[OpenSpecFun](https://github.com/JuliaLang/openspecfun) and [OpenLibm](https://github.com/JuliaLang/openlibm)
libraries.

| Function                                                      | Description                                                                                                                                                     |
|:------------------------------------------------------------- |:--------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| [`erf(x)`](@ref)                                              | [error function](https://en.wikipedia.org/wiki/Error_function) at `x`                                                                                           |
| [`erfc(x)`](@ref)                                             | complementary error function, i.e. the accurate version of `1-erf(x)` for large `x`                                                                             |
| [`erfinv(x)`](@ref)                                           | inverse function to [`erf()`](@ref)                                                                                                                             |
| `erfcinv(x)`                                                  | inverse function to [`erfc()`](@ref)                                                                                                                            |
| [`erfi(x)`](@ref)                                             | imaginary error function defined as `-im * erf(x * im)`, where [`im`](@ref) is the imaginary unit                                                               |
| [`erfcx(x)`](@ref)                                            | scaled complementary error function, i.e. accurate `exp(x^2) * erfc(x)` for large `x`                                                                           |
| [`dawson(x)`](@ref)                                           | scaled imaginary error function, a.k.a. Dawson function, i.e. accurate `exp(-x^2) * erfi(x) * sqrt(pi) / 2` for large `x`                                       |
| [`digamma(x)`](@ref)                                          | [digamma function](https://en.wikipedia.org/wiki/Digamma_function) (i.e. the derivative of [`lgamma()`](@ref)) at `x`                                           |
| [`eta(x)`](@ref)                                              | [Dirichlet eta function](https://en.wikipedia.org/wiki/Dirichlet_eta_function) at `x`                                                                           |
| [`zeta(x)`](@ref)                                             | [Riemann zeta function](https://en.wikipedia.org/wiki/Riemann_zeta_function) at `x`                                                                             |
| [`airyai(z)`](@ref)                                           | [Airy Ai function](https://en.wikipedia.org/wiki/Airy_function) at `z`                                                                                          |
| [`airyaiprime(z)`](@ref)                                      | derivative of the Airy Ai function at `z`                                                                                                                       |
| [`airybi(z)`](@ref)                                           | [Airy Bi function](https://en.wikipedia.org/wiki/Airy_function) at `z`                                                                                          |
| [`airybiprime(z)`](@ref)                                      | derivative of the Airy Bi function at `z`                                                                                                                       |
| [`airyaix(z)`](@ref), [`airyaiprimex(z)`](@ref), [`airybix(z)`](@ref), [`airybiprimex(z)`](@ref) | scaled Airy AI function and `k` th derivatives at `z`                                                                                                           |
| [`besselj(nu,z)`](@ref)                                       | [Bessel function](https://en.wikipedia.org/wiki/Bessel_function) of the first kind of order `nu` at `z`                                                         |
| [`besselj0(z)`](@ref)                                         | `besselj(0,z)`                                                                                                                                                  |
| [`besselj1(z)`](@ref)                                         | `besselj(1,z)`                                                                                                                                                  |
| [`besseljx(nu,z)`](@ref)                                      | scaled Bessel function of the first kind of order `nu` at `z`                                                                                                   |
| [`bessely(nu,z)`](@ref)                                       | [Bessel function](https://en.wikipedia.org/wiki/Bessel_function) of the second kind of order `nu` at `z`                                                        |
| [`bessely0(z)`](@ref)                                         | `bessely(0,z)`                                                                                                                                                  |
| [`bessely1(z)`](@ref)                                         | `bessely(1,z)`                                                                                                                                                  |
| [`besselyx(nu,z)`](@ref)                                      | scaled Bessel function of the second kind of order `nu` at `z`                                                                                                  |
| [`besselh(nu,k,z)`](@ref)                                     | [Bessel function](https://en.wikipedia.org/wiki/Bessel_function) of the third kind (a.k.a. Hankel function) of order `nu` at `z`; `k` must be either `1` or `2` |
| [`hankelh1(nu,z)`](@ref)                                      | `besselh(nu, 1, z)`                                                                                                                                             |
| [`hankelh1x(nu,z)`](@ref)                                     | scaled `besselh(nu, 1, z)`                                                                                                                                      |
| [`hankelh2(nu,z)`](@ref)                                      | `besselh(nu, 2, z)`                                                                                                                                             |
| [`hankelh2x(nu,z)`](@ref)                                     | scaled `besselh(nu, 2, z)`                                                                                                                                      |
| [`besseli(nu,z)`](@ref)                                       | modified [Bessel function](https://en.wikipedia.org/wiki/Bessel_function) of the first kind of order `nu` at `z`                                                |
| [`besselix(nu,z)`](@ref)                                      | scaled modified Bessel function of the first kind of order `nu` at `z`                                                                                          |
| [`besselk(nu,z)`](@ref)                                       | modified [Bessel function](https://en.wikipedia.org/wiki/Bessel_function) of the second kind of order `nu` at `z`                                               |
| [`besselkx(nu,z)`](@ref)                                      | scaled modified Bessel function of the second kind of order `nu` at `z`                                                                                         |

## Installation

The package is available for Julia versions 0.5 and up. To install it, run

```julia
Pkg.add("SpecialFunctions")
```

from the Julia REPL.

## Note

Prior to Julia 0.6, these functions were available in Julia's Base module.
Because of this, the symbols from this package are not exported on Julia 0.5
to avoid name conflicts.
In this case, the symbols will need to be explicitly imported or called
with the prefix `SpecialFunctions`.
This is not necessary for Julia versions 0.6 and later.
