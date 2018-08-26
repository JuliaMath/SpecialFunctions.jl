# Special Functions

This package provides a comprehensive collection of special functions based on the
[OpenSpecFun](https://github.com/JuliaLang/openspecfun) and [OpenLibm](https://github.com/JuliaLang/openlibm)
libraries.

| Function                                                      | Description                                                                                                                                                     |
|:------------------------------------------------------------- |:--------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| [`erf(x)`](@ref SpecialFunctions.erf)                         | [error function](https://en.wikipedia.org/wiki/Error_function) at `x`                                                                                           |
| [`erfc(x)`](@ref SpecialFunctions.erfc)                       | complementary error function, i.e. the accurate version of `1-erf(x)` for large `x`                                                                             |
| [`erfinv(x)`](@ref SpecialFunctions.erfinv)                   | inverse function to [`erf()`](@ref SpecialFunctions.erf)                                                                                                        |
| [`erfcinv(x)`](@ref SpecialFunctions.erfcinv)                 | inverse function to [`erfc()`](@ref SpecialFunctions.erfc)                                                                                                       |
| [`erfi(x)`](@ref SpecialFunctions.erfi)                       | imaginary error function defined as `-im * erf(x * im)`, where `im` is the imaginary unit                                                                       |
| [`erfcx(x)`](@ref SpecialFunctions.erfcx)                     | scaled complementary error function, i.e. accurate `exp(x^2) * erfc(x)` for large `x`                                                                           |
| [`dawson(x)`](@ref SpecialFunctions.dawson)                   | scaled imaginary error function, a.k.a. Dawson function, i.e. accurate `exp(-x^2) * erfi(x) * sqrt(pi) / 2` for large `x`                                       |
| [`sinint(x)`](@ref SpecialFunctions.sinint)                   | [sine integral](https://en.wikipedia.org/wiki/Trigonometric_integral) at `x` |
| [`cosint(x)`](@ref SpecialFunctions.cosint)                   | [cosine integral](https://en.wikipedia.org/wiki/Trigonometric_integral) at `x` |
| [`digamma(x)`](@ref SpecialFunctions.digamma)                 | [digamma function](https://en.wikipedia.org/wiki/Digamma_function) (i.e. the derivative of `lgamma` at `x`)                                                     |
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
| [`lgamma(x)`](@ref SpecialFunctions.lgamma)                                           | accurate `log(gamma(x))` for large `x`                                                                                                                          |
| [`lfact(x)`](@ref SpecialFunctions.lfact)                                            | accurate `log(factorial(x))` for large `x`; same as `lgamma(x+1)` for `x > 1`, zero otherwise                                                                   |
| [`beta(x,y)`](@ref SpecialFunctions.beta)                                           | [beta function](https://en.wikipedia.org/wiki/Beta_function) at `x,y`                                                                                           |
| [`lbeta(x,y)`](@ref SpecialFunctions.lbeta)                                          | accurate `log(beta(x,y))` for large `x` or `y`      |

## Installation

The package is available for Julia versions 0.5 and up. To install it, run

```julia
Pkg.add("SpecialFunctions")
```

from the Julia REPL.

## Note

Prior to Julia 0.6, most of these functions were available in Julia's Base module.
Because of this, the symbols from this package are not exported on Julia 0.5
to avoid name conflicts.
In this case, the symbols will need to be explicitly imported or called
with the prefix `SpecialFunctions`.
This is not necessary for Julia versions 0.6 and later.

On Julia 0.7, [openspecfun](https://github.com/JuliaLang/openspecfun) is not build as
part of Julia.
Thus for Julia versions 0.7 and later, installing this package downloads openspecfun.
Binaries of openspecfun are available for macOS, Windows, and Linux (glibc >= 2.6).
Other systems will need to build the library from source.
You can force a build from source by setting an environment variable called
`JULIA_SPECIALFUNCTIONS_BUILD_SOURCE` equal to `true` before running `Pkg.build`.
This ensures that the library is built locally from source, even if binaries are
available.
Doing this requires a C compiler (Clang on macOS and FreeBSD, GCC elsewhere) and
gfortran.
If you always want to build this library from source, consider adding

```julia
ENV["JULIA_SPECIALFUNCTIONS_BUILD_SOURCE"] = "true"
```

to your .juliarc.jl file.
