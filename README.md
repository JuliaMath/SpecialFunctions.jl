# SpecialFunctions.jl

Special mathematical functions in Julia, include Bessel, Hankel, Airy, error, Dawson, exponential (or sine and cosine) integrals,
eta, zeta, digamma, inverse digamma, trigamma, and polygamma functions.
Most of these functions were formerly part of Base in early versions of Julia.

CI (Linux, macOS, FreeBSD, Windows):
[![CI](https://github.com/JuliaMath/SpecialFunctions.jl/workflows/CI/badge.svg)](https://github.com/JuliaMath/SpecialFunctions.jl/actions?query=workflow%3ACI)
[![codecov](https://codecov.io/gh/JuliaMath/SpecialFunctions.jl/branch/master/graph/badge.svg?token=qIKzX2I5ZK)](https://codecov.io/gh/JuliaMath/SpecialFunctions.jl)

Documentation:
[![Documentation](https://img.shields.io/badge/docs-stable-blue.svg)](https://specialfunctions.JuliaMath.org/stable)
[![Documentation](https://img.shields.io/badge/docs-dev-blue.svg)](https://specialfunctions.JuliaMath.org/dev)

Test status (most recent release):
[![PkgEval](https://juliaci.github.io/NanosoldierReports/pkgeval_badges/S/SpecialFunctions.svg)](https://juliaci.github.io/NanosoldierReports/pkgeval_badges/report.html)

## Upgrading from SpecialFunctions 1

SpecialFunctions 2 has only a [single breaking change](https://github.com/JuliaMath/SpecialFunctions.jl/pull/297):
The removal of the type piracy `Base.factorial(x::Number) = gamma(x + 1)`.
For most users this change will not break anything but for users of `factorial` it might.
If you want to upgrade from SpecialFunctions 1 to SpecialFunctions 2 we recommend:

- If your code does not use `factorial` then update the compat entry for SpecialFunctions to e.g. `"1.8.1, 2"`.
- If your code does use `factorial` then check for all occurrences of `factorial`:

  - If `factorial` is called on an `Integer`, keep `factorial`,
  - Otherwise replace `factorial(x)` with a call to `gamma(x + 1)`.

  Afterwards update the compat entry for SpecialFunctions and check that your package works with SpecialFunctions 2.

As the previous overload of `factorial` was type piratical ([added 4 years ago when code was moved out of Base](https://github.com/fredrikekre/SpecialFunctions.jl/blame/148574086f3da1d9f7e05d4eb538f91a73775d96/src/gamma.jl#L757-L758)), it is possible that you used it without a direct dependency on SpecialFunctions as long as SpecialFunctions was loaded.
The package ecosystem was analyzed and this only impacted a couple of packages. However, it is possible that private packages that depend on this may need updating, or stay with the older release of SpecialFunctions.jl.
