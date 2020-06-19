# Special Functions

This package provides a comprehensive collection of special functions based on the
[OpenSpecFun](https://github.com/JuliaLang/openspecfun) and [OpenLibm](https://github.com/JuliaLang/openlibm)
libraries.

## Installation

The package is available for Julia versions 0.5 and up. To install it, run

```julia
Pkg.add("SpecialFunctions")
```

from the Julia REPL.

## Note

Prior to Julia 0.6, most of these functions were available in Julia's `Base` module.
Because of this, the symbols from this package are not exported on Julia 0.5
to avoid name conflicts.
In this case, the symbols will need to be explicitly imported or called
with the prefix `SpecialFunctions`.
This is not necessary for Julia versions 0.6 and later.

On Julia 0.7, [openspecfun](https://github.com/JuliaLang/openspecfun) is not built as
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

to your `.juliarc.jl` file.
