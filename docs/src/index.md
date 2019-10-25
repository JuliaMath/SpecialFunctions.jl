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

On Julia 0.7, [openspecfun](https://github.com/JuliaLang/openspecfun) is not build as
part of Julia.
Thus for Julia versions 0.7 and later, installing this package downloads openspecfun.
