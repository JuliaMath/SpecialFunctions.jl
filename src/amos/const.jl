# SPDX-License-Identifier: BSD-3-Clause OR MIT

# TODO: Use NamedTuple?
"""
    I1_MACH :: Int32

integer machine dependent constants.

`I1_MACH` can be used to obtain machine-dependent parameters
for the local machine environment.

# fortran comments
```fortran
C  I/O unit numbers.
C    I1MACH( 1) = the standard input unit.
C    I1MACH( 2) = the standard output unit.
C    I1MACH( 3) = the standard punch unit.
C    I1MACH( 4) = the standard error message unit.
C
C  Words.
C    I1MACH( 5) = the number of bits per integer storage unit.
C    I1MACH( 6) = the number of characters per integer storage unit.
C  Integers.
C    I1MACH( 7) = A, the base.
C    I1MACH( 8) = S, the number of base-A digits.
C    I1MACH( 9) = A**S - 1, the largest magnitude.
C    I1MACH(10) = B, the base.
C
C  Single-Precision
C    I1MACH(11) = T, the number of base-B digits.
C    I1MACH(12) = EMIN, the smallest exponent E.
C    I1MACH(13) = EMAX, the largest exponent E.
C  Double-Precision
C    I1MACH(14) = T, the number of base-B digits.
C    I1MACH(15) = EMIN, the smallest exponent E.
C    I1MACH(16) = EMAX, the largest exponent E.
```

# Impl Ref
- [`openspecfun/amos/i1mach.f`](https://github.com/JuliaMath/openspecfun/blob/v0.5.6/amos/i1mach.f)
- [`scipy/scipy/special/_amos.c`](https://github.com/scipy/scipy/blob/b882f1b7ebe55e534f29a8d68a54e4ecd30aeb1a/scipy/special/_amos.c#L134-L151)
"""
const I1_MACH = Int32[
# I/O unit numbers.
    5,           # [ 1] standard input
    6,           # [ 2] standard output
    7,           # [ 3] standard punch
    0,           # [ 4] standard error

# Words.
    32,          # [ 5] bits per integer
    4,           # [ 6] sizeof(Int32)
# Integers.
    2,           # [ 7] base for integers
    31,          # [ 8] digits of integer base
    2147483647,  # [ 9] typemax(Int32)

# Floating-Point Numbers.
    2,           # [10] FLT_RADIX;
    # Single-Precision :: Float32
    24,          # [11] FLT_MANT_DIG;
    -126,        # [12] FLT_MIN_EXP;
    128,         # [13] FLT_MAX_EXP;
    # Double-Precision :: Float64
    53,          # [14] DBL_MANT_DIG;
    -1021,       # [15] DBL_MIN_EXP;
    1024         # [16] DBL_MAX_EXP;
] # I1_MACH

"""
    D1_MACH :: Float64

double precision machine dependent constants.

`D1_MACH` can be used to obtain machine-dependent parameters
for the local machine environment.

# fortran comments
```fortran
C   D1MACH( 1) = B**(EMIN-1), 
C   D1MACH( 2) = B**EMAX*(1 - B**(-T)), the largest magnitude.
C   D1MACH( 3) = B**(-T), the smallest relative spacing.
C   D1MACH( 4) = B**(1-T), the largest relative spacing.
C   D1MACH( 5) = LOG10(B)
C
C   I1MACH(10) = B, the base.
C   I1MACH(14) = T, the number of base-B digits.
C   I1MACH(15) = EMIN, the smallest exponent E.
C   I1MACH(16) = EMAX, the largest exponent E.
```

# Impl Ref
- [`openspecfun/amos/d1mach.f`](https://github.com/JuliaMath/openspecfun/blob/v0.5.6/amos/d1mach.f)
- [`scipy/scipy/special/_amos.c`](https://github.com/scipy/scipy/blob/b882f1b7ebe55e534f29a8d68a54e4ecd30aeb1a/scipy/special/_amos.c#L126-L132)
"""
const D1_MACH = Float64[
# [1] smallest positive magnitude： np.finfo(np.float64).tiny
    2.2250738585072014e-308,
# [2] largest magnitude:  np.finfo(np.float64).max
    1.7976931348623157e+308,
# [3] smallest relative spacing:  eps(Float64)/2
    1.1102230246251565e-16,
# [4] largest relative spacing:  eps(Float64)
    2.220446049250313e-16,
# [5] log10(2)
    0.3010299956639812,
] # D1_MACH
