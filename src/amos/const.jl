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


"""
    GAMMALN_GLN :: Vector{Float64}

LNGAMMA(N), N=1,100

# Impl Ref
- [`openspecfun/amos/dgamln.f`](https://github.com/JuliaMath/openspecfun/blob/v0.5.6/amos/dgamln.f)
- [`scipy/scipy/special/_amos.c:dgamln_cf[22]`](https://github.com/scipy/scipy/blob/b882f1b7ebe55e534f29a8d68a54e4ecd30aeb1a/scipy/special/_amos.c#L345-L371)
"""
const GAMMALN_GLN = Float64[
    0.00000000000000000e+00, 0.00000000000000000e+00, 6.93147180559945309e-01, 1.79175946922805500e+00,     #   0
    3.17805383034794562e+00, 4.78749174278204599e+00, 6.57925121201010100e+00, 8.52516136106541430e+00,     #   4
    1.06046029027452502e+01, 1.28018274800814696e+01, 1.51044125730755153e+01, 1.75023078458738858e+01,     #   8
    1.99872144956618861e+01, 2.25521638531234229e+01, 2.51912211827386815e+01, 2.78992713838408916e+01,     #  12
    3.06718601060806728e+01, 3.35050734501368889e+01, 3.63954452080330536e+01, 3.93398841871994940e+01,     #  16
    4.23356164607534850e+01, 4.53801388984769080e+01, 4.84711813518352239e+01, 5.16066755677643736e+01,     #  20
    5.47847293981123192e+01, 5.80036052229805199e+01, 6.12617017610020020e+01, 6.45575386270063311e+01,     #  24
    6.78897431371815350e+01, 7.12570389671680090e+01, 7.46582363488301644e+01, 7.80922235533153106e+01,     #  28
    8.15579594561150372e+01, 8.50544670175815174e+01, 8.85808275421976788e+01, 9.21361756036870925e+01,     #  32
    9.57196945421432025e+01, 9.93306124547874269e+01, 1.02968198614513813e+02, 1.06631760260643459e+02,     #  36
    1.10320639714757395e+02, 1.14034211781461703e+02, 1.17771881399745072e+02, 1.21533081515438634e+02,     #  40
    1.25317271149356895e+02, 1.29123933639127215e+02, 1.32952575035616310e+02, 1.36802722637326368e+02,     #  44
    1.40673923648234259e+02, 1.44565743946344886e+02, 1.48477766951773032e+02, 1.52409592584497358e+02,     #  48
    1.56360836303078785e+02, 1.60331128216630907e+02, 1.64320112263195181e+02, 1.68327445448427652e+02,     #  52
    1.72352797139162802e+02, 1.76395848406997352e+02, 1.80456291417543771e+02, 1.84533828861449491e+02,     #  56
    1.88628173423671591e+02, 1.92739047287844902e+02, 1.96866181672889994e+02, 2.01009316399281527e+02,     #  60
    2.05168199482641199e+02, 2.09342586752536836e+02, 2.13532241494563261e+02, 2.17736934113954227e+02,     #  64
    2.21956441819130334e+02, 2.26190548323727593e+02, 2.30439043565776952e+02, 2.34701723442818268e+02,     #  68
    2.38978389561834323e+02, 2.43268849002982714e+02, 2.47572914096186884e+02, 2.51890402209723194e+02,     #  72
    2.56221135550009525e+02, 2.60564940971863209e+02, 2.64921649798552801e+02, 2.69291097651019823e+02,     #  76
    2.73673124285693704e+02, 2.78067573440366143e+02, 2.82474292687630396e+02, 2.86893133295426994e+02,     #  80
    2.91323950094270308e+02, 2.95766601350760624e+02, 3.00220948647014132e+02, 3.04686856765668715e+02,     #  84
    3.09164193580146922e+02, 3.13652829949879062e+02, 3.18152639620209327e+02, 3.22663499126726177e+02,     #  88
    3.27185287703775217e+02, 3.31717887196928473e+02, 3.36261181979198477e+02, 3.40815058870799018e+02,     #  92
    3.45379407062266854e+02, 3.49954118040770237e+02, 3.54539085519440809e+02, 3.59134205369575399e+02,     #  96
] # GAMMALN_GLN

"""
    GAMMALN_CF :: Vector{Float64}

COEFFICIENTS OF ASYMPTOTIC EXPANSION [22]

# Impl Ref
- [`openspecfun/amos/dgamln.f`](https://github.com/JuliaMath/openspecfun/blob/v0.5.6/amos/dgamln.f)
- [`scipy/scipy/special/_amos.c:dgamln_cf[22]`](https://github.com/scipy/scipy/blob/b882f1b7ebe55e534f29a8d68a54e4ecd30aeb1a/scipy/special/_amos.c#L373-L380)
"""
const GAMMALN_CF = Float64[
    8.33333333333333333e-02, -2.77777777777777778e-03, 7.93650793650793651e-04, -5.95238095238095238e-04,   #  0
    8.41750841750841751e-04, -1.91752691752691753e-03, 6.41025641025641026e-03, -2.95506535947712418e-02,   #  4
    1.79644372368830573e-01, -1.39243221690590112e+00, 1.34028640441683920e+01, -1.56848284626002017e+02,   #  8
    2.19310333333333333e+03, -3.61087712537249894e+04, 6.91472268851313067e+05, -1.52382215394074162e+07,   # 12
    3.82900751391414141e+08, -1.08822660357843911e+10, 3.47320283765002252e+11, -1.23696021422692745e+13,   # 16
    4.88788064793079335e+14, -2.13203339609193739e+16
] # GAMMALN_CF
