# This file contains code that was formerly a part of Julia. License is MIT: http://julialang.org/license

using Base.Math: @horner
using Base.MPFR: ROUNDING_MODE

for f in (:erf, :erfc)
    internalf = Symbol(:_, f)
    libopenlibmf = QuoteNode(f)
    libopenlibmf0 = QuoteNode(Symbol(f, :f))
    openspecfunf = QuoteNode(Symbol(:Faddeeva_, f))
    mpfrf = QuoteNode(Symbol(:mpfr_, f))
    @eval begin
        $f(x::Number) = $internalf(float(x))

        $internalf(x::Float64) = ccall(($libopenlibmf, libopenlibm), Float64, (Float64,), x)
        $internalf(x::Float32) = ccall(($libopenlibmf0, libopenlibm), Float32, (Float32,), x)
        $internalf(x::Float16) = Float16($internalf(Float32(x)))

        $internalf(z::Complex{Float64}) = Complex{Float64}(ccall(($openspecfunf, libopenspecfun), Complex{Float64}, (Complex{Float64}, Float64), z, zero(Float64)))
        $internalf(z::Complex{Float32}) = Complex{Float32}(ccall(($openspecfunf, libopenspecfun), Complex{Float64}, (Complex{Float64}, Float64), Complex{Float64}(z), Float64(eps(Float32))))
        $internalf(z::Complex{Float16}) = Complex{Float16}($internalf(Complex{Float32}(z)))

        function $internalf(x::BigFloat)
            z = BigFloat()
            ccall(($mpfrf, :libmpfr), Int32, (Ref{BigFloat}, Ref{BigFloat}, Int32), z, x, ROUNDING_MODE[])
            return z
        end
    end
end

for f in (:erfcx, :erfi, :dawson, :faddeeva)
    internalf = Symbol(:_, f)
    openspecfunfsym = Symbol(:Faddeeva_, f === :dawson ? :Dawson : f === :faddeeva ? :w : f)
    openspecfunfF64 = QuoteNode(Symbol(openspecfunfsym, :_re))
    openspecfunfCF64 = QuoteNode(openspecfunfsym)
    @eval begin
        $f(x::Number) = $internalf(float(x))

        $internalf(x::Float64) = ccall(($openspecfunfF64, libopenspecfun), Float64, (Float64,), x)
        $internalf(x::Float32) = Float32($internalf(Float64(x)))
        $internalf(x::Float16) = Float16($internalf(Float64(x)))

        $internalf(z::Complex{Float64}) = Complex{Float64}(ccall(($openspecfunfCF64, libopenspecfun), Complex{Float64}, (Complex{Float64}, Float64), z, zero(Float64)))
        $internalf(z::Complex{Float32}) = Complex{Float32}(ccall(($openspecfunfCF64, libopenspecfun), Complex{Float64}, (Complex{Float64}, Float64), Complex{Float64}(z), Float64(eps(Float32))))
        $internalf(z::Complex{Float16}) = Complex{Float16}($internalf(Complex{Float32}(z)))
    end
end

faddeeva(x::Real) = faddeeva(complex(x))

# MPFR has an open TODO item for this function
# until then, we use [DLMF 7.12.1](https://dlmf.nist.gov/7.12.1) for the tail
function _erfcx(x::BigFloat)
    if x <= (Clong == Int32 ? 0x1p15 : 0x1p30)
        # any larger gives internal overflow
        return exp(x^2)*erfc(x)
    elseif !isfinite(x)
        return 1/x
    else
        # asymptotic series
        # starts to diverge at iteration i = 2^30 or 2^60
        # final term will be < Γ(2*i+1)/(2^i * Γ(i+1)) / (2^(i+1))
        # so good to (lgamma(2*i+1) - lgamma(i+1))/log(2) - 2*i - 1
        #            ≈ 3.07e10 or 6.75e19 bits
        # which is larger than the memory of the respective machines
        ϵ = eps(BigFloat)/4
        v = 1/(2*x*x)
        k = 1
        s = w = -k*v
        while abs(w) > ϵ
            k += 2
            w *= -k*v
            s += w
        end
        return (1+s)/(x*sqrtπ)
    end
end

@doc raw"""
    erf(x)

Compute the error function of ``x``, defined by

```math
\operatorname{erf}(x) = \frac{2}{\sqrt{\pi}} \int_0^x \exp(-t^2) \; \mathrm{d}t
\quad \text{for} \quad x \in \mathbb{C} \, .
```

    erf(x, y)

Accurate version of `erf(y) - erf(x)` (for real arguments only).

External links: [DLMF](https://dlmf.nist.gov/7.2.1), [Wikipedia](https://en.wikipedia.org/wiki/Error_function).

See also: [`erfc(x)`](@ref erfc), [`erfcx(x)`](@ref erfcx),
[`erfi(x)`](@ref erfi), [`dawson(x)`](@ref dawson),
[`erfinv(x)`](@ref erfinv), [`erfcinv(x)`](@ref erfcinv).

# Implementation by
- `Float32`/`Float64`: C standard math library
    [libm](https://en.wikipedia.org/wiki/C_mathematical_functions#libm).
- `BigFloat`: C library for multiple-precision floating-point [MPFR](https://www.mpfr.org/)
"""
erf

function erf(x::Real, y::Real)
    if abs(x) ≤ 1/√2 && abs(y) ≤ 1/√2
        erf(y) - erf(x)
    elseif 0 ≤ x && 0 ≤ y
        erfc(x) - erfc(y)
    elseif x ≤ 0 && y ≤ 0
        erfc(-y) - erfc(-x)
    else
        erf(y) - erf(x)
    end
end

@doc raw"""
    erfc(x)

Compute the complementary error function of ``x``, defined by

```math
\operatorname{erfc}(x)
= 1 - \operatorname{erf}(x)
= \frac{2}{\sqrt{\pi}} \int_x^\infty \exp(-t^2) \; \mathrm{d}t
\quad \text{for} \quad x \in \mathbb{C} \, .
```

This is the accurate version of `1-erf(x)` for large ``x``.

External links: [DLMF](https://dlmf.nist.gov/7.2.2),
[Wikipedia](https://en.wikipedia.org/wiki/Error_function#Complementary_error_function).

See also: [`erf(x)`](@ref erf).

# Implementation by
- `Float32`/`Float64`: C standard math library
    [libm](https://en.wikipedia.org/wiki/C_mathematical_functions#libm).
- `BigFloat`: C library for multiple-precision floating-point [MPFR](https://www.mpfr.org/)
"""
erfc

@doc raw"""
    erfcx(x)

Compute the scaled complementary error function of ``x``, defined by

```math
\operatorname{erfcx}(x)
= e^{x^2} \operatorname{erfc}(x)
\quad \text{for} \quad x \in \mathbb{C} \, .
```

This is the accurate version of ``e^{x^2} \operatorname{erfc}(x)`` for large ``x``.
Note also that ``\operatorname{erfcx}(-ix)`` computes the Faddeeva function `w(x)`.

External links: [DLMF](https://dlmf.nist.gov/7.2.3),
[Wikipedia](https://en.wikipedia.org/wiki/Error_function#Complementary_error_function).

See also: [`erfc(x)`](@ref erfc).

# Implementation by
- `Float32`/`Float64`: C standard math library
    [libm](https://en.wikipedia.org/wiki/C_mathematical_functions#libm).
- `BigFloat`: MPFR has an open TODO item for this function until then, we use
    [DLMF 7.12.1](https://dlmf.nist.gov/7.12.1) for the tail.
"""
erfcx

@doc raw"""
    erfi(x)

Compute the imaginary error function of ``x``, defined by

```math
\operatorname{erfi}(x)
= -i \operatorname{erf}(ix)
\quad \text{for} \quad x \in \mathbb{C} \, .
```

External links:
[Wikipedia](https://en.wikipedia.org/wiki/Error_function#Imaginary_error_function).

See also: [`erf(x)`](@ref erf).

# Implementation by
- `Float32`/`Float64`: C standard math library
    [libm](https://en.wikipedia.org/wiki/C_mathematical_functions#libm).
"""
erfi

@doc raw"""
    dawson(x)

Compute the Dawson function (scaled imaginary error function) of ``x``, defined by

```math
\operatorname{dawson}(x)
= \frac{\sqrt{\pi}}{2} e^{-x^2} \operatorname{erfi}(x)
\quad \text{for} \quad x \in \mathbb{C} \, .
```

This is the accurate version of ``\frac{\sqrt{\pi}}{2} e^{-x^2} \operatorname{erfi}(x)``
for large ``x``.

External links: [DLMF](https://dlmf.nist.gov/7.2.5),
[Wikipedia](https://en.wikipedia.org/wiki/Dawson_function).

See also: [`erfi(x)`](@ref erfi).

# Implementation by
- `Float32`/`Float64`: C standard math library
    [libm](https://en.wikipedia.org/wiki/C_mathematical_functions#libm).
"""
dawson

"""
    faddeeva(z)

Compute the Faddeeva function of complex `z`, defined by
``e^{-z^2} \\operatorname{erfc}(-iz)``.
Note that this function, also named `w` (original Faddeeva package) or `wofz` (Scilab package),
is equivalent to``\\operatorname{erfcx}(-iz)``.
"""
faddeeva

@doc raw"""
    erfinv(x)

Compute the inverse error function of a real ``x``, that is

```math
\operatorname{erfinv}(x) = \operatorname{erf}^{-1}(x)
\quad \text{for} \quad x \in \mathbb{R} \, .
```

External links:
[Wikipedia](https://en.wikipedia.org/wiki/Error_function#Inverse_functions).

See also: [`erf(x)`](@ref erf).

# Implementation
Using the rational approximants tabulated in:
> J. M. Blair, C. A. Edwards, and J. H. Johnson,
> "Rational Chebyshev approximations for the inverse of the error function",
> Math. Comp. 30, pp. 827--830 (1976).
> <https://doi.org/10.1090/S0025-5718-1976-0421040-7>,
> <http://www.jstor.org/stable/2005402>
combined with Newton iterations for `BigFloat`.
"""
erfinv(x::Real) = _erfinv(float(x))

function _erfinv(x::Float64)
    a = abs(x)
    if a >= 1.0
        if x == 1.0
            return Inf
        elseif x == -1.0
            return -Inf
        end
        throw(DomainError(a, "`abs(x)` cannot be greater than 1."))
    elseif a <= 0.75 # Table 17 in Blair et al.
        t = x*x - 0.5625
        return x * @horner(t, 0.16030_49558_44066_229311e2,
                             -0.90784_95926_29603_26650e2,
                              0.18644_91486_16209_87391e3,
                             -0.16900_14273_46423_82420e3,
                              0.65454_66284_79448_7048e2,
                             -0.86421_30115_87247_794e1,
                              0.17605_87821_39059_0) /
                   @horner(t, 0.14780_64707_15138_316110e2,
                             -0.91374_16702_42603_13936e2,
                              0.21015_79048_62053_17714e3,
                             -0.22210_25412_18551_32366e3,
                              0.10760_45391_60551_23830e3,
                             -0.20601_07303_28265_443e2,
                              0.1e1)
    elseif a <= 0.9375 # Table 37 in Blair et al.
        t = x*x - 0.87890625
        return x * @horner(t, -0.15238_92634_40726_128e-1,
                               0.34445_56924_13612_5216,
                              -0.29344_39867_25424_78687e1,
                               0.11763_50570_52178_27302e2,
                              -0.22655_29282_31011_04193e2,
                               0.19121_33439_65803_30163e2,
                              -0.54789_27619_59831_8769e1,
                               0.23751_66890_24448) /
                   @horner(t, -0.10846_51696_02059_954e-1,
                               0.26106_28885_84307_8511,
                              -0.24068_31810_43937_57995e1,
                               0.10695_12997_33870_14469e2,
                              -0.23716_71552_15965_81025e2,
                               0.24640_15894_39172_84883e2,
                              -0.10014_37634_97830_70835e2,
                               0.1e1)
    else # Table 57 in Blair et al.
        t = inv(sqrt(-log1p(-a)))
        return @horner(t, 0.10501_31152_37334_38116e-3,
                          0.10532_61131_42333_38164_25e-1,
                          0.26987_80273_62432_83544_516,
                          0.23268_69578_89196_90806_414e1,
                          0.71678_54794_91079_96810_001e1,
                          0.85475_61182_21678_27825_185e1,
                          0.68738_08807_35438_39802_913e1,
                          0.36270_02483_09587_08930_02e1,
                          0.88606_27392_96515_46814_9) /
              (copysign(t, x) *
               @horner(t, 0.10501_26668_70303_37690e-3,
                          0.10532_86230_09333_27531_11e-1,
                          0.27019_86237_37515_54845_553,
                          0.23501_43639_79702_53259_123e1,
                          0.76078_02878_58012_77064_351e1,
                          0.11181_58610_40569_07827_3451e2,
                          0.11948_78791_84353_96667_8438e2,
                          0.81922_40974_72699_07893_913e1,
                          0.40993_87907_63680_15361_45e1,
                          0.1e1))
    end
end

function _erfinv(x::Float32)
    a = abs(x)
    if a >= 1.0f0
        if x == 1.0f0
            return Inf32
        elseif x == -1.0f0
            return -Inf32
        end
        throw(DomainError(a, "`abs(x)` cannot be greater than 1."))
    elseif a <= 0.75f0 # Table 10 in Blair et al.
        t = x*x - 0.5625f0
        return x * @horner(t, -0.13095_99674_22f2,
                               0.26785_22576_0f2,
                              -0.92890_57365f1) /
                   @horner(t, -0.12074_94262_97f2,
                               0.30960_61452_9f2,
                              -0.17149_97799_1f2,
                               0.1f1)
    elseif a <= 0.9375f0 # Table 29 in Blair et al.
        t = x*x - 0.87890625f0
        return x * @horner(t, -0.12402_56522_1f0,
                               0.10688_05957_4f1,
                              -0.19594_55607_8f1,
                               0.42305_81357f0) /
                   @horner(t, -0.88276_97997f-1,
                               0.89007_43359f0,
                              -0.21757_03119_6f1,
                               0.1f1)
    else # Table 50 in Blair et al.
        t = inv(sqrt(-log1p(-a)))
        return @horner(t, 0.15504_70003_116f0,
                          0.13827_19649_631f1,
                          0.69096_93488_87f0,
                         -0.11280_81391_617f1,
                          0.68054_42468_25f0,
                         -0.16444_15679_1f0) /
              (copysign(t, x) *
               @horner(t, 0.15502_48498_22f0,
                          0.13852_28141_995f1,
                          0.1f1))
    end
end

function _erfinv(y::BigFloat)
    xfloat = erfinv(Float64(y))
    if isfinite(xfloat)
        x = BigFloat(xfloat)
    else
        # Float64 overflowed, use asymptotic estimate instead
        # from erfc(x) ≈ exp(-x²)/x√π ≈ y  ⟹  -log(yπ) ≈ x² + log(x) ≈ x²
        x = copysign(sqrt(-log((1-abs(y))*sqrtπ)), y)
        isfinite(x) || return x
    end
    sqrtπhalf = sqrtπ * big(0.5)
    tol = 2eps(abs(x))
    while true # Newton iterations
        Δx = sqrtπhalf * (erf(x) - y) * exp(x^2)
        x -= Δx
        abs(Δx) < tol && break
    end
    return x
end

@doc raw"""
    erfcinv(x)

Compute the inverse error complementary function of a real ``x``, that is

```math
\operatorname{erfcinv}(x) = \operatorname{erfc}^{-1}(x)
\quad \text{for} \quad x \in \mathbb{R} \, .
```

External links:
[Wikipedia](https://en.wikipedia.org/wiki/Error_function#Inverse_functions).

See also: [`erfc(x)`](@ref erfc).

# Implementation
Using the rational approximants tabulated in:
> J. M. Blair, C. A. Edwards, and J. H. Johnson,
> "Rational Chebyshev approximations for the inverse of the error function",
> Math. Comp. 30, pp. 827--830 (1976).
> <https://doi.org/10.1090/S0025-5718-1976-0421040-7>,
> <http://www.jstor.org/stable/2005402>
combined with Newton iterations for `BigFloat`.
"""
erfcinv(x::Real) = _erfcinv(float(x))

function _erfcinv(y::Float64)
    if y > 0.0625
        return erfinv(1.0 - y)
    elseif y <= 0.0
        if y == 0.0
            return Inf
        end
        throw(DomainError(y, "`y` must be nonnegative."))
    elseif y >= 1e-100 # Table 57 in Blair et al.
        t = 1.0 / sqrt(-log(y))
        return @horner(t, 0.10501_31152_37334_38116e-3,
                          0.10532_61131_42333_38164_25e-1,
                          0.26987_80273_62432_83544_516,
                          0.23268_69578_89196_90806_414e1,
                          0.71678_54794_91079_96810_001e1,
                          0.85475_61182_21678_27825_185e1,
                          0.68738_08807_35438_39802_913e1,
                          0.36270_02483_09587_08930_02e1,
                          0.88606_27392_96515_46814_9) /
              (t *
               @horner(t, 0.10501_26668_70303_37690e-3,
                          0.10532_86230_09333_27531_11e-1,
                          0.27019_86237_37515_54845_553,
                          0.23501_43639_79702_53259_123e1,
                          0.76078_02878_58012_77064_351e1,
                          0.11181_58610_40569_07827_3451e2,
                          0.11948_78791_84353_96667_8438e2,
                          0.81922_40974_72699_07893_913e1,
                          0.40993_87907_63680_15361_45e1,
                          0.1e1))
    else # Table 80 in Blair et al.
        t = 1.0 / sqrt(-log(y))
        return @horner(t, 0.34654_29858_80863_50177e-9,
                          0.25084_67920_24075_70520_55e-6,
                          0.47378_13196_37286_02986_534e-4,
                          0.31312_60375_97786_96408_3388e-2,
                          0.77948_76454_41435_36994_854e-1,
                          0.70045_68123_35816_43868_271e0,
                          0.18710_42034_21679_31668_683e1,
                          0.71452_54774_31351_45428_3e0) /
          (t * @horner(t, 0.34654_29567_31595_11156e-9,
                          0.25084_69079_75880_27114_87e-6,
                          0.47379_53129_59749_13536_339e-4,
                          0.31320_63536_46177_68848_0813e-2,
                          0.78073_48906_27648_97214_733e-1,
                          0.70715_04479_95337_58619_993e0,
                          0.19998_51543_49112_15105_214e1,
                          0.15072_90269_27316_80008_56e1,
                          0.1e1))
    end
end

function _erfcinv(y::Float32)
    if y > 0.0625f0
        return erfinv(1.0f0 - y)
    elseif y <= 0.0f0
        if y == 0.0f0
            return Inf32
        end
        throw(DomainError(y, "`y` must be nonnegative."))
    else # Table 50 in Blair et al.
        t = 1.0f0 / sqrt(-log(y))
        return @horner(t, 0.15504_70003_116f0,
                          0.13827_19649_631f1,
                          0.69096_93488_87f0,
                         -0.11280_81391_617f1,
                          0.68054_42468_25f0,
                         -0.16444_15679_1f0) /
        (t * @horner(t, 0.15502_48498_22f0,
                        0.13852_28141_995f1,
                        0.1f1))
    end
end

function _erfcinv(y::BigFloat)
    yfloat = Float64(y)
    xfloat = erfcinv(yfloat)
    if isfinite(xfloat)
        x = BigFloat(xfloat)
    else
        # Float64 overflowed, use asymptotic estimate instead
        # from erfc(x) ≈ exp(-x²)/x√π ≈ y  ⟹  -log(yπ) ≈ x² + log(x) ≈ x²
        if yfloat < 1
            x = sqrt(-log(y*sqrtπ))
        else # y must be close to 2
            x = -sqrt(-log((2-y)*sqrtπ))
        end
        # TODO: Newton convergence is slow near y=0 singularity; accelerate?
        isfinite(x) || return x
    end
    sqrtπhalf = sqrtπ * big(0.5)
    tol = 2eps(abs(x))
    while true # Newton iterations
        Δx = sqrtπhalf * (erfc(x) - y) * exp(x^2)
        x += Δx
        abs(Δx) < tol && break
    end
    return x
end

@doc raw"""
    logerfc(x)

Compute the natural logarithm of the complementary error function of ``x``, that is

```math
\operatorname{logerfc}(x) = \operatorname{ln}(\operatorname{erfc}(x))
\quad \text{for} \quad x \in \mathbb{R} \, .
```

This is the accurate version of ``\operatorname{ln}(\operatorname{erfc}(x))`` for large ``x``.

External links: [Wikipedia](https://en.wikipedia.org/wiki/Error_function).

See also: [`erfcx(x)`](@ref erfcx).

# Implementation
Based on the [`erfc(x)`](@ref erfc) and [`erfcx(x)`](@ref erfcx) functions.
Currently only implemented for `Float32`, `Float64`, and `BigFloat`.
"""
logerfc(x::Real) = _logerfc(float(x))

function _logerfc(x::Union{Float32, Float64, BigFloat})
    # Don't include Float16 in the Union, otherwise logerfc would currently work for x <= 0.0, but not x > 0.0
    if x > 0.0
        return log(erfcx(x)) - x^2
    else
        return log(erfc(x))
    end
end

@doc raw"""
    logerfcx(x)

Compute the natural logarithm of the scaled complementary error function of ``x``, that is

```math
\operatorname{logerfcx}(x) = \operatorname{ln}(\operatorname{erfcx}(x))
\quad \text{for} \quad x \in \mathbb{R} \, .
```

This is the accurate version of ``\operatorname{ln}(\operatorname{erfcx}(x))`` for large and negative ``x``.

External links: [Wikipedia](https://en.wikipedia.org/wiki/Error_function).

See also: [`erfcx(x)`](@ref erfcx).

# Implementation
Based on the [`erfc(x)`](@ref erfc) and [`erfcx(x)`](@ref erfcx) functions.
Currently only implemented for `Float32`, `Float64`, and `BigFloat`.
"""
logerfcx(x::Real) = _logerfcx(float(x))

function _logerfcx(x::Union{Float32, Float64, BigFloat})
    # Don't include Float16 in the Union, otherwise logerfc would currently work for x <= 0.0, but not x > 0.0
    if x < 0.0
        return log(erfc(x)) + x^2
    else
        return log(erfcx(x))
    end
end

@doc raw"""
    logerf(x, y)

Compute the natural logarithm of two-argument error function. This is an accurate version of
 `log(erf(x, y))`, which works for large `x, y`.

External links: [Wikipedia](https://en.wikipedia.org/wiki/Error_function).

See also: [`erf(x,y)`](@ref erf).
"""
function logerf(a::Real, b::Real)
    if abs(a) ≤ invsqrt2 && abs(b) ≤ invsqrt2
        return log(erf(a, b))
    elseif b > a > 0
        return logerfc(a) + LogExpFunctions.log1mexp(logerfc(b) - logerfc(a))
    elseif a < b < 0
        return logerfc(-b) + LogExpFunctions.log1mexp(logerfc(-a) - logerfc(-b))
    else
        return log(erf(a, b))
    end
end
