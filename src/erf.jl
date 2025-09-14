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

    _erfc(x::Float64) = ccall((:erfc, libopenlibm), Float64, (Float64,), x)
    _erfc(x::Float32) = ccall((:erfcf, libopenlibm), Float32, (Float32,), x)
    _erfc(x::Float16) = Float16(_erfc(Float32(x)))

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

External links:
[DLMF 7.2.1](https://dlmf.nist.gov/7.2.1),
[Wikipedia](https://en.wikipedia.org/wiki/Error_function).

See also:
[`erfc(x)`](@ref erfc), [`erfcx(x)`](@ref erfcx),
[`erfi(x)`](@ref erfi), [`dawson(x)`](@ref dawson),
[`erfinv(x)`](@ref erfinv), [`erfcinv(x)`](@ref erfcinv).

# Implementation by
- `Float32`: polynomial approximations of erf 
- `Float64`: Julia implementation of https://github.com/ARM-software/optimized-routines/blob/master/math/erf.c 
- `BigFloat`: C library for multiple-precision floating-point [MPFR](https://www.mpfr.org/)
"""

#    Fast erf implementation using a mix of
#    rational and polynomial approximations.
#    Highest measured error is 1.01 ULPs at 0x1.39956ac43382fp+0.  
function _erf(x::Float64)

    # # top 32 bits 
    ix::UInt32=reinterpret(UInt64,x)>>32
    # # top 32, without sign bit 
    ia::UInt32=ix & 0x7fffffff

    if (ia < 0x3feb0000)
    #  a = |x| < 0.84375.  

        x2 = x * x

        if (ia < 0x3fe00000)
        ## a < 0.5  - Use polynomial approximation.  

            # Minimax approximation of erf of the form x*P(x^2) approximately on the interval [0;0.5]
            PA=(0.12837916709551256, -0.3761263890318287, 0.11283791670896592, -0.026866170630075903, 0.005223977428649761, -0.0008548312229989974, 0.00012054654502151287, -1.4906315067498891e-5, 1.6126444349070824e-6, -1.3074259056679966e-7)
        
            r= fma(x,evalpoly(x2,PA),x) ## This fma is crucial for accuracy.  
            return r
        else
        ## 0.5 <= a < 0.84375 - Use rational approximation.  

            # Rational approximation on [0x1p-28, 0.84375] 
            NA=(0x1.06eba8214db68p-3, -0x1.4cd7d691cb913p-2, -0x1.d2a51dbd7194fp-6,-0x1.7a291236668e4p-8, -0x1.8ead6120016acp-16)
            DA=(1,0x1.97779cddadc09p-2, 0x1.0a54c5536cebap-4, 0x1.4d022c4d36b0fp-8,0x1.15dc9221c1a1p-13, -0x1.09c4342a2612p-18)

            P=evalpoly(x2,NA)
            Q=evalpoly(x2,DA)

            return fma(P / Q, x, x)
        end
    elseif (ia < 0x3ff40000)
    ## 0.84375 <= |x| < 1.25.  

        # Rational approximation on [0.84375, 1.25] 
        NB=( -0x1.359b8bef77538p-9, 0x1.a8d00ad92b34dp-2, -0x1.7d240fbb8c3f1p-2, 0x1.45fca805120e4p-2, -0x1.c63983d3e28ecp-4, 0x1.22a36599795ebp-5, -0x1.1bf380a96073fp-9 )
        DB=( 1, 0x1.b3e6618eee323p-4, 0x1.14af092eb6f33p-1, 0x1.2635cd99fe9a7p-4, 0x1.02660e763351fp-3, 0x1.bedc26b51dd1cp-7, 0x1.88b545735151dp-7 )

        C = 0x1.b0ac16p-1

        a = abs(x) - 1.0

        P=evalpoly(a,NB)
        Q=evalpoly(a,DB)

        r= C + P / Q
        return copysign(r,x)

    elseif (ia < 0x40000000)
    ## 1.25 <= |x| < 2.0.  
        a = abs(x)
        a = a - 1.25
        
        # erfc polynomial approximation
        # Generated using Sollya::remez(f(c*x+d), deg, [(a-d)/c;(b-d)/c], 1, 1e-16), [|D ...|] with deg=15 a=1.25 b=2 c=1 d=1.25 
        PC=( 0x1.3bcd133aa0ffcp-4, -0x1.e4652fadcb702p-3, 0x1.2ebf3dcca0446p-2, -0x1.571d01c62d66p-3, 0x1.93a9a8f5b3413p-8, 0x1.8281cbcc2cd52p-5, -0x1.5cffd86b4de16p-6, -0x1.db4ccf595053ep-9, 0x1.757cbf8684edap-8, -0x1.ce7dfd2a9e56ap-11, -0x1.99ee3bc5a3263p-11, 0x1.3c57cf9213f5fp-12, 0x1.60692996bf254p-14, -0x1.6e44cb7c1fa2ap-14, 0x1.9d4484ac482b2p-16, -0x1.578c9e375d37p-19)

        # Obtains erfc of |x|
        r=1.0-evalpoly(a,PC)
        return copysign(r,x)

    elseif (ia < 0x400a0000)
    ## 2 <= |x| < 3.25.  
        a = abs(x)
        a = fma(0.5, a, -1.0)

        # Generated using Sollya::remez(f(c*x+d), deg, [(a-d)/c;(b-d)/c], 1, 1e-16), [|D ...|] with deg=17 a=2 b=3.25 c=2 d=2 
        PD=( 0x1.328f5ec350e5p-8, -0x1.529b9e8cf8e99p-5, 0x1.529b9e8cd9e71p-3, -0x1.8b0ae3a023bf2p-2, 0x1.1a2c592599d82p-1, -0x1.ace732477e494p-2, -0x1.e1a06a27920ffp-6, 0x1.bae92a6d27af6p-2, -0x1.a15470fcf5ce7p-2, 0x1.bafe45d18e213p-6, 0x1.0d950680d199ap-2, -0x1.8c9481e8f22e3p-3, -0x1.158450ed5c899p-4, 0x1.c01f2973b44p-3, -0x1.73ed2827546a7p-3, 0x1.47733687d1ff7p-4, -0x1.2dec70d00b8e1p-6, 0x1.a947ab83cd4fp-10 )

        # Obtains erfc of |x|
        r=1.0-evalpoly(a,PD)
        return copysign(r,x)

    elseif (ia < 0x40100000)
    ## 3.25 <= |x| < 4.0.  
        a = abs(x)
        a = a - 3.25

        # Generated using Sollya::remez(f(c*x+d), deg, [(a-d)/c;(b-d)/c], 1, 1e-16), [|D ...|] with deg=13 a=3.25 b=4 c=1 d=3.25 
        PE=( 0x1.20c13035539e4p-18, -0x1.e9b5e8d16df7ep-16, 0x1.8de3cd4733bf9p-14, -0x1.9aa48beb8382fp-13, 0x1.2c7d713370a9fp-12, -0x1.490b12110b9e2p-12, 0x1.1459c5d989d23p-12, -0x1.64b28e9f1269p-13, 0x1.57c76d9d05cf8p-14, -0x1.bf271d9951cf8p-16, 0x1.db7ea4d4535c9p-19, 0x1.91c2e102d5e49p-20, -0x1.e9f0826c2149ep-21, 0x1.60eebaea236e1p-23 )

        r=1.0-evalpoly(a,PE)
        return copysign(r,x)

    elseif (ia < 0x4017a000)
    ## 4 <= |x| < 5.90625.  
        a = abs(x)
        a = fma(0.5, a, -2.0)

        # Generated using Sollya::remez(f(c*x+d), deg, [(a-d)/c;(b-d)/c], 1, 1e-16), [|D ...|] with deg=16 a=4 b=5.90625 c=2 d=4 
        PF=( 0x1.08ddd130d1fa6p-26, -0x1.10b146f59ff06p-22, 0x1.10b135328b7b2p-19, -0x1.6039988e7575fp-17, 0x1.497d365e19367p-15, -0x1.da48d9afac83ep-14, 0x1.1024c9b1fbb48p-12, -0x1.fc962e7066272p-12, 0x1.87297282d4651p-11, -0x1.f057b255f8c59p-11, 0x1.0228d0eee063p-10, -0x1.b1b21b84ec41cp-11, 0x1.1ead8ae9e1253p-11, -0x1.1e708fba37fccp-12, 0x1.9559363991edap-14, -0x1.68c827b783d9cp-16, 0x1.2ec4adeccf4a2p-19 )

        r=1.0-evalpoly(a,PF)
        return copysign(r,x)
    else
        if(isnan(x))
             return NaN
        end
        copysign(1.0,x)
    end
end

#    Fast erf implementation using 
#    polynomial approximations of erf and erfc.
#    Highest measured error is 1.12 ULPs at x = 1.232469  
function _erf(x::Float32)
    xabs=abs(x)

  if (xabs< 0.5)
    # range [0;0.5]
    # # erf approximation using erf(x)=x+x*P(x^2) with degree 6
    # Sollya::remez(erf(sqrt(x))/sqrt(x)-1,6,[1e-32;0.5],1,1e-32);

    p=(0.12837917f0, -0.37612638f0, 0.11283784f0, -0.026865287f0, 0.005218856f0, -0.0008394848f0)
    return copysign(fma(evalpoly(x^2,p),x,x),x)

  elseif(xabs<1.25)
    # range [0.5;1.25]
    # # erfc approximation with degree 11
    # Sollya::remez(erfc(x+0.5),11,[0;1.25-0.5],1,1e-32);
    
    p=(0.47950011f0, -0.8787826f0, 0.4393913f0, 0.14646415f0, -0.18308467f0, -0.007286422f0, 0.04987047f0, -0.0048868246f0, -0.011067663f0, 0.003422347f0, 0.00073027064f0, -0.0003758171f0)
    return copysign(1f0-evalpoly(xabs-0.5f0,p),x)

  elseif(xabs<2.5)
    # range [1.25;2.5]
    # erfc approximation with degree 13
    # Sollya::remez(erfc(x+1.25),13,[1.25-1.25;2.5-1.25],1,1e-32);

    p=(0.077099875f0, -0.23652112f0, 0.2956514f0, -0.16753574f0, 0.006158593f0, 0.04718712f0, -0.021331023f0, -0.0035262543f0, 0.005461831f0, -0.00047858673f0, -0.0012763853f0, 0.00073386944f0, -0.00017831658f0, 1.7451624f-5)
    return copysign(1f0-evalpoly(xabs-1.25f0,p),x)

  elseif (xabs<4.0)
    # range [2.5;4.0]
    # erfc approximation with degree 13
    # Sollya::remez(erfc(x+2.5),13,[0;4-2.5],1,1e-32);

    p=(0.00040695202f0, -0.002178284f0, 0.0054457085f0, -0.008350053f0, 0.008622011f0, -0.006115167f0, 0.0027899458f0, -0.000519395f0, -0.00030461047f0, 0.00031068458f0, -0.00013866898f0, 3.6909692f-5, -5.682889f-6, 3.929763f-7)
    return copysign(1f0-evalpoly(xabs-2.5f0,p),x)
    
  else
    # range [4.0,inf)
    r=copysign(1f0,x)
    return r
  end
end

_erf(x::Float16)=Float16(_erf(Float32(x)))


function erf end
"""
    erf(x, y)

Accurate version of `erf(y) - erf(x)` (for real arguments only).
"""
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

External links: [DLMF 7.2.2](https://dlmf.nist.gov/7.2.2),
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

External links: [DLMF 7.2.3](https://dlmf.nist.gov/7.2.3),
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

External links: [DLMF 7.2.5](https://dlmf.nist.gov/7.2.5),
[Wikipedia](https://en.wikipedia.org/wiki/Dawson_function).

See also: [`erfi(x)`](@ref erfi).

# Implementation by
- `Float32`/`Float64`: C standard math library
    [libm](https://en.wikipedia.org/wiki/C_mathematical_functions#libm).
"""
dawson

@doc raw"""
    faddeeva(z)

Compute the Faddeeva function of complex `z`, defined by
``e^{-z^2} \operatorname{erfc}(-iz)``.
Note that this function, also named `w` (original Faddeeva package) or `wofz` (Scilab package),
is equivalent to ``\operatorname{erfcx}(-iz)``.
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
Using the rational approximants tabulated in [Blair (1976)](@cite blair_1976) combined with
Newton iterations for `BigFloat`.
"""
erfinv(x::Real) = _erfinv(float(x))

function _erfinv(x::Float64)
    a = abs(x)
    if a > 1.0
        throw(DomainError(a, "`abs(x)` cannot be greater than 1."))
    elseif a == 1.0
        return copysign(Inf, x)
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
    if a > 1f0
        throw(DomainError(a, "`abs(x)` cannot be greater than 1."))
    elseif a == 1f0
        return copysign(Inf32, x)
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

function _erfinv(x::Float16)
    a = abs(x)
    if a > Float16(1)
        throw(DomainError(a, "`abs(x)` cannot be greater than 1."))
    elseif a == Float16(1)
        return copysign(Inf16, x)
    else
        # Perform calculations with `Float32`
        x32 = Float32(x)
        a32 = Float32(a)
        if a32 <= 0.75f0
            # Simpler and more accurate alternative to Table 7 in Blair et al.
            # Ref: https://github.com/JuliaMath/SpecialFunctions.jl/pull/372#discussion_r1592832735
            t = muladd(-6.73815f1, x32, 1f0) / muladd(-4.18798f0, x32, 4.54263f0)
            y = copysign(muladd(0.88622695f0, x32, t), x32)
        elseif a32 <= 0.9375f0 # Table 26 in Blair et al.
            t = x32^2 - 0.87890625f0
            y = x32 * @horner(t, 0.10178_950f1,
                              -0.32827_601f1) /
                      @horner(t, 0.72455_99f0,
                              -0.33871_553f1,
                              0.1f1)
        else
            # Simpler alternative to Table 47 in Blair et al.
            # because of the reduced accuracy requirement
            # (it turns out that this branch only covers 128 values).
            # Note that the use of log(1-x) rather than log1p is intentional since it will be
            # slightly faster and 1-x is exact.
            # Ref: https://github.com/JuliaMath/SpecialFunctions.jl/pull/372#discussion_r1592710586
            t = sqrt(-log(1-a32))
            y = copysign(@horner(t, -0.429159f0, 1.04868f0), x32)
        end
        return Float16(y)
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
Using the rational approximants tabulated in [Blair (1976)](@cite blair_1976)
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

function _erfcinv(y::Float16)
    if y > Float16(0.0625)
        return erfinv(Float16(1) - y)
    elseif y <= Float16(0)
        if y == Float16(0)
            return Inf16
        end
        throw(DomainError(y, "`y` must be nonnegative."))
    else # Table 47 in Blair et al.
        t = 1.0f0 / sqrt(-log(Float32(y)))
        x = @horner(t, 0.98650_088f0,
                       0.92601_777f0) /
            (t * @horner(t, 0.98424_719f0,
                            0.10074_7432f0,
                            0.1f0))
        return Float16(x)
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
\operatorname{logerfc}(x) = \ln(\operatorname{erfc}(x))
\quad \text{for} \quad x \in \mathbb{R} \, .
```

This is the accurate version of ``\ln(\operatorname{erfc}(x))`` for large ``x``.

External links: [Wikipedia](https://en.wikipedia.org/wiki/Error_function).

See also: [`erfcx(x)`](@ref erfcx).

# Implementation
Based on the [`erfc(x)`](@ref erfc) and [`erfcx(x)`](@ref erfcx) functions.
"""
function logerfc(x::Real)
    if x > zero(x)
        return log(erfcx(x)) - x^2
    else
        return log(erfc(x))
    end
end

@doc raw"""
    logerfcx(x)

Compute the natural logarithm of the scaled complementary error function of ``x``, that is

```math
\operatorname{logerfcx}(x) = \ln(\operatorname{erfcx}(x))
\quad \text{for} \quad x \in \mathbb{R} \, .
```

This is the accurate version of ``\ln(\operatorname{erfcx}(x))`` for large and negative ``x``.

External links: [Wikipedia](https://en.wikipedia.org/wiki/Error_function).

See also: [`erfcx(x)`](@ref erfcx).

# Implementation
Based on the [`erfc(x)`](@ref erfc) and [`erfcx(x)`](@ref erfcx) functions.
"""
function logerfcx(x::Real)
    if x < zero(x)
        return log(erfc(x)) + x^2
    else
        return log(erfcx(x))
    end
end

@doc raw"""
    logerf(x, y)

Compute the natural logarithm of two-argument error function. This is an accurate version of
`log(erf(x, y))`, which works for large `x`, `y`.

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
