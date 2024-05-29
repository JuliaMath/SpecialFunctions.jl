using Base.Math: @horner

function sinint(x::Float64)
    t = x*x
    if t ≤ 36.0
        return x * @horner(t, 1.00000000000000000000E0,
                             -0.44663998931312457298E-1,
                              0.11209146443112369449E-2,
                             -0.13276124407928422367E-4,
                              0.85118014179823463879E-7,
                             -0.29989314303147656479E-9,
                              0.55401971660186204711E-12,
                             -0.42406353433133212926E-15) /
                   @horner(t, 1.00000000000000000000E0,
                              0.10891556624243098264E-1,
                              0.59334456769186835896E-4,
                              0.21231112954641805908E-6,
                              0.54747121846510390750E-9,
                              0.10378561511331814674E-11,
                              0.13754880327250272679E-14,
                              0.10223981202236205703E-17)
    elseif t ≤ 144.0
        invt = inv(t)
        return copysign(π/2, x) - cos(x) *
                                         @horner(invt, 0.99999999962173909991E0,
                                                       0.36451060338631902917E3,
                                                       0.44218548041288440874E5,
                                                       0.22467569405961151887E7,
                                                       0.49315316723035561922E8,
                                                       0.43186795279670283193E9,
                                                       0.11847992519956804350E10,
                                                       0.45573267593795103181E9) /
                                    (x * @horner(invt, 1.00000000000000000000E0,
                                                       0.36651060273229347594E3,
                                                       0.44927569814970692777E5,
                                                       0.23285354882204041700E7,
                                                       0.53117852017228262911E8,
                                                       0.50335310667241870372E9,
                                                       0.16575285015623175410E10,
                                                       0.11746532837038341076E10)) -
                           sin(x)*invt * @horner(invt, 0.99999999920484901956E0,
                                                       0.51385504875307321394E3,
                                                       0.92293483452013810811E5,
                                                       0.74071341863359841727E7,
                                                       0.28142356162841356551E9,
                                                       0.49280890357734623984E10,
                                                       0.35524762685554302472E11,
                                                       0.79194271662085049376E11,
                                                       0.17942522624413898907E11) /
                                         @horner(invt, 1.00000000000000000000E0,
                                                       0.51985504708814870209E3,
                                                       0.95292615508125947321E5,
                                                       0.79215459679762667578E7,
                                                       0.31977567790733781460E9,
                                                       0.62273134702439012114E10,
                                                       0.54570971054996441467E11,
                                                       0.18241750166645704670E12,
                                                       0.15407148148861454434E12)
    elseif t < Inf
        invt = inv(t)
        return copysign(π/2, x) - cos(x) / x * (1.0 -
                                                @horner(invt, 0.19999999999999978257E1,
                                                              0.22206119380434958727E4,
                                                              0.84749007623988236808E6,
                                                              0.13959267954823943232E9,
                                                              0.10197205463267975592E11,
                                                              0.30229865264524075951E12,
                                                              0.27504053804288471142E13,
                                                              0.21818989704686874983E13) /
                                                @horner(invt, 1.00000000000000000000E0,
                                                              0.11223059690217167788E4,
                                                              0.43685270974851313242E6,
                                                              0.74654702140658116258E8,
                                                              0.58580034751805687471E10,
                                                              0.20157980379272098841E12,
                                                              0.26229141857684496445E13,
                                                              0.87852907334918467516E13)*invt) -
                         sin(x) * invt * (1.0 - @horner(invt, 0.59999999999999993089E1,
                                                              0.96527746044997139158E4,
                                                              0.56077626996568834185E7,
                                                              0.15022667718927317198E10,
                                                              0.19644271064733088465E12,
                                                              0.12191368281163225043E14,
                                                              0.31924389898645609533E15,
                                                              0.25876053010027485934E16,
                                                              0.12754978896268878403E16) /
                                                @horner(invt, 1.00000000000000000000E0,
                                                              0.16287957674166143196E4,
                                                              0.96636303195787870963E6,
                                                              0.26839734750950667021E9,
                                                              0.37388510548029219241E11,
                                                              0.26028585666152144496E13,
                                                              0.85134283716950697226E14,
                                                              0.11304079361627952930E16,
                                                              0.42519841479489798424E16)*invt)
    elseif isnan(x)
        return NaN
    else
        return copysign(π/2, x)
    end
end

function cosint(x::Float64)
    t, r0, r1   = x*x, 0.616505485620716233797110404100, 3.384180422551186426397851146402
    r01, r02    = 0.6162109375, 0.29454812071623379711E-3
    r11, r12    = 3.3837890625, 0.39136005118642639785E-3
    if x < 0.0
        throw(DomainError(x, "`x` must be nonnegative."))
    elseif x ≤ 3.0
        return log(x/r0) + ((x - r01) - r02) * (x + r0) *
                                       @horner(t, -0.24607411378767540707E0,
                                                   0.72113492241301534559E-2,
                                                  -0.11867127836204767056E-3,
                                                   0.90542655466969866243E-6,
                                                  -0.34322242412444409037E-8,
                                                   0.51950683460656886834E-11) /
                                        @horner(t, 1.00000000000000000000E0,
                                                   0.12670095552700637845E-1,
                                                   0.78168450570724148921E-4,
                                                   0.29959200177005821677E-6,
                                                   0.73191677761328838216E-9,
                                                   0.94351174530907529061E-12)
    elseif x ≤ 6.0
        return log(x/r1) + ((x - r11) - r12) * (x + r1) *
                                       @horner(t, -0.15684781827145408780E0,
                                                   0.66253165609605468916E-2,
                                                  -0.12822297297864512864E-3,
                                                   0.12360964097729408891E-5,
                                                  -0.66450975112876224532E-8,
                                                   0.20326936466803159446E-10,
                                                  -0.33590883135343844613E-13,
                                                   0.23686934961435015119E-16) /
                                        @horner(t, 1.00000000000000000000E0,
                                                   0.96166044388828741188E-2,
                                                   0.45257514591257035006E-4,
                                                   0.13544922659627723233E-6,
                                                   0.27715365686570002081E-9,
                                                   0.37718676301688932926E-12,
                                                   0.27706844497155995398E-15)
    elseif x ≤ 12.0
        invt = inv(t)
        return sin(x) * @horner(invt, 0.99999999962173909991E0,
                                      0.36451060338631902917E3,
                                      0.44218548041288440874E5,
                                      0.22467569405961151887E7,
                                      0.49315316723035561922E8,
                                      0.43186795279670283193E9,
                                      0.11847992519956804350E10,
                                      0.45573267593795103181E9) /
                   (x * @horner(invt, 1.00000000000000000000E0,
                                      0.36651060273229347594E3,
                                      0.44927569814970692777E5,
                                      0.23285354882204041700E7,
                                      0.53117852017228262911E8,
                                      0.50335310667241870372E9,
                                      0.16575285015623175410E10,
                                      0.11746532837038341076E10)) -
        cos(x) * invt * @horner(invt, 0.99999999920484901956E0,
                                      0.51385504875307321394E3,
                                      0.92293483452013810811E5,
                                      0.74071341863359841727E7,
                                      0.28142356162841356551E9,
                                      0.49280890357734623984E10,
                                      0.35524762685554302472E11,
                                      0.79194271662085049376E11,
                                      0.17942522624413898907E11) /
                        @horner(invt, 1.00000000000000000000E0,
                                      0.51985504708814870209E3,
                                      0.95292615508125947321E5,
                                      0.79215459679762667578E7,
                                      0.31977567790733781460E9,
                                      0.62273134702439012114E10,
                                      0.54570971054996441467E11,
                                      0.18241750166645704670E12,
                                      0.15407148148861454434E12)
    elseif x < Inf
        invt = inv(t)
        return sin(x)/x * (1.0 - @horner(invt, 0.19999999999999978257E1,
                                               0.22206119380434958727E4,
                                               0.84749007623988236808E6,
                                               0.13959267954823943232E9,
                                               0.10197205463267975592E11,
                                               0.30229865264524075951E12,
                                               0.27504053804288471142E13,
                                               0.21818989704686874983E13) /
                                 @horner(invt, 1.00000000000000000000E0,
                                               0.11223059690217167788E4,
                                               0.43685270974851313242E6,
                                               0.74654702140658116258E8,
                                               0.58580034751805687471E10,
                                               0.20157980379272098841E12,
                                               0.26229141857684496445E13,
                                               0.87852907334918467516E13)*invt) -
            cos(x)*invt * (1.0 - @horner(invt, 0.59999999999999993089E1,
                                               0.96527746044997139158E4,
                                               0.56077626996568834185E7,
                                               0.15022667718927317198E10,
                                               0.19644271064733088465E12,
                                               0.12191368281163225043E14,
                                               0.31924389898645609533E15,
                                               0.25876053010027485934E16,
                                               0.12754978896268878403E16) /
                                 @horner(invt, 1.00000000000000000000E0,
                                               0.16287957674166143196E4,
                                               0.96636303195787870963E6,
                                               0.26839734750950667021E9,
                                               0.37388510548029219241E11,
                                               0.26028585666152144496E13,
                                               0.85134283716950697226E14,
                                               0.11304079361627952930E16,
                                               0.42519841479489798424E16)*invt)
    elseif isnan(x)
        return NaN
    else
        return 0.0
    end
end

for f in (:sinint, :cosint)
    @eval begin
        ($f)(x::Float32)        = Float32(($f)(Float64(x)))
        ($f)(x::Float16)        = Float16(($f)(Float64(x)))
        ($f)(x::Real)           = ($f)(float(x))
        ($f)(x::AbstractFloat)  = error("not implemented for ", typeof(x))
    end
end


@doc raw"""
    sinint(x)

Compute the sine integral function of ``x``, defined by

```math
\operatorname{Si}(x)
:= \int_0^x \frac{\sin t}{t} \, dt
\quad \text{for} \quad
x \in \mathbb{R} \,.
```

External links: [DLMF](https://dlmf.nist.gov/6.2.9),
[Wikipedia](https://en.wikipedia.org/wiki/Trigonometric_integral#Sine_integral).

See also: [`cosint(x)`](@ref SpecialFunctions.cosint).

# Implementation
Using the rational approximants tabulated in:
> A.J. MacLeod,
> "Rational approximations, software and test methods for sine and cosine integrals",
> Numer. Algor. 12, pp. 259--272 (1996).
> <https://doi.org/10.1007/BF02142806>,
> <https://link.springer.com/article/10.1007/BF02142806>.

Note: the second zero of ``\text{Ci}(x)`` has a typo that is fixed:
``r_1 = 3.38418 0422\mathbf{8} 51186 42639 78511 46402`` in the article, but is in fact:
``r_1 = 3.38418 0422\mathbf{5} 51186 42639 78511 46402``.
"""
sinint

@doc raw"""
    cosint(x)

Compute the cosine integral function of ``x``, defined by

```math
\operatorname{Ci}(x)
:= \gamma + \log x + \int_0^x \frac{\cos (t) - 1}{t} \, dt
\quad \text{for} \quad
x > 0 \,,
```
where ``\gamma`` is the Euler-Mascheroni constant.

External links: [DLMF](https://dlmf.nist.gov/6.2.13),
[Wikipedia](https://en.wikipedia.org/wiki/Trigonometric_integral#Cosine_integral).

See also: [`sinint(x)`](@ref SpecialFunctions.sinint).

# Implementation
Using the rational approximants tabulated in:
> A.J. MacLeod,
> "Rational approximations, software and test methods for sine and cosine integrals",
> Numer. Algor. 12, pp. 259--272 (1996).
> <https://doi.org/10.1007/BF02142806>,
> <https://link.springer.com/article/10.1007/BF02142806>.

Note: the second zero of ``\text{Ci}(x)`` has a typo that is fixed:
``r_1 = 3.38418 0422\mathbf{8} 51186 42639 78511 46402`` in the article, but is in fact:
``r_1 = 3.38418 0422\mathbf{5} 51186 42639 78511 46402``.
"""
cosint
