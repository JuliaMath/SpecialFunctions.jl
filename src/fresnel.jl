# Code ported from https://github.com/kiranshila/FresnelIntegrals.jl


"""
    fresnelc(z::Number)

Calculate the normalized Fresnel cosine integral
```math
C(z) = \\int_{0}^{z} \\cos{\\left(\\frac{\\pi t^2}{2}\\right)} \, \mathrm{d}t
```
for the number ``z``.
"""
function fresnelc(z::Number)
    x = (z * sqrtπ) / 2
    re_x, im_x = reim(x)
    a = (re_x - im_x) + (re_x - im_x) * im
    b = (re_x + im_x) + (im_x - re_x) * im
    re_erf_a, im_erf_a = reim(erf(a))
    re_erf_b, im_erf_b = reim(erf(b))
    re_y = (re_erf_a + im_erf_a + re_erf_b - im_erf_b) / 4
    im_y = (im_erf_a - re_erf_a + re_erf_b + im_erf_b) / 4
    y = re_y + im_y * im
    return y
end
function fresnelc(z::Real)
    x = (z * sqrtπ) / 2
    a = x + x * im
    re_erf_a, im_erf_a = reim(erf(a))
    y = (re_erf_a + im_erf_a) / 2
    return y
end

"""
    fresnels(z::Number)

Calculate the normalized Fresnel sine integral
```math
S(z) = \\int_{0}^{z} \\sin{\\left(\\frac{\\pi t^2}{2}\\right)} \, \mathrm{d}t
```
for the number ``z``.
"""
function fresnels(z::Number)
    x = (z * sqrtπ) / 2
    re_x, im_x = reim(x)
    a = (re_x - im_x) + (re_x - im_x) * im
    b = (re_x + im_x) + (im_x - re_x) * im
    re_erf_a, im_erf_a = reim(erf(a))
    re_erf_b, im_erf_b = reim(erf(b))
    re_y = (re_erf_a - im_erf_a + re_erf_b + im_erf_b) / 4
    im_y = (re_erf_a + im_erf_a - re_erf_b + im_erf_b) / 4
    y = re_y + im_y * im
    return y
end
function fresnels(z::Real)
    x = (z * sqrtπ) / 2
    a = x + x * im
    re_erf_a, im_erf_a = reim(erf(a))
    y = (re_erf_a - im_erf_a) / 2
    return y
end

"""
    fresnel(z::Number)

Calculate the normalized sine and cosine fresnel integrals.

See also [`fresnels`](@ref), [`fresnelc`](@ref).
"""
function fresnel(z::Number)
    x = (z * sqrtπ) / 2
    re_x, im_x = reim(x)
    a = (re_x - im_x) + (re_x - im_x) * im
    b = (re_x + im_x) + (im_x - re_x) * im
    re_erf_a, im_erf_a = reim(erf(a))
    re_erf_b, im_erf_b = reim(erf(b))
    re_y_sin = (re_erf_a - im_erf_a + re_erf_b + im_erf_b) / 4
    im_y_sin = (re_erf_a + im_erf_a - re_erf_b + im_erf_b) / 4
    re_y_cos = (re_erf_a + im_erf_a + re_erf_b - im_erf_b) / 4
    im_y_cos = (im_erf_a - re_erf_a + re_erf_b + im_erf_b) / 4
    y_sin = re_y_sin + im_y_sin * im
    y_cos = re_y_cos + im_y_cos * im
    return (y_sin, y_cos)
end
function fresnel(z::Real)
    x = (z * sqrtπ) / 2
    a = x + x * im
    re_erf_a, im_erf_a = reim(erf(a))
    y_sin = (re_erf_a - im_erf_a) / 2
    y_cos = (re_erf_a + im_erf_a) / 2
    return (y_sin, y_cos)
end

