

# Code ported from https://github.com/kiranshila/FresnelIntegrals.jl


"""
        fresnelc(z)
Calculates the Fresnel cosine integral for the number z for
    ``C(z) = \\int_{0}^{z} \\cos{\\left(\\frac{\\pi t^2}{2}\\right)}dt``
"""
fresnelc(z::Number) = 0.25*(1-1im)*(1im*erf(0.5*(1-1im)*z*√(π)) + erf(0.5*(1+1im)*z*√(π)))

"""
        fresnels(z)
Calculates the Fresnel sine integral for the number z for
    ``S(z) = \\int_{0}^{z} \\sin{\\left(\\frac{\\pi t^2}{2}\\right)}dt``
"""
fresnels(z::Number) = 0.25*(1+1im)*(-1im*erf(0.5*(1-1im)*z*√(π)) + erf(0.5*(1+1im)*z*√(π)))

"""
        fresnel(z)
Calculates the cosine and sine fresnel integrals
"""
fresnel(z::Number) = (fresnelc(z),fresnels(z))

