"""
    ErrFunApprox(x,iserfc,result)

Compute the erf and erfc for x::Float64, translated from apple's libm.
"""
function ErrFunApprox(x::Float64, iserfc::Bool, result::Float64)
    if isinf(x)
        if x > 0
            return iserfc ? 0.0 : 1.0
        else
            return iserfc ? 2.0 : -1.0
        end
    end
    if isnan(x)
        return NaN
    end

    # the largest representable double value
    _HUGE = 6.71e+7
    
    InvSqrtPI = 5.6418958354775628695e-1

    a = [3.16112374387056560e+0,
    1.13864154151050156e+2,
    3.77485237685302021e+2,
    3.20937758913846947e+3,
    1.85777706184603153e-1]

    b = [2.36012909523441209e+1,
    2.44024637934444173e+2,
    1.28261652607737228e+3,
    2.84423683343917062e+3]

    ccc = [5.64188496988670089e-1,
    8.88314979438837594e+0,
    6.61191906371416295e+1,
    2.98635138197400131e+2,
    8.81952221241769090e+2,
    1.71204761263407058e+3,
    2.05107837782607147e+3,
    1.23033935479799725e+3,
    2.15311535474403846e-8]

    d = [1.57449261107098347e+1,
    1.17693950891312499e+2,
    5.37181101862009858e+2,
    1.62138957456669019e+3,
    3.29079923573345963e+3,
    4.36261909014324716e+3,
    3.43936767414372164e+3,
    1.23033935480374942e+3]

    pp = [3.05326634961232344e-1,
    3.60344899949804439e-1,
    1.25781726111229246e-1,
    1.60837851487422766e-2,
    6.58749161529837803e-4,
    1.63153871373020978e-2]

    qq = [2.56852019228982242e+0,
    1.87295284992346047e+0,
    5.27905102951428412e-1,
    6.05183413124413191e-2,
    2.33520497626869185e-3]

    y = abs(x)

    if y <= 0.46875e+0
        if y > 1.11e-16
            ysquared = y^2
            numerator=@horner(ysquared,0.,a[3],a[2],a[1],a[5])
            denominator=@horner(ysquared,0.,b[3],b[2],b[1],1.)

            result = y * (numerator + a[4]) / (denominator + b[4]);
        else
            result = y * a[4] / b[4]
        end
        if iserfc
            result = 1.0 - result
        end
        return result
    elseif y <= 4.0
        numerator = ccc[9] * y
        denominator = y
        for i in 1:7
            numerator = (numerator + ccc[i]) * y
            denominator = (denominator + d[i]) * y
        end
        # refer to apple's libm, I don't exactly know the algorithm
        result = (numerator + ccc[8]) / (denominator + d[8])
        ysquared = trunc(y * 16.0) / 16.0
        del = (y - ysquared) * (y + ysquared)
        result = exp(-ysquared^2) * exp(-del) * result
    else
        if y >= _HUGE
            result = InvSqrtPI / y 
            return result
        end
        ysquared = 1.0 / (y^2)
        numerator = pp[6] * ysquared
        denominator = ysquared
        for i in 1:4
            numerator = (numerator + pp[i]) * ysquared
            denominator = (denominator + qq[i]) * ysquared
        end
        result = ysquared * (numerator + pp[5]) / (denominator + qq[5])
        result = (InvSqrtPI - result) / y
        ysquared = trunc(y * 16.0) / 16.0
        del = (y - ysquared) * (y + ysquared)
        result = exp(-ysquared^2) * exp(-del) * result
    end
    return iserfc ? result : (0.5 - result) + 0.5
end

erf(x::Float64) = x >= 0 ? ErrFunApprox(x, false, 1.0) : -ErrFunApprox(x, false, 1.0)
erfc(x::Float64) = x >= 0 ? ErrFunApprox(x, true, 0.0) : 2.0 - ErrFunApprox(x, true, 0.0)
