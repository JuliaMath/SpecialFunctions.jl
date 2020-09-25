"""
    owens_t(x, a)

Owen's T function

Reference:
    MA Porter, DJ Winstanley, Remark AS R30: A Remark on Algorithm AS76: An Integral Useful in 
    Calculating Noncentral T and Bivariate Normal Probabilities, Applied Statistics, Volume 28, 
    Number 1, 1979, page 113.

    JC Young, Christoph Minder, Algorithm AS 76: An Algorithm Useful in Calculating Non-Central T and
    Bivariate Normal Distributions, Applied Statistics, Volume 23, Number 3, 1974, pages 455-457.
"""
function owens_t(x::Real, a::Real)
    ng = 5
    r = [0.1477621, 0.1346334, 0.1095432, 0.0747257, 0.0333357]
    tp = 0.159155
    tv1 = 1.0E-35
    tv2 = 15.0
    tv3 = 15.0
    tv4 = 1.0E-05
    u = [0.0744372, 0.2166977, 0.3397048, 0.4325317, 0.4869533]
    # Test for X near zero.
    if abs(x) < tv1
        return tp * atan(a)
    end
    # Test for large values of abs(X).
    if tv2 < abs(x)
        return 0.
    end
    # Test for `a` near zero.
    if abs(a) < tv1
        return 0.
    end
    # Test whether abs(`a`) is so large that it must be truncated.
    xs = -0.5 * x * x
    x2 = a
    fxs = a * a
    #  Computation of truncation point by Newton iteration.  
    if tv3 <= log(1.0 + fxs) - xs * fxs
        x1 = 0.5 * a
        fxs = 0.25 * fxs
        while true
            rt = fxs + 1.
            x2 = x1 + (xs * fxs + tv3 - log(rt)) / (2. * x1 * ( 1. / rt - xs))
            fxs = x2 * x2
            if abs(x2 - x1) < tv4
                break
            end
            x1 = x2
        end
    end
    # Gaussian quadrature.
    rt = 0.0
    for i in 1:ng
        r1 = 1.0 + fxs * (0.5 + u[i])^2
        r2 = 1.0 + fxs * (0.5 - u[i])^2
        rt = rt + r[i] * (exp(xs * r1) / r1 + exp(xs * r2) / r2)
    end
    return rt * x2 * tp
end
