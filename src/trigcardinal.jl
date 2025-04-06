# u corresponds to unnormalized

# sinc/sincu is zero when the real part is Inf and imag is finite
isinf_real(x::Real) = isinf(x)
isinf_real(x::Number) = isinf(real(x)) && isfinite(imag(x))

# sinhc/sinhcu is zero when the imag part is Inf and real is finite
isinf_imag(x::Real) = false
isinf_imag(x::Number) = isfinite(real(x)) && isinf(imag(x))

# sincu copied exactly from Boost library
# https://www.boost.org/doc/libs/1_87_1/boost/math/special_functions/sinc.hpp
sincu(x) = _sinc(float(x))
function _sinc(x::Union{T,Complex{T}}) where {T}
    if isinf_real(x)
        return zero(x)
    end
    
    nrm = fastabs(x)
    if nrm >= 3.3*sqrt(sqrt(eps(T)))
        return sin(x)/x
    else
        # |x| < (eps*120)^(1/4)
        return 1 - x*x/6
    end
end

# sinhcu copied exactly from Boost library 
# https://www.boost.org/doc/libs/1_87_1/boost/math/special_functions/sinhc.hpp
sinhcu(x) = _sinhcu(float(x))
function _sinhcu(x::Union{T,Complex{T}}) where {T}
    taylor_0_bound = eps(T)
    taylor_2_bound = sqrt(taylor_0_bound)
    taylor_n_bound = sqrt(taylor_2_bound)

    if isinf_imag(x) 
        return zero(x)
    end
    
    nrm = fastabs(x)

    if nrm >= taylor_n_bound || isnan(nrm)
        return sinh(x)/x
    else
        # approximation by taylor series in x at 0 up to order 0
        res = one(x)
        if nrm >= taylor_0_bound
            x2 = x*x
            # approximation by taylor series in x at 0 up to order 2
            res += x2/6
            if nrm >= taylor_2_bound
                # approximation by taylor series in x at 0 up to order 4
                res += (x2*x2)/120
            end
        end
        return res
    end
end
