import Polynomials
using Base.MathConstants

function E₁_cfpoly_approx(n::Integer, pstart::Polynomials.Polynomial{T}, ::Type{T}=BigInt) where {T<:Real}
    q = Polynomials.Polynomial(T[1])
    p = pstart
    x = Polynomials.Polynomial(T[0,1])
    for i = n:-1:1
        p, q = x*p+(1+i)*q, p # from cf = x + (1+i)/cf = x + (1+i)*q/p
        p, q = p + i*q, p     # from cf = 1 + i/cf = 1 + i*q/p
    end
    # do final 1/(x + inv(cf)) = 1/(x + q/p) = p/(x*p + q)
    return p, x*p + q
end

macro E₁_cf64(x, n::Integer, start)
    pstart = Polynomials.Polynomial(eval(start))
    # consider using BigFloat?
    p, q = E₁_cfpoly_approx(n, pstart, Float64)
    xesc = esc(x)
    
    num_expr =  :(@evalpoly $xesc)
    append!(num_expr.args, Float64.(Polynomials.coeffs(p)))
    den_expr = :(@evalpoly $xesc)
    append!(den_expr.args, Float64.(Polynomials.coeffs(q)))
    :( exp(-$xesc) * $num_expr / $den_expr )
end


function E₁_taylor_coefficients(::Type{T}, n::Integer) where {T<:Number}
    n < 0 && throw(ArgumentError("$n ≥ 0 is required"))
    n == 0 && return T[]
    n == 1 && return T[-eulergamma]
    # iteratively compute the terms in the series, starting with k=1
    term::T = 1
    terms = T[-eulergamma, term]
    for k=2:n
        term = -term * (k-1) / (k * k)
        push!(terms, term)
    end
    return terms
end

# inline the Taylor expansion for a given order n, in double precision
macro E₁_taylor64(z, n::Integer)
    c = E₁_taylor_coefficients(Float64, n)
    zesc = esc(z)
    taylor = :(@evalpoly $zesc)
    append!(taylor.args, c)
    :( $taylor - log($zesc) )
end

E₁(x::Real) = E₁(float(x))

function E₁(x::Float64)
    x < 0 && throw(DomainError(x, "negative argument"))
    x == 0 && return Inf
    if x > 2.15
        # specially chosen approximants for faster convergence
        x < 3.0   && return @E₁_cf64(x, 18, [6.267445506556548, -2.278962735947262, 0.5577231261815463, -0.05413049191473329])
        x < 4.0   && return @E₁_cf64(x, 16, [5.114292670961982, -1.2789140459431323, 0.22066200334871455, -0.015067049382830766])
        x < 6.1   && return @E₁_cf64(x, 14, [4.194988480897909, -0.7263593325667503, 0.08956574399359891, -0.00434973529065973])
        x < 8.15  && return @E₁_cf64(x, 9,  [3.0362016309948228, -0.33793806630590445, 0.029410409377178114, -0.0010060498260648586])
        x < 25.0  && return @E₁_cf64(x, 8,  [2.5382065303376895, -0.18352177433259526, 0.011141562002742184, -0.0002634921890930066])
        x < 200.0 && return @E₁_cf64(x, 8,  [0.0, 1.0])
        return x < 740.0 ? @E₁_cf64(x, 4, [0.0, 1.0]) : 0.0 # underflow
    else
        # crossover point to taylor should be tuned more
        return x ≤ 0.6 ? (x ≤ 0.053 ? (x ≤ 4.4e-3 ? @E₁_taylor64(x,4) :
                                                       @E₁_taylor64(x,8)) :
                                       @E₁_taylor64(x,15)) :
                          @E₁_taylor64(x,37)
    end
end


function En_cf(x::Number, ν::Number, n::Int=1000)
    B = float(x + ν)
    Bprev::typeof(B) = x
    A::typeof(B) = 1
    Aprev::typeof(B) = 1
    
    # two recurrence steps / loop
    iters = 0
    for i = 2:n
        iters += 1

        A′ = A
        A = x*A + (i-1) * Aprev
        Aprev = A′
        B′ = B
        B = x*B + (i-1) * Bprev
        Bprev = B′
        
        A′ = A
        A = A + (ν+i-1) * Aprev
        Aprev = A′
        B′ = B
        B = B + (ν+i-1) * Bprev
        Bprev = B′
        
        conv = abs(Aprev*B - A*Bprev) < 1e-15 * abs(B*Bprev)
        conv && break
        
        # rescale 
        if max(abs(real(A)), abs(imag(A))) > 1e50
            A /= 1e50
            Aprev /= 1e50
            B /= 1e50
            Bprev /= 1e50
        end
    end
    
    return A/B * exp(-x), iters
end


# Compute expint(ν, z₀+Δ) given start = expint(ν, z₀)
function En_taylor(ν::Number, start::Number, z₀::Number, Δ::Number)
    a = exp(z₀) * start
    k, iters = 0, 0
    asum = a
    Δ_prod_fact = -Δ
    
    while iters < 100 # perhaps
        a_pre = Δ_prod_fact + a*Δ*(ν - k - 1)/(k+1)
        a = a_pre / z₀
        asum_prev = asum
        asum += a
        
        if abs(asum_prev - asum) < 1e-15
            break
        end
        
        Δ_prod_fact *= -Δ / (k+2)
        
        iters += 1
        k += 1
    end

    res = exp(-z₀) * asum
    return res
end

# series about origin
# https://functions.wolfram.com/GammaBetaErf/ExpIntegralE/06/01/04/01/01/0003/
function En_expand_origin(x::Number, ν::Number)
    if isreal(ν) && floor(ν) == ν
        # go to special case for integer ν
        return En_expand_origin(x, Int(ν))
    end
    gammaterm = gamma(1 - ν) * x^(ν-1)
    frac = 1
    sumterm = frac / (1 - ν)
    k, maxiter = 1, 100
    while k < maxiter
        frac *= -x / k
        prev = sumterm
        sumterm += frac / (k + 1 - ν)
        if abs(sumterm - prev) < 1e-15 * abs(prev)
            break
        end
        k += 1
    end
    
    return gammaterm - sumterm
end

# series about the origin, special case for integer n
function En_expand_origin(x::Number, n::Integer)
    gammaterm = 1
    for i = 1:n-1
        gammaterm *= -x / i
    end

    gammaterm *= digamma(n) - log(x)
    sumterm = n == 1 ? 0 : 1 / (1 - n)
    frac = 1
    k, maxiter = 1, 100
    while k < maxiter
        frac *= -x / k
        # skip term with zero denominator
        if k != n-1
            prev = sumterm
            sumterm += frac / (k + 1 - n)
            if abs(sumterm - prev) < 1e-15 * abs(prev)
                break
            end
        end
        k += 1
    end
    return gammaterm - sumterm
end

# can find imaginary part of E_ν(x) for x on negative real axis analytically
function imagbranchcut(x::Number, ν::Number)
    a = real(x)
    impart = π * im * exp(-π*im*ν) * a^(ν-1) / gamma(ν)
    # exp(n*log(z) - loggamma(n))
    return imag(impart) * im # get rid of any real error
end

const ORIGIN_EXPAND_THRESH = 3
function En(ν::Number, x::Number, niter::Int=1000, debug=false)
    if x == 0.0
        if real(ν) > 0
            return 1.0 / (ν - 1)
        else
            return Inf
        end
    end
    if ν == 0
        return exp(-x) / x
    end
    
    if abs(x) < ORIGIN_EXPAND_THRESH
        return En_expand_origin(x, ν)
    end
    if real(x) > 0
        res, i = En_cf(x, ν, niter)
        return res
    elseif real(x) < 0
        doconj = imag(x) < 0
        rex, imx = real(x), abs(imag(x))
        x = doconj ? conj(x) : x
        ν = doconj ? conj(ν) : ν
        
        # empirical boundary for 500 iterations
        boundary = min(1 + 0.5*abs(ν), 50)
        if imx > boundary
            res, i = En_cf(x, ν, niter)
            return doconj ? conj(res) : res
        else
            # iterate with taylor
            # first find starting point
            # TODO: switch back to CF for large ν
            imstart = boundary
            z₀ = rex + imstart*im
            E_start, i = En_cf(z₀, ν, niter)
            while i >= niter
                imstart *= 2
                z₀ = rex + imstart*im
                E_start, i = En_cf(z₀, ν, niter)
            end
            
            # nsteps chosen so |Δ| ≤ 0.5
            nsteps = ceil(2 * (imstart - imx))
            Δ = (imx - imstart)*im / nsteps

            for j = 1:nsteps
                E_start = En_taylor(ν, E_start, z₀, Δ)
                z₀ += Δ
            end
            
            # more exact imaginary part available for non-integer ν
            if imx == 0
                E_start = real(E_start) + imagbranchcut(x, ν)
            end
            
            return doconj ? conj(E_start) : E_start
        end
    end
end
