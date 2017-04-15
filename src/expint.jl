########################################################################
# Inlined, optimized code for the exponential integral E₁ in double precison.
# For more explanations, see course notes by Steven G. Johnson at
#     https://github.com/stevengj/18S096-iap17/blob/master/pset3/pset3-solutions.ipynb

# n coefficients of the Taylor series of E₁(z) + log(z), in type T:
function E₁_taylor_coefficients{T<:Number}(::Type{T}, n::Integer)
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
    taylor = Expr(:macrocall, Symbol("@evalpoly"), :t, c...)
    quote
        let t = $(esc(z))
            $taylor - log(t)
        end
    end
end

# for numeric-literal coefficients: simplify to a ratio of two polynomials:
import Polynomials
# return (p,q): the polynomials p(x) / q(x) corresponding to E₁_cf(x, a...),
# but without the exp(-x) term
function E₁_cfpoly{T<:Real}(n::Integer, ::Type{T}=BigInt)
    q = Polynomials.Poly(T[1])
    p = x = Polynomials.Poly(T[0,1])
    for i = n:-1:1
        p, q = x*p+(1+i)*q, p # from cf = x + (1+i)/cf = x + (1+i)*q/p
        p, q = p + i*q, p     # from cf = 1 + i/cf = 1 + i*q/p
    end
    # do final 1/(x + inv(cf)) = 1/(x + q/p) = p/(x*p + q)
    return p, x*p + q
end
macro E₁_cf64(x, n::Integer)
    p,q = E₁_cfpoly(n, BigInt)
    evalpoly = Symbol("@evalpoly")
    num_expr = Expr(:macrocall, evalpoly, :t, Float64.(Polynomials.coeffs(p))...)
    den_expr = Expr(:macrocall, evalpoly, :t, Float64.(Polynomials.coeffs(q))...)
    quote
        let t = $(esc(x))
            exp(-t) * $num_expr / $den_expr
        end
    end
end

# exponential integral function E₁(z)
function expint(z::Union{Float64,Complex{Float64}})
    x² = real(z)^2
    y² = imag(z)^2
    if x² + 0.233*y² ≥ 7.84 # use cf expansion, ≤ 30 terms
        if (x² ≥ 546121) & (real(z) > 0) # underflow
            return zero(z)
        elseif x² + 0.401*y² ≥ 58.0 # ≤ 15 terms
            if x² + 0.649*y² ≥ 540.0 # ≤ 8 terms
                x² + y² ≥ 4e4 && return @E₁_cf64 z 4
                return @E₁_cf64 z 8
            end
            return @E₁_cf64 z 15
        end
        return @E₁_cf64 z 30
    else # use Taylor expansion, ≤ 37 terms
        r² = x² + y²
        return r² ≤ 0.36 ? (r² ≤ 2.8e-3 ? (r² ≤ 2e-7 ? @E₁_taylor64(z,4) :
                                                       @E₁_taylor64(z,8)) :
                                         @E₁_taylor64(z,15)) :
                          @E₁_taylor64(z,37)
    end
end
expint{T<:Integer}(z::Union{T,Complex{T},Rational{T},Complex{Rational{T}}}) = expint(float(z))

# TODO: in principle, we could implement more-efficient single- and/or half-precision expint
expint(z::Union{Float16,Float32}) = oftype(x, expint(Float64(z)))
expint(z::Union{Complex{Float16},Complex{Float32}}) = oftype(x, expint(Complex128(z)))


######################################################################
######################################################################
# exponential integral Eₙ(z)

# evaluate continued-fraction expansion of E_ν(z) to nterms > 1,
# not including the e⁻ᶻ factor.
function expintcf(ν::Number,z::Number,nterms)
    # (truncation the cf expansion with z seems to minimize the
    # number of terms required for intermediate |z|~3.)
    cf = 1 + nterms*z/(z*z + (ν+nterms))
    for i = nterms-1:-1:1
        cf = 1 + i*cf/(z*cf + (ν+i))
    end
    return cf / (z*cf + ν)
end

function expint(n::Integer, z::Union{Float64,Complex{Float64}})
    if n == 1
        return expint(z)
    elseif n < 1
        # backwards recurrence from E₀ = e⁻ᶻ/z
        zinv = inv(z)
        e⁻ᶻ = exp(-z)
        Eᵢ = zinv * e⁻ᶻ
        for i = 1:-n
            Eᵢ = zinv * (e⁻ᶻ + i * Eᵢ)
        end
        return Eᵢ
    elseif n > 1
        # the recurrence Eᵢ = (e⁻ᶻ - z*Eᵢ) / (i - 1) is inaccurate
        # for large |z| because Eᵢ ≈ e⁻ᶻ/z for large |z|
        x² = real(z)^2
        y² = imag(z)^2
        if x² + 0.233*y² ≥ 7.84 # use cf expansion, ≤ 30 terms
            # TODO: estimate n-dependence of number of terms
            if (x² ≥ 546121) & (real(z) > 0) # underflow
                return zero(z)
            elseif x² + 0.401*y² ≥ 58.0 # ≤ 15 terms
                if x² + 0.649*y² ≥ 540.0 # ≤ 8 terms
                    x² + y² ≥ 4e4 && return expintcf(n,z,4)*exp(-z)
                    return expintcf(n,z,8)*exp(-z)
                end
                return expintcf(n,z,15)*exp(-z)
            end
            return expintcf(n,z,30)*exp(-z)
        else
            # forwards recurrence from E₁
            e⁻ᶻ = exp(-z)
            Eᵢ = expint(z)
            for i = 2:n
                Eᵢ = (e⁻ᶻ - z*Eᵢ) / (i - 1)
            end
            return Eᵢ
        end
    end
end

expint{T<:Integer}(n::Integer, z::Union{T,Complex{T},Rational{T},Complex{Rational{T}}}) = expint(n, float(z))

# TODO: in principle, we could implement more-efficient single- and/or half-precision expint
expint(n::Integer, z::Union{Float16,Float32}) = oftype(x, expint(n, Float64(z)))
expint(n::Integer, z::Union{Complex{Float16},Complex{Float32}}) = oftype(x, expint(n, Complex128(z)))
