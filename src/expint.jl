using Base.MathConstants

# barebones polynomial type supporting +, -, *
struct SimplePoly{T}
    coeffs::Vector{T}
end

Base.:*(p::SimplePoly, c) = SimplePoly(p.coeffs * c)
Base.:*(c, p::SimplePoly) = p * c

function Base.:+(p::SimplePoly{S}, q::SimplePoly{T}) where {S, T}
    n, m = length(p.coeffs), length(q.coeffs)
    ext1, ext2 = max(0, m - n), max(0, n - m)
    c1 = [p.coeffs ; zeros(S, ext1)]
    c2 = [q.coeffs ; zeros(T, ext2)]
    return SimplePoly(c1 + c2)
end

Base.:-(p::SimplePoly) = -1 * p
Base.:-(p::SimplePoly, q::SimplePoly) = p + (-q)

function Base.:*(p::SimplePoly{S}, q::SimplePoly{T}) where {S, T}
    n, m = length(p.coeffs) - 1, length(q.coeffs) - 1
    c = zeros(promote_type(S, T), m + n + 1)
    for i in 0:n, j in 0:m
        c[i + j + 1] += p.coeffs[i + 1] * q.coeffs[j + 1]
    end
    return SimplePoly(c)
end

function E₁_cfpoly_approx(n::Integer, ::Type{T}=BigInt) where {T<:Real}
    q = SimplePoly(T[1])
    p = x = SimplePoly(T[0,1])
    for i = n:-1:1
        p, q = x*p+(1+i)*q, p # from cf = x + (1+i)/cf = x + (1+i)*q/p
        p, q = p + i*q, p     # from cf = 1 + i/cf = 1 + i*q/p
    end
    # do final 1/(x + inv(cf)) = 1/(x + q/p) = p/(x*p + q)
    return p, x*p + q
end

macro E₁_cf64(x, n::Integer)
    # consider using BigFloat?
    p, q = E₁_cfpoly_approx(n, Float64)
    xesc = esc(x)

    num_expr =  :(@evalpoly $xesc)
    append!(num_expr.args, Float64.(p.coeffs))
    den_expr = :(@evalpoly $xesc)
    append!(den_expr.args, Float64.(q.coeffs))
    :( $num_expr / $den_expr )
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

# adapted from Steven G Johnson's initial implementation: issue #19
function expint_opt(x::Float64, ::Val{expscaled}=Val{false}()) where {expscaled}
    x < 0 && throw(DomainError(x, "negative argument, convert to complex first"))
    x == 0 && return Inf
    # E_1(inf) * exp(inf) is indeterminate
    x == Inf && return expscaled ? NaN : zero(x)
    if x > 2.15
        mult = expscaled ? 1 : exp(-x)
        # minimax rational approximations to E_1(x)/exp(-x):
        x < 4.0     && return mult *
            @evalpoly(x, 3.600530862438501481559423277418128014798, 28.73031134165011013783185685393062481126, 46.04314409968065653003548224846863877481, 21.47189493062368074985000918414086604187, 2.719957622921613321844755385973197500235, 1.508750885580864752293599048121678438568e-6) /
            @evalpoly(x, 1.0, 18.06743589038646654075831055159865459831, 61.19456872238615922515377354442679999566, 64.81772518730116701207299231777089576859, 24.19034591054828198408354214931112970741, 2.720026796991567940556850921390829046015)
        x < 10.0    && return mult *
            @evalpoly(x, 3.149019890512432117647119992448352099575, 14.77395058090815888166486507990452627136, 14.78214309058953358717796744960600201013, 4.559401130686434886620102186841739864936, 0.4027858394909585103775445204576054721422, 2.302781920509468929446800686773538387432e-9) /
            @evalpoly(x, 1.0, 11.65960373479520446458792926669115987821, 26.20023773503894444535165299118172674849, 18.93899893550582921168134979000319186841, 4.962178168140565906794561070524079194193, 0.4027860481050182948737116109665534358805)
        x < 20.0    && return mult *
            @evalpoly(x, 2.564801308922428705318296668924369454617, 5.482252510134574167659359513298970499768, 2.379528224853089764405551768869103662657, 0.2523431276282591480166881146593592650031, 1.444719769329975045925053905197199934930e-9, -8.977332061106791470409502623202677045468e-12) /
            @evalpoly(x, 1.0, 6.421214405318108272004472721910023284626, 7.609584052290707052877352911548283916247, 2.631866613851763364839413823973711355890, 0.2523432345539146248733539983749722854603)
        x < 200.0 && return mult * @E₁_cf64(x, 8) # fallback continued fraction
        return mult * @E₁_cf64(x, 4)
    else
        # crossover point to taylor should be tuned more
        mult = expscaled ? exp(x) : 1
        return mult * (x ≤ 0.6 ? (x ≤ 0.053 ? (x ≤ 4.4e-3 ? @E₁_taylor64(x,4) :
                                                       @E₁_taylor64(x,8)) :
                                       @E₁_taylor64(x,15)) :
                          @E₁_taylor64(x,37))
    end
end

function expint_opt(z::Complex{Float64}, ::Val{expscaled}=Val{false}()) where {expscaled}
    x² = real(z)^2
    y² = imag(z)^2
    if x² + 0.233*y² ≥ 7.84 # use cf expansion, ≤ 30 terms
        mult = expscaled ? 1 : exp(-z)
        #if (x² ≥ 546121) & (real(z) > 0) # underflow
        #    return zero(z)
        if x² + 0.401*y² ≥ 58.0 # ≤ 15 terms
            if x² + 0.649*y² ≥ 540.0 # ≤ 8 terms
                x² + y² ≥ 4e4 && return mult * @E₁_cf64(z, 4)
                return mult * @E₁_cf64(z, 8)
            end
            return mult * @E₁_cf64(z, 15)
        end
        return mult * @E₁_cf64(z, 30)
    else # use Taylor expansion, ≤ 37 terms
        r² = x² + y²
        mult = expscaled ? exp(z) : 1
        return mult * (r² ≤ 0.36 ? (r² ≤ 2.8e-3 ? (r² ≤ 2e-7 ? @E₁_taylor64(z,4) :
                                                       @E₁_taylor64(z,8)) :
                                         @E₁_taylor64(z,15)) :
                          @E₁_taylor64(z,37))
    end
end

function _expint(z::Complex{Float64}, ::Val{expscaled}=Val{false}()) where {expscaled}
    if real(z) < 0
        return _expint(1, z, 1000, Val{expscaled}())
    else
        return expint_opt(z, Val{expscaled}())
    end
end
expint(z::Complex{Float64}) = _expint(z)
expintx(z::Complex{Float64}) = _expint(z, Val{true}())
expint(x::Float64) = expint_opt(x)
expintx(x::Float64) = expint_opt(x, Val{true}())

function expint(x::BigFloat)
    iszero(x) && return Inf
    signbit(x) && throw(DomainError(x, "negative argument, convert to complex first"))
    return -expinti(-x)
end

expint(x::Number) = expint(1, x)
expint(z::Union{T,Complex{T},Rational{T},Complex{Rational{T}}}) where {T<:Integer} = expint(float(z))
for f in (:expint, :expintx)
    for T in (Float16, Float32)
        @eval $f(x::$T) = $T($f(Float64(x)))
    end
    for CT in (ComplexF16, ComplexF32)
        @eval $f(x::$CT) = $CT($f(ComplexF64(x)))
    end
end

# Continued fraction for En(ν, z) that doesn't use a term with
# the gamma function: https://functions.wolfram.com/GammaBetaErf/ExpIntegralE/10/0001/
function En_cf_nogamma(ν::Number, z::Number, n::Int=1000)
    B = float(z + ν)
    Bprev::typeof(B) = z
    A::typeof(B) = 1
    Aprev::typeof(B) = 1
    ϵ = 10*eps(real(B))
    scale = sqrt(floatmax(real(A)))

    # two recurrence steps / loop
    iters = 1
    for i = 2:n
        iters += 1

        A, Aprev = z*A + (i-1) * Aprev, A
        B, Bprev = z*B + (i-1) * Bprev, B

        (isinf(A) | isinf(B)) && return Aprev/Bprev, iters-1

        A, Aprev = A + (ν+i-1) * Aprev, A
        B, Bprev = B + (ν+i-1) * Bprev, B

        i > 4 && abs(Aprev*B - A*Bprev) < ϵ*abs(B*Bprev) && break

        # rescale
        if fastabs(A) > scale
            A     /= scale
            Aprev /= scale
            B     /= scale
            Bprev /= scale
        end
    end

    return A/B, iters
end

# Calculate Γ(1 - ν) * z^(ν-1) safely
function En_safe_gamma_term(ν::Number, z::Number)
    ν1 = 1 - oftype(z, ν)
    lgamma, lgammasign = ν1 isa Real ? logabsgamma(ν1) : (loggamma(ν1), 1)
    return lgammasign * exp((ν - 1)*log(z) + lgamma)
end
En_safe_gamma_term(ν::Integer, z::Real) = (z ≥ 0 || isodd(ν) ? 1 : -1) * exp((ν - 1)*log(abs(z)) + loggamma(1 - oftype(z, ν)))

# continued fraction for En(ν, z) that uses the gamma function:
# https://functions.wolfram.com/GammaBetaErf/ExpIntegralE/10/0005/
# returns the two terms from the above equation separately
function En_cf_gamma(ν::Number, z::Number, n::Int=1000)
    A = float(1 - ν)
    B::typeof(A) = 1
    Bprev::typeof(A) = 0
    Aprev::typeof(A) = 1
    ϵ = 10*eps(real(B))
    scale = sqrt(floatmax(real(A)))

    iters = 0
    for i = 1:n
        iters += 1

        a = iseven(i) ? (i÷2)*z : -((i + 1)÷2 - ν)*z
        b = (i + 1 - ν)

        A, Aprev = b*A + a*Aprev, A
        B, Bprev = b*B + a*Bprev, B

        if abs(Aprev*B - A*Bprev) < ϵ*abs(Aprev*B)
            break
        end

        if max(abs(real(A)), abs(imag(A))) > scale
            A     /= scale
            Aprev /= scale
            B     /= scale
            Bprev /= scale
        end
    end

    gammapart = En_safe_gamma_term(ν, z)
    return gammapart, -B/A, iters
end

# picks between continued fraction representations in
# En_cf_nogamma and En_cf_gamma
# returns (evaluated result, # iterations used, whether En_cf_gamma was chosen)
function En_cf(ν::Number, z::Number, niter::Int=1000)
    if real(1-ν) > 0
        gammapart, cfpart, iters = En_cf_gamma(ν, z, niter)
        gammaabs, cfabs = abs(gammapart), abs(cfpart)
        if gammaabs != Inf && gammaabs > 1.0 && gammaabs > cfabs
            # significant gamma part, use this
            return gammapart, cfpart, iters, true
        end
    end
    return zero(z), En_cf_nogamma(ν, z, niter)..., false
end

# Compute expint(ν, z₀+Δ) given start = expint(ν, z₀), as described by [Amos 1980].
# This is used to incrementally approach the negative real axis.
function En_taylor(ν::Number, start::Number, z₀::Number, Δ::Number)
    a = exp(z₀) * start
    k, iters = 0, 0
    asum = a
    Δ_prod_fact = -Δ
    ϵ = 10*eps(real(asum))

    for k = 0:100
        a_pre = Δ_prod_fact + a*Δ*(ν - k - 1)/(k + 1)
        a = a_pre / z₀
        asum_prev = asum
        asum += a

        if abs(asum_prev - asum) < ϵ
            break
        end

        Δ_prod_fact *= -Δ / (k + 2)
    end

    res = exp(-z₀) * asum
    return res
end

# series about origin, general ν
# https://functions.wolfram.com/GammaBetaErf/ExpIntegralE/06/01/04/01/01/0003/
function En_expand_origin_general(ν::Number, z::Number, niter::Integer)
    # gammaterm = En_safe_gamma_term(ν, z)
    gammaterm = gamma(1-ν)*z^(ν-1)
    frac = one(z)
    blowup  = abs(1 - ν) < 0.5 ? frac / (1 - ν) : zero(z)
    sumterm = abs(1 - ν) < 0.5 ? zero(z) : frac / (1 - ν)
    k = 1
    ε = 10*eps(typeof(abs(frac)))
    while k < niter
        frac *= -z / k
        prev = sumterm
        if abs(k + 1 - ν) < 0.5
            blowup += frac / (k + 1 - ν)
        else
            sumterm += frac / (k + 1 - ν)
            if abs(sumterm - prev) < ε*abs(prev)
                break
            end
        end
        k += 1
    end

    if real(ν+z) isa Union{Float64, Float32} && abs(gammaterm - blowup) < 1e-3 * abs(blowup)
        δ = round(ν) - ν
        n = real(round(ν)) - 1

        # (1 - z^δ)/δ series
        logz = log(z)
        series1 = -logz - logz^2*δ/2 - logz^3*δ^2/6 - logz^4*δ^3/24 - logz^5*δ^4/120

        # due to https://functions.wolfram.com/GammaBetaErf/Gamma/06/01/05/01/0004/
        # expressions for higher order terms found using:
        # https://gist.github.com/augustt198/348e8f9ba33c0248f1548309c47c6d0e
        ψ₀, ψ₁, ψ₂, ψ₃, ψ₄ = polygamma.((0,1,2,3,4), n+1)
        series2 = ψ₀ + (3*ψ₀^2 + π^2 - 3*ψ₁)*δ/6 + (ψ₀^3 + (π^2 - 3ψ₁)*ψ₀ + ψ₂)δ^2/6
        series2 += (7π^4 + 15*(ψ₀^4 + 2ψ₀^2 * (π^2 - 3ψ₁) + ψ₁*(-2π^2 + 3ψ₁) + 4ψ₀*ψ₂) - 15ψ₃)*δ^3/360
        series2 += (3ψ₀^5 + ψ₀^3*(10π^2 - 30ψ₁) + 30ψ₀^2*ψ₂ + ψ₀*(45ψ₁^2 - 30π^2*ψ₁ - 15ψ₃ + 7π^4) - 30ψ₁*ψ₂ + 10π^2*ψ₂ + 3ψ₄)*δ^4/360

        return (series1 + series2) * En_safe_expfact(n, z) * z^(ν-n-1) - sumterm
    end
    return gammaterm - (blowup + sumterm)
end

# compute (-z)^n / n!, avoiding overflow if possible, where n is an integer ≥ 0 (but not necessarily an Integer)
function En_safe_expfact(n::Real, z::Number)
    if n < 100
        powerterm = one(z)
        for i = 1:Int(n)
            powerterm *= -z/i
        end
        return powerterm
    else
        if z isa Real
            sgn = z ≤ 0 ? one(n) : (n <= typemax(Int) ? (isodd(Int(n)) ? -one(n) : one(n)) : (-1)^n)
            return sgn * exp(n * log(abs(z)) - loggamma(n+1))
        else
            return exp(n * log(-z) - loggamma(n+1))
        end
    end
end

# series about the origin, special case for integer n > 0
# https://functions.wolfram.com/GammaBetaErf/ExpIntegralE/06/01/04/01/02/0005/
function En_expand_origin_posint(n, z::Number, niter::Integer)
    gammaterm = En_safe_expfact(n-1, z) # (-z)^(n-1) / (n-1)!
    frac = one(real(z))
    gammaterm *= digamma(oftype(frac,n)) - log(z)
    sumterm = n == 1 ? zero(frac) : frac / (1 - n)
    k = 1
    ϵ = 10*eps(real(sumterm))
    while k < niter
        frac *= -z / k
        # skip term with zero denominator
        if k != n-1
            prev = sumterm
            sumterm += frac / (k + 1 - n)
            if abs(sumterm - prev) < ϵ
                break
            end
        end
        k += 1
    end
    return gammaterm - sumterm
end

function En_expand_origin(ν::Number, z::Number, niter::Integer)
    if isinteger(ν) && real(ν) > 0
        return real(ν) < (typemax(Int)>>2) ? En_expand_origin_posint(Int(real(ν)), z, niter) : En_expand_origin_posint(real(ν), z, niter)
    else
        return En_expand_origin_general(ν, z, niter)
    end
end

# can find imaginary part of E_ν(x) for x on negative real axis analytically
# https://functions.wolfram.com/GammaBetaErf/ExpIntegralE/04/05/01/0003/
function En_imagbranchcut(ν::Number, z::Number)
    a = real(z)
    e1 = exp(oftype(a, π) * imag(ν))
    e2 = Complex(cospi(real(ν)), -sinpi(real(ν)))
    lgamma, lgammasign = ν isa Real ? logabsgamma(ν) : (loggamma(ν), 1)
    return -2 * lgammasign * e1 * π * e2 * exp((ν-1)*log(complex(a)) - lgamma) * im
end

function En_safeexpmult(z, a)
    zexp = exp(z)
    if isinf(zexp) || iszero(zexp)
        return a isa Real ? sign(a) * exp(z + log(abs(a))) : exp(z + log(a))
    else
        return zexp*a
    end
end

function _expint(ν::Number, z::Number, niter::Int=1000, ::Val{expscaled}=Val{false}()) where {expscaled}
    if abs(ν) > 50 && !(isreal(ν) && real(ν) > 0)
        throw(ArgumentError("Unsupported order |ν| > 50 off the positive real axis"))
    end

    z, = promote(float(z), ν)
    if isnan(ν) || isnan(z)
        return oftype(z, NaN) * z
    end

    if z isa Real && (isinteger(ν) ? (z < 0 && ν > 0) : z < 0)
        throw(DomainError(z, "En will only return a complex result if called with a complex argument"))
    end

    if z == 0
        typ = typeof(real(z))
        return oftype(z, real(ν) > 0 ? one(typ)/(ν-1) : convert(typ, Inf))
    end

    if ν == 0
        return expscaled ? inv(z) : exp(-z)/z
    elseif ν == 1 && real(z) > 0 && z isa Union{Float64, Complex{Float64}}
        return expint_opt(z, Val{expscaled}())
    end
    # asymptotic test for underflow when Re z → ∞
    # https://functions.wolfram.com/GammaBetaErf/ExpIntegralE/06/02/0003/
    if !expscaled && real(z) > -log(nextfloat(zero(real(z))))+1 # exp(-z) is zero
        return zero(z)
    end

    if abs2(z) < 9
        # use Taylor series about the origin for small z
        mult = expscaled ? exp(z) : 1
        return mult * En_expand_origin(ν, z, niter)
    end

    if z isa Real || real(z) > 0
        gmult = expscaled ? exp(z) : 1
        if real(z) > 0
            g, cf, _ = En_cf(ν, z, niter)
        else
            g, cf, _ = zero(z), En_cf_nogamma(ν, z, niter)...
        end

        g != 0 && (g *= gmult)
        cf = expscaled ? cf : En_safeexpmult(-z, cf)
        return g + cf
    else
        # Procedure for z near the negative real axis based on
        # (without looking at the accompanying source code):
        #   Amos, D. E. (1990). Computation of exponential integrals
        #   of a complex argument. ACM Transactions on Mathematical Software,
        #   16(2), 169–177. https://doi.org/10.1145/78928.78933
        doconj = imag(z) < 0
        rez, imz = real(z), abs(imag(z))
        zorig = z
        z = doconj ? conj(z) : z
        ν = doconj ? conj(ν) : ν

        quick_niter = niter >> 4
        # start with small imaginary part if exactly on negative real axis
        imstart = (imz == 0) ? abs(z)*sqrt(eps(typeof(real(z)))) : imz
        z₀ = rez + imstart*im
        g_start, cf_start, i, _ = En_cf(ν, z₀, quick_niter)
        E_start = g_start + En_safeexpmult(-z₀, cf_start)
        isnan(E_start) && return E_start
        if imz > 0 && i < quick_niter
            # didn't need to take any steps
            if expscaled
                E_start = En_safeexpmult(z₀, g_start) + cf_start
            end
            return doconj ? conj(E_start) : E_start
        end
        while i == quick_niter
            # double imaginary part until in region with fast convergence
            imstart *= 2
            z₀ = rez + imstart*im
            g_start, cf_start, i, _ = En_cf(ν, z₀, quick_niter)
            E_start = g_start + En_safeexpmult(-z₀, cf_start)
        end

        # nsteps chosen so |Δ| ≤ 0.5
        nsteps = ceil(2 * (imstart - imz))
        Δ = (imz - imstart)*im / nsteps

        for j = 1:nsteps
            # take Δ sized steps towards the desired z
            E_start = En_taylor(ν, E_start, z₀, Δ)
            z₀ += Δ
        end

        En = doconj ? conj(E_start) : E_start

        # handle branch cut
        if imz == 0
            bc = En_imagbranchcut(ν, z)
            bit = !signbit(imag(z))
            sign = bit ? 1 : -1
            if isreal(ν)
                # can separate real/im in case of real ν
                En = real(En) - sign * imag(bc)/2 * im
            else
                En = bit ? En : En + bc
            end
        end
        return expscaled ? En_safeexpmult(zorig, En) : En
    end
end

"""
    expint(z)
    expint(ν, z)

Computes the exponential integral ``\\operatorname{E}_\\nu(z) = \\int_1^\\infty \\frac{e^{-zt}}{t^\\nu} dt``.
If ``\\nu`` is not specified, ``\\nu=1`` is used. Arbitrary complex ``\\nu`` and ``z`` are supported.

External links: [DLMF](https://dlmf.nist.gov/8.19), [Wikipedia](https://en.wikipedia.org/wiki/Exponential_integral)
"""
expint(ν::Number, z::Number, niter::Int=1000) = _expint(ν, z, niter, Val{false}())


"""
    expintx(z)
    expintx(ν, z)

Computes the scaled exponential integral ``\\exp(z) \\operatorname{E}_\\nu(z) = e^z \\int_1^\\infty \\frac{e^{-zt}}{t^\\nu} dt``.
If ``\\nu`` is not specified, ``\\nu=1`` is used. Arbitrary complex ``\\nu`` and ``z`` are supported.

See also: [`expint(ν, z)`](@ref SpecialFunctions.expint)
"""
expintx(ν::Number, z::Number, niter::Int=1000) = _expint(ν, z, niter, Val{true}())

##############################################################################
# expinti function Ei

"""
    expinti(x::Real)

Computes the exponential integral function ``\\operatorname{Ei}(x) = \\int_{-\\infty}^x \\frac{e^t}{t} dt``,
which is equivalent to ``-\\Re[\\operatorname{E}_1(-x)]`` where ``\\operatorname{E}_1`` is the `expint` function.
"""
expinti(x::Real) = x ≤ 0 ? -expint(-x) : -real(expint(complex(-x)))

# inline the Taylor expansion of Ei for a given order n, in double precision
macro Ei_taylor64(z, n::Integer)
    c = .- E₁_taylor_coefficients(Float64, n) .* (-1).^(0:n)
    zesc = esc(z)
    taylor = :(@evalpoly $zesc)
    append!(taylor.args, c)
    :( $taylor + log($zesc) )
end

# optimized double-precision version
function expinti(x::Float64)
    x ≤ 0 && return -expint(-x)

    if x > 2.15
        # minimax rational approximants:
        x < 4 && return @evalpoly(x,-2.43791466332154621,3.09402100064798205,9.35202477109609955,0.152659977028953397,0.0157273683896079142,0.0345566671015011426,-0.000421531433157416203) /
                        @evalpoly(x,1.0,4.28055563991564399,0.537599625698465573,-0.511064414527643313,0.0867748262262250088,-0.00623913330836521800,0.000172066498182538260)
        x < 10 && return exp(x) *
            @evalpoly(x,-1.58447018083420958,4.71806833998906997,-0.587691572500210206,0.125012472861504555,-0.00178055441724967428,0.000633648975971195928,0.0000147213934578379204,2.12391754244415544e-6) /
            @evalpoly(x,1.0,1.93297600031287800,0.660790440069106542,0.198322636197663277,0.0272447293513279631,0.00399501571688512611,0.000362510989191199243,0.0000182930089855534336,2.06800780072894204e-6)
        x < 20 && return exp(x) *
            @evalpoly(x,-1.77183291754640123,0.795659966861260409,-0.221223333413388642,0.0328877243243796815,-0.00331846947191676458,0.000180945604349930285,-5.97641401680304362e-6,2.42151808626299747e-11) /
            @evalpoly(x,1.0,-2.10926998628216150,0.933357955421497965,-0.245433884954174394,0.0356954809772243699,-0.00348034743685382360,0.000186615220328647350,-5.97232581353392099e-6)

        x > 710 && return Inf # overflow
        x⁻¹ = inv(x)

        # minimax rational approximant in x⁻¹ for x ∈ [20,200]
        x < 200 && return exp(x) *
            @evalpoly(x⁻¹, -5.29842699621003563e-14, +1.00000000004732488, -60.4361334939888359, +1327.83891720487710, -6810.63668974273961, -177755.383525765400,+3.00773484037048848e6, -1.53642380695372707e7, +2.08174653368702692e7) /
            @evalpoly(x⁻¹, 1.0, -61.4361334756161381, +1387.27504658395142, -8081.03888544858393, -172104.333927401741, +3.18903576285551101e6, -1.81873890267574206e7, +3.37312131843327704e7, -1.22198734384213631e7)

        # asymptotic series eˣ/x ∑k!/xᵏ for x ≥ 200:
        return exp(x)*x⁻¹ * @evalpoly(x⁻¹, 1,1,2,6,24,120,720,5040)
    end

    root = 0.37250741078136663446 # OEIS A091723
    dx = x - root
    if abs(dx) < 0.03
        # taylor series around real root of Ei
        return dx * @evalpoly(dx, 3.896215733907167310, -3.281607866398561671, 6.52237614543892570, -12.96969738353651704, 27.88629796294204998, -62.3788015289154187, 143.5349488096750988, -337.155827178746892, 804.531839982138251, -1943.79664572349884, 4743.76565040243084, -11673.46399116716364, 28926.9553054354509)
    else
        # TODO: crossover points and degrees could be tuned more
        return x ≤ 0.6 ? (x ≤ 0.053 ? (x ≤ 4.4e-3 ? @Ei_taylor64(x,4) :
                                                    @Ei_taylor64(x,8)) :
                                       @Ei_taylor64(x,15)) :
                          @Ei_taylor64(x,37)
    end
end

expinti(z::Union{T,Rational{T}}) where {T<:Integer} = expinti(float(z))
expinti(z::Float32) = Float32(expinti(Float64(z)))
expinti(z::Float16) = Float16(expinti(Float64(z)))

function expinti(x::BigFloat)
    z = BigFloat()
    ccall((:mpfr_eint, :libmpfr), Int32, (Ref{BigFloat}, Ref{BigFloat}, Base.MPFR.MPFRRoundingMode), z, x, Base.MPFR.ROUNDING_MODE[])
    return z
end
