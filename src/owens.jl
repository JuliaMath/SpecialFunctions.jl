function owens_t(h::T, a::T; rtol=1e-12, atol=1e-15) where {T<:Real}
    isnan(h) && return h
    isnan(a) && return a

    (iszero(a) || isinf(h)) && return zero(T)

    s = sign(a)
    θmax = atan(abs(a))

    # Exact result when h == 0
    if iszero(h)
        return s * θmax / (2π)
    end

    h2_over_2 = 0.5 * h^2
    integrand(θ) = exp(-h2_over_2 / cos(θ)^2)

    integral,_ = QuadGK.quadgk(
        integrand,
        0.0,
        θmax,
        rtol=rtol,
        atol=atol,
    )

    return s * integral / (2π)
end
