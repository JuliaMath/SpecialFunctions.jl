#=
/* e_lgammaf_r.c -- float version of e_lgamma_r.c.
* Conversion to float by Ian Lance Taylor, Cygnus Support, ian@cygnus.com.
*/

/*
* ====================================================
* Copyright (C) 1993 by Sun Microsystems, Inc. All rights reserved.
*
* Developed at SunPro, a Sun Microsystems, Inc. business.
* Permission to use, copy, modify, and distribute this
* software is freely granted, provided that this notice
* is preserved.
* ====================================================
*/
=#

# Matches OpenLibm behavior exactly, including return of sign
function _logabsgamma(x::Float32)
    hx = reinterpret(Int32, x)

    #= purge off +-inf, NaN, +-0, tiny and negative arguments =#
    signgam = 1
    isneg = hx < Int32(0)
    ix = hx & 0x7fffffff
    ix ≥ 0x7f800000 && return x * x, signgam
    ix == 0x00000000 && return Inf32, signgam
    if ix < 0x35000000 #= |x|<2**-21, return -log(|x|) =#
        if isneg
            signgam = -1
            return -log(-x), signgam
        else
            return -log(x), signgam
        end
    end
    if isneg
        # ix ≥ 0x4b000000 && return Inf32, signgam #= |x|>=2**23, must be -integer =#
        t = sinpi(x)
        t == 0.0f0 && return Inf32, signgam #= -integer =#
        nadj = logπ - log(abs(t * x))
        if t < 0.0f0; signgam = -1; end
        x = -x
    end

    if ix < 0x40000000 #= x < 2.0 =#
        i = round(x, RoundToZero)
        f = x - i
        if f == 0.0f0 #= purge off 1; 2 handled by x < 8.0 branch =#
            return 0.0f0, signgam
        elseif i == 1.0f0
            r = 0.0f0
            c = 1.0f0
        else
            r = -log(x)
            c = 0.0f0
        end
        if f ≥ 0.7315998f0
            y = 1.0f0 + c - x
            z = y * y
            p1 = @evalpoly(z, 7.7215664089f-2, 6.7352302372f-2, 7.3855509982f-3, 1.1927076848f-3, 2.2086278477f-4, 2.5214456400f-5)
            p2 = z * @evalpoly(z, 3.2246702909f-1, 2.0580807701f-2, 2.8905137442f-3, 5.1006977446f-4, 1.0801156895f-4, 4.4864096708f-5)
            p = muladd(p1, y, p2)
            r += muladd(y, -0.5f0, p)
        elseif f ≥ 0.23163998f0 # or, the lb? 0.2316322f0
            y = x - 0.46163213f0 - c
            z = y * y
            w = z * y
            p1 = @evalpoly(w, 4.8383611441f-1, -3.2788541168f-2, 6.1005386524f-3, -1.4034647029f-3, 3.1563205994f-4)
            p2 = @evalpoly(w, -1.4758771658f-1, 1.7970675603f-2, -3.6845202558f-3, 8.8108185446f-4, -3.1275415677f-4)
            p3 = @evalpoly(w, 6.4624942839f-2, -1.0314224288f-2, 2.2596477065f-3, -5.3859531181f-4, 3.3552918467f-4)
            p = muladd(z, p1, -muladd(w, -muladd(p3, y, p2), 6.6971006518f-9))
            r += p - 1.2148628384f-1
        else
            y = x - c
            p1 = y * @evalpoly(y, -7.7215664089f-2, 6.3282704353f-1, 1.4549225569f0, 9.7771751881f-1, 2.2896373272f-1, 1.3381091878f-2)
            p2 = @evalpoly(y, 1.0f0, 2.4559779167f0, 2.1284897327f0, 7.6928514242f-1, 1.0422264785f-1, 3.2170924824f-3)
            r += muladd(y, -0.5f0, p1 / p2)
        end
    elseif ix < 0x41000000 #= x < 8.0 =#
        i = round(x, RoundToZero)
        y = x - i
        z = 1.0f0
        p = 0.0f0
        u = x
        while u ≥ 3.0f0
            p -= 1.0f0
            u = x + p
            z *= u
        end
        p = y * @evalpoly(y, -7.7215664089f-2, 2.1498242021f-1, 3.2577878237f-1, 1.4635047317f-1, 2.6642270386f-2, 1.8402845599f-3, 3.1947532989f-5)
        q = @evalpoly(y, 1.0f0, 1.3920053244f0, 7.2193557024f-1, 1.7193385959f-1, 1.8645919859f-2, 7.7794247773f-4, 7.3266842264f-6)
	    r = log(z) + muladd(0.5f0, y, p / q)
    elseif ix < 0x5c800000 #= 8.0 ≤ x < 2^58 =#
        z = 1.0f0 / x
        y = z * z
        w = muladd(z, @evalpoly(y, 8.3333335817f-2, -2.7777778450f-3, 7.9365057172f-4, -5.9518753551f-4, 8.3633989561f-4, -1.6309292987f-3), 4.1893854737f-1)
        r = muladd(x - 0.5f0, log(x) - 1.0f0, w)
    else
        #= 2^58 ≤ x ≤ Inf =#
        r = muladd(x, log(x), -x)
    end
    if isneg
        r = nadj - r
    end
    return r, signgam
end
