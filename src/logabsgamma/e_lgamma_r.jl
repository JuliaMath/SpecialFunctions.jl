#=
/* @(#)e_lgamma_r.c 1.3 95/01/18 */
/*
* ====================================================
* Copyright (C) 1993 by Sun Microsystems, Inc. All rights reserved.
*
* Developed at SunSoft, a Sun Microsystems, Inc. business.
* Permission to use, copy, modify, and distribute this
* software is freely granted, provided that this notice
* is preserved.
* ====================================================
*
*/
=#

#=

/* __ieee754_lgamma_r(x, signgamp)
 * Reentrant version of the logarithm of the Gamma function
 * with user provide pointer for the sign of Gamma(x).
 *
 * Method:
 *   1. Argument Reduction for 0 < x <= 8
 * 	Since gamma(1+s)=s*gamma(s), for x in [0,8], we may
 * 	reduce x to a number in [1.5,2.5] by
 * 		lgamma(1+s) = log(s) + lgamma(s)
 *	for example,
 *		lgamma(7.3) = log(6.3) + lgamma(6.3)
 *			    = log(6.3*5.3) + lgamma(5.3)
 *			    = log(6.3*5.3*4.3*3.3*2.3) + lgamma(2.3)
 *   2. Polynomial approximation of lgamma around its
 *	minimun ymin=1.461632144968362245 to maintain monotonicity.
 *	On [ymin-0.23, ymin+0.27] (i.e., [1.23164,1.73163]), use
 *		Let z = x-ymin;
 *		lgamma(x) = -1.214862905358496078218 + z^2*poly(z)
 *	where
 *		poly(z) is a 14 degree polynomial.
 *   2. Rational approximation in the primary interval [2,3]
 *	We use the following approximation:
 *		s = x-2.0;
 *		lgamma(x) = 0.5*s + s*P(s)/Q(s)
 *	with accuracy
 *		|P/Q - (lgamma(x)-0.5s)| < 2**-61.71
 *	Our algorithms are based on the following observation
 *
 *                             zeta(2)-1    2    zeta(3)-1    3
 * lgamma(2+s) = s*(1-Euler) + --------- * s  -  --------- * s  + ...
 *                                 2                 3
 *
 *	where Euler = 0.5771... is the Euler constant, which is very
 *	close to 0.5.
 *
 *   3. For x>=8, we have
 *	lgamma(x)~(x-0.5)log(x)-x+0.5*log(2pi)+1/(12x)-1/(360x**3)+....
 *	(better formula:
 *	   lgamma(x)~(x-0.5)*(log(x)-1)-.5*(log(2pi)-1) + ...)
 *	Let z = 1/x, then we approximation
 *		f(z) = lgamma(x) - (x-0.5)(log(x)-1)
 *	by
 *	  			    3       5             11
 *		w = w0 + w1*z + w2*z  + w3*z  + ... + w6*z
 *	where
 *		|w - f(z)| < 2**-58.74
 *
 *   4. For negative x, since (G is gamma function)
 *		-x*G(-x)*G(x) = pi/sin(pi*x),
 * 	we have
 * 		G(x) = pi/(sin(pi*x)*(-x)*G(-x))
 *	since G(-x) is positive, sign(G(x)) = sign(sin(pi*x)) for x<0
 *	Hence, for x<0, signgam = sign(sin(pi*x)) and
 *		lgamma(x) = log(|Gamma(x)|)
 *			  = log(pi/(|x*sin(pi*x)|)) - lgamma(-x);
 *	Note: one should avoid compute pi*(-x) directly in the
 *	      computation of sin(pi*(-x)).
 *
 *   5. Special Cases
 *		lgamma(2+s) ~ s*(1-Euler) for tiny s
 *		lgamma(1) = lgamma(2) = 0
 *		lgamma(x) ~ -log(|x|) for tiny x
 *		lgamma(0) = lgamma(neg.integer) = inf and raise divide-by-zero
 *		lgamma(inf) = inf
 *		lgamma(-inf) = inf (bug for bug compatible with C99!?)
 *
 */

=#

# Matches OpenLibm behavior (except commented out |x|≥2^52 early exit)
function _logabsgamma(x::Float64)
    ux = reinterpret(UInt64, x)
    hx = ux >>> 32 % Int32
    lx = ux % UInt32

    #= purge off +-inf, NaN, +-0, tiny and negative arguments =#
    signgam = 1
    isneg = hx < Int32(0)
    ix = hx & 0x7fffffff
    ix ≥ 0x7ff00000 && return x * x, signgam
    ix | lx == 0x00000000 && return Inf, signgam
    if ix < 0x3b900000 #= |x|<2**-70, return -log(|x|) =#
        if isneg
            signgam = -1
            return -log(-x), signgam
        else
            return -log(x), signgam
        end
    end
    if isneg
        # ix ≥ 0x43300000 && return Inf, signgam #= |x|>=2**52, must be -integer =#
        t = sinpi(x)
        iszero(t) && return Inf, signgam #= -integer =#
        nadj = logπ - log(abs(t * x))
        if t < 0.0; signgam = -1; end
        x = -x
    end
    if ix < 0x40000000     #= x < 2.0 =#
        i = round(x, RoundToZero)
        f = x - i
        if f == 0.0 #= purge off 1; 2 handled by x < 8.0 branch =#
            return 0.0, signgam
        elseif i == 1.0
            r = 0.0
            c = 1.0
        else
            r = -log(x)
            c = 0.0
        end
        if f ≥ 0.7315998077392578
            y = 1.0 + c - x
            z = y * y
            p1 = @evalpoly(z, 7.72156649015328655494e-02, 6.73523010531292681824e-02, 7.38555086081402883957e-03, 1.19270763183362067845e-03, 2.20862790713908385557e-04, 2.52144565451257326939e-05)
            p2 = z * @evalpoly(z, 3.22467033424113591611e-01, 2.05808084325167332806e-02, 2.89051383673415629091e-03, 5.10069792153511336608e-04, 1.08011567247583939954e-04, 4.48640949618915160150e-05)
            p = muladd(p1, y, p2)
            r += muladd(y, -0.5, p)
        elseif f ≥ 0.2316399812698364 # or, the lb? 0.2316322326660156
            y = x - 0.46163214496836225 - c
            z = y * y
            w = z * y
            p1 = @evalpoly(w, 4.83836122723810047042e-01, -3.27885410759859649565e-02, 6.10053870246291332635e-03, -1.40346469989232843813e-03, 3.15632070903625950361e-04)
            p2 = @evalpoly(w, -1.47587722994593911752e-01, 1.79706750811820387126e-02, -3.68452016781138256760e-03, 8.81081882437654011382e-04, -3.12754168375120860518e-04)
            p3 = @evalpoly(w, 6.46249402391333854778e-02, -1.03142241298341437450e-02, 2.25964780900612472250e-03, -5.38595305356740546715e-04, 3.35529192635519073543e-04)
            p = muladd(z, p1, -muladd(w, -muladd(p3, y, p2), -3.63867699703950536541e-18))
            r += p - 1.21486290535849611461e-1
        else
            y = x - c
            p1 = y * @evalpoly(y, -7.72156649015328655494e-02, 6.32827064025093366517e-01, 1.45492250137234768737, 9.77717527963372745603e-01, 2.28963728064692451092e-01, 1.33810918536787660377e-02)
            p2 = @evalpoly(y, 1.0, 2.45597793713041134822, 2.12848976379893395361, 7.69285150456672783825e-01, 1.04222645593369134254e-01, 3.21709242282423911810e-03)
		    r += muladd(y, -0.5, p1 / p2)
        end
    elseif ix < 0x40200000              #= x < 8.0 =#
        i = round(x, RoundToZero)
        y = x - i
	    z = 1.0
        p = 0.0
        u = x
        while u ≥ 3.0
            p -= 1.0
            u = x + p
            z *= u
        end
        p = y * @evalpoly(y, -7.72156649015328655494e-2, 2.14982415960608852501e-1, 3.25778796408930981787e-1, 1.46350472652464452805e-1, 2.66422703033638609560e-2, 1.84028451407337715652e-3, 3.19475326584100867617e-5)
        q = @evalpoly(y, 1.0, 1.39200533467621045958, 7.21935547567138069525e-1, 1.71933865632803078993e-1, 1.86459191715652901344e-2, 7.77942496381893596434e-4, 7.32668430744625636189e-6)
        r = log(z) + muladd(0.5, y, p / q)
    elseif ix < 0x43900000              #= 8.0 ≤ x < 2^58 =#
        z = 1.0 / x
        y = z * z
        w = muladd(z, @evalpoly(y, 8.33333333333329678849e-2, -2.77777777728775536470e-3, 7.93650558643019558500e-4, -5.95187557450339963135e-4, 8.36339918996282139126e-4, -1.63092934096575273989e-3), 4.18938533204672725052e-1)
	    r = muladd(x - 0.5, log(x) - 1.0, w)
    else #= 2^58 ≤ x ≤ Inf =#
        r = muladd(x, log(x), -x)
    end
    if isneg
        r = nadj - r
    end
    return r, signgam
end
