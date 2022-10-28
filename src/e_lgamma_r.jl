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

const a0  =  7.72156649015328655494e-02 #= 0x3FB3C467, 0xE37DB0C8 =#
const a1  =  3.22467033424113591611e-01 #= 0x3FD4A34C, 0xC4A60FAD =#
const a2  =  6.73523010531292681824e-02 #= 0x3FB13E00, 0x1A5562A7 =#
const a3  =  2.05808084325167332806e-02 #= 0x3F951322, 0xAC92547B =#
const a4  =  7.38555086081402883957e-03 #= 0x3F7E404F, 0xB68FEFE8 =#
const a5  =  2.89051383673415629091e-03 #= 0x3F67ADD8, 0xCCB7926B =#
const a6  =  1.19270763183362067845e-03 #= 0x3F538A94, 0x116F3F5D =#
const a7  =  5.10069792153511336608e-04 #= 0x3F40B6C6, 0x89B99C00 =#
const a8  =  2.20862790713908385557e-04 #= 0x3F2CF2EC, 0xED10E54D =#
const a9  =  1.08011567247583939954e-04 #= 0x3F1C5088, 0x987DFB07 =#
const a10 =  2.52144565451257326939e-05 #= 0x3EFA7074, 0x428CFA52 =#
const a11 =  4.48640949618915160150e-05 #= 0x3F07858E, 0x90A45837 =#
const tc  =  1.46163214496836224576e+00 #= 0x3FF762D8, 0x6356BE3F =#
const tf  = -1.21486290535849611461e-01 #= 0xBFBF19B9, 0xBCC38A42 =#
#= tt = -(tail of tf) =#
const tt  = -3.63867699703950536541e-18 #= 0xBC50C7CA, 0xA48A971F =#
const t0  =  4.83836122723810047042e-01 #= 0x3FDEF72B, 0xC8EE38A2 =#
const t1  = -1.47587722994593911752e-01 #= 0xBFC2E427, 0x8DC6C509 =#
const t2  =  6.46249402391333854778e-02 #= 0x3FB08B42, 0x94D5419B =#
const t3  = -3.27885410759859649565e-02 #= 0xBFA0C9A8, 0xDF35B713 =#
const t4  =  1.79706750811820387126e-02 #= 0x3F9266E7, 0x970AF9EC =#
const t5  = -1.03142241298341437450e-02 #= 0xBF851F9F, 0xBA91EC6A =#
const t6  =  6.10053870246291332635e-03 #= 0x3F78FCE0, 0xE370E344 =#
const t7  = -3.68452016781138256760e-03 #= 0xBF6E2EFF, 0xB3E914D7 =#
const t8  =  2.25964780900612472250e-03 #= 0x3F6282D3, 0x2E15C915 =#
const t9  = -1.40346469989232843813e-03 #= 0xBF56FE8E, 0xBF2D1AF1 =#
const t10 =  8.81081882437654011382e-04 #= 0x3F4CDF0C, 0xEF61A8E9 =#
const t11 = -5.38595305356740546715e-04 #= 0xBF41A610, 0x9C73E0EC =#
const t12 =  3.15632070903625950361e-04 #= 0x3F34AF6D, 0x6C0EBBF7 =#
const t13 = -3.12754168375120860518e-04 #= 0xBF347F24, 0xECC38C38 =#
const t14 =  3.35529192635519073543e-04 #= 0x3F35FD3E, 0xE8C2D3F4 =#
const u0  = -7.72156649015328655494e-02 #= 0xBFB3C467, 0xE37DB0C8 =#
const u1  =  6.32827064025093366517e-01 #= 0x3FE4401E, 0x8B005DFF =#
const u2  =  1.45492250137234768737e+00 #= 0x3FF7475C, 0xD119BD6F =#
const u3  =  9.77717527963372745603e-01 #= 0x3FEF4976, 0x44EA8450 =#
const u4  =  2.28963728064692451092e-01 #= 0x3FCD4EAE, 0xF6010924 =#
const u5  =  1.33810918536787660377e-02 #= 0x3F8B678B, 0xBF2BAB09 =#
const v1  =  2.45597793713041134822e+00 #= 0x4003A5D7, 0xC2BD619C =#
const v2  =  2.12848976379893395361e+00 #= 0x40010725, 0xA42B18F5 =#
const v3  =  7.69285150456672783825e-01 #= 0x3FE89DFB, 0xE45050AF =#
const v4  =  1.04222645593369134254e-01 #= 0x3FBAAE55, 0xD6537C88 =#
const v5  =  3.21709242282423911810e-03 #= 0x3F6A5ABB, 0x57D0CF61 =#
const s0  = -7.72156649015328655494e-02 #= 0xBFB3C467, 0xE37DB0C8 =#
const s1  =  2.14982415960608852501e-01 #= 0x3FCB848B, 0x36E20878 =#
const s2  =  3.25778796408930981787e-01 #= 0x3FD4D98F, 0x4F139F59 =#
const s3  =  1.46350472652464452805e-01 #= 0x3FC2BB9C, 0xBEE5F2F7 =#
const s4  =  2.66422703033638609560e-02 #= 0x3F9B481C, 0x7E939961 =#
const s5  =  1.84028451407337715652e-03 #= 0x3F5E26B6, 0x7368F239 =#
const s6  =  3.19475326584100867617e-05 #= 0x3F00BFEC, 0xDD17E945 =#
const r1  =  1.39200533467621045958e+00 #= 0x3FF645A7, 0x62C4AB74 =#
const r2  =  7.21935547567138069525e-01 #= 0x3FE71A18, 0x93D3DCDC =#
const r3  =  1.71933865632803078993e-01 #= 0x3FC601ED, 0xCCFBDF27 =#
const r4  =  1.86459191715652901344e-02 #= 0x3F9317EA, 0x742ED475 =#
const r5  =  7.77942496381893596434e-04 #= 0x3F497DDA, 0xCA41A95B =#
const r6  =  7.32668430744625636189e-06 #= 0x3EDEBAF7, 0xA5B38140 =#
const w0  =  4.18938533204672725052e-01 #= 0x3FDACFE3, 0x90C97D69 =#
const w1  =  8.33333333333329678849e-02 #= 0x3FB55555, 0x5555553B =#
const w2  = -2.77777777728775536470e-03 #= 0xBF66C16C, 0x16B02E5C =#
const w3  =  7.93650558643019558500e-04 #= 0x3F4A019F, 0x98CF38B6 =#
const w4  = -5.95187557450339963135e-04 #= 0xBF4380CB, 0x8C0FE741 =#
const w5  =  8.36339918996282139126e-04 #= 0x3F4B67BA, 0x4CDAD5D1 =#
const w6  = -1.63092934096575273989e-03 #= 0xBF5AB89D, 0x0B9E43E4 =#

# Matches OpenLibm behavior exactly, including return of sign
function _lgamma_r(x::Float64)
    u = reinterpret(UInt64, x)
    hx = (u >>> 32) % Int32
    lx = u % Int32

    #= purge off +-inf, NaN, +-0, tiny and negative arguments =#
    signgamp = Int32(1)
    ix = signed(hx & 0x7fffffff)
    ix ≥ 0x7ff00000 && return x * x, signgamp
    ix | lx == 0 && return 1.0 / 0.0, signgamp
    if ix < 0x3b900000 #= |x|<2**-70, return -log(|x|) =#
        if hx < 0
            signgamp = Int32(-1)
            return -log(-x), signgamp
        else
            return -log(x), signgamp
        end
    end
    if hx < 0
        ix ≥ 0x43300000 && return 1.0 / 0.0, signgamp #= |x|>=2**52, must be -integer =#
        t = sinpi(x)
        t == 0.0 && return 1.0 / 0.0, signgamp #= -integer =#
        nadj = log(π / abs(t * x))
        if t < 0.0; signgamp = Int32(-1); end
        x = -x
    end

    #= purge off 1 and 2 =#
    if ((ix - 0x3ff00000) | lx) == 0 || ((ix - 0x40000000) | lx) == 0
        r = 0.0
        #= for x < 2.0 =#
    elseif ix < 0x40000000
        if ix ≤ 0x3feccccc #= lgamma(x) = lgamma(x+1)-log(x) =#
            r = -log(x)
            if ix ≥ 0x3FE76944
                y = 1.0 - x
                i = Int8(0)
            elseif ix ≥ 0x3FCDA661
                y = x - (tc - 1.0)
                i = Int8(1)
            else
                y = x
                i = Int8(2)
            end
        else
            r = 0.0
            if ix ≥ 0x3FFBB4C3 #= [1.7316,2] =#
                y = 2.0 - x
                i = Int8(0)
            elseif ix ≥ 0x3FF3B4C4 #= [1.23,1.73] =#
                y = x - tc
                i = Int8(1)
            else
                y = x - 1.0
                i = Int8(2)
            end
        end
        if i == Int8(0)
            z = y*y;
		    p1 = a0+z*(a2+z*(a4+z*(a6+z*(a8+z*a10))));
		    p2 = z*(a1+z*(a3+z*(a5+z*(a7+z*(a9+z*a11)))));
		    p  = y*p1+p2;
		    r  += (p-0.5*y);
        elseif i == Int8(1)
            z = y*y;
		    w = z*y;
		    p1 = t0+w*(t3+w*(t6+w*(t9 +w*t12)));	#= parallel comp =#
		    p2 = t1+w*(t4+w*(t7+w*(t10+w*t13)));
		    p3 = t2+w*(t5+w*(t8+w*(t11+w*t14)));
		    p  = z*p1-(tt-w*(p2+y*p3));
		    r += (tf + p)
        elseif i == Int8(2)
            p1 = y*(u0+y*(u1+y*(u2+y*(u3+y*(u4+y*u5)))));
		    p2 = 1.0+y*(v1+y*(v2+y*(v3+y*(v4+y*v5))));
		    r += (-0.5*y + p1/p2);
        end
    elseif ix < 0x40200000                #= x < 8.0 =#
        i = Base.unsafe_trunc(Int8, x)
        y = x - float(i)
        # If performed here, performance is 2x worse; hence, move it below.
        # p = y*(s0+y*(s1+y*(s2+y*(s3+y*(s4+y*(s5+y*s6))))));
	    # q = 1.0+y*(r1+y*(r2+y*(r3+y*(r4+y*(r5+y*r6)))));
	    # r = 0.5*y+p/q;
	    z = 1.0;	#= lgamma(1+s) = log(s) + lgamma(s) =#
        if i == Int8(7)
            z *= (y + 6.0)
            @goto case6
        elseif i == Int8(6)
            @label case6
            z *= (y + 5.0)
            @goto case5
        elseif i == Int8(5)
            @label case5
            z *= (y + 4.0)
            @goto case4
        elseif i == Int8(4)
            @label case4
            z *= (y + 3.0)
            @goto case3
        elseif i == Int8(3)
            @label case3
            z *= (y + 2.0)
        end
        # r += log(z)
        p = y*(s0+y*(s1+y*(s2+y*(s3+y*(s4+y*(s5+y*s6))))));
        q = 1.0+y*(r1+y*(r2+y*(r3+y*(r4+y*(r5+y*r6)))));
        r = log(z) + 0.5*y+p/q;
        #= 8.0 ≤ x < 2^58 =#
    elseif ix < 0x43900000
        t = log(x)
        z = 1.0 / x
        y = z * z
        w = w0+z*(w1+y*(w2+y*(w3+y*(w4+y*(w5+y*w6)))));
	    r = (x-0.5)*(t-1.0)+w;
    else
        #= 2^58 ≤ x ≤ Inf =#
        r = x * (log(x) - 1.0)
    end
    if hx < 0
        r = nadj - r
    end
    return r, signgamp
end

# Deviates from OpenLibm: throws instead of returning negative sign; approximately 25% faster
# when sign is not needed in subsequent computations.
function _loggamma_r(x::Float64)
    u = reinterpret(UInt64, x)
    hx = (u >>> 32) % Int32
    lx = u % Int32

    #= purge off +-inf, NaN, +-0, tiny and negative arguments =#
    ix = signed(hx & 0x7fffffff)
    ix ≥ 0x7ff00000 && return x * x
    ix | lx == 0 && return 1.0 / 0.0
    if ix < 0x3b900000 #= |x|<2**-70, return -log(|x|) =#
        if hx < 0
            # return -log(-x)
            throw(DomainError(x, "`gamma(x)` must be non-negative"))
        else
            return -log(x)
        end
    end
    if hx < 0
        ix ≥ 0x43300000 && return 1.0 / 0.0 #= |x|>=2**52, must be -integer =#
        t = sinpi(x)
        t == 0.0 && return 1.0 / 0.0 #= -integer =#
        nadj = log(π / abs(t * x))
        t < 0.0 && throw(DomainError(x, "`gamma(x)` must be non-negative"))
        x = -x
    end

    #= purge off 1 and 2 =#
    if ((ix - 0x3ff00000) | lx) == 0 || ((ix - 0x40000000) | lx) == 0
        r = 0.0
        #= for x < 2.0 =#
    elseif ix < 0x40000000
        if ix ≤ 0x3feccccc #= lgamma(x) = lgamma(x+1)-log(x) =#
            r = -log(x)
            if ix ≥ 0x3FE76944
                y = 1.0 - x
                i = Int8(0)
            elseif ix ≥ 0x3FCDA661
                y = x - (tc - 1.0)
                i = Int8(1)
            else
                y = x
                i = Int8(2)
            end
        else
            r = 0.0
            if ix ≥ 0x3FFBB4C3 #= [1.7316,2] =#
                y = 2.0 - x
                i = Int8(0)
            elseif ix ≥ 0x3FF3B4C4 #= [1.23,1.73] =#
                y = x - tc
                i = Int8(1)
            else
                y = x - 1.0
                i = Int8(2)
            end
        end
        if i == Int8(0)
            z = y*y;
		    p1 = a0+z*(a2+z*(a4+z*(a6+z*(a8+z*a10))));
		    p2 = z*(a1+z*(a3+z*(a5+z*(a7+z*(a9+z*a11)))));
		    p  = y*p1+p2;
		    r  += (p-0.5*y);
        elseif i == Int8(1)
            z = y*y;
		    w = z*y;
		    p1 = t0+w*(t3+w*(t6+w*(t9 +w*t12)));	#= parallel comp =#
		    p2 = t1+w*(t4+w*(t7+w*(t10+w*t13)));
		    p3 = t2+w*(t5+w*(t8+w*(t11+w*t14)));
		    p  = z*p1-(tt-w*(p2+y*p3));
		    r += (tf + p)
        elseif i == Int8(2)
            p1 = y*(u0+y*(u1+y*(u2+y*(u3+y*(u4+y*u5)))));
		    p2 = 1.0+y*(v1+y*(v2+y*(v3+y*(v4+y*v5))));
		    r += (-0.5*y + p1/p2);
        end
    elseif ix < 0x40200000                #= x < 8.0 =#
        i = Base.unsafe_trunc(Int8, x)
        y = x - float(i)
        # If performed here, performance is 2x worse; hence, move it below.
        # p = y*(s0+y*(s1+y*(s2+y*(s3+y*(s4+y*(s5+y*s6))))));
	    # q = 1.0+y*(r1+y*(r2+y*(r3+y*(r4+y*(r5+y*r6)))));
	    # r = 0.5*y+p/q;
	    z = 1.0;	#= lgamma(1+s) = log(s) + lgamma(s) =#
        if i == Int8(7)
            z *= (y + 6.0)
            @goto case6
        elseif i == Int8(6)
            @label case6
            z *= (y + 5.0)
            @goto case5
        elseif i == Int8(5)
            @label case5
            z *= (y + 4.0)
            @goto case4
        elseif i == Int8(4)
            @label case4
            z *= (y + 3.0)
            @goto case3
        elseif i == Int8(3)
            @label case3
            z *= (y + 2.0)
        end
        # r += log(z)
        p = y*(s0+y*(s1+y*(s2+y*(s3+y*(s4+y*(s5+y*s6))))));
        q = 1.0+y*(r1+y*(r2+y*(r3+y*(r4+y*(r5+y*r6)))));
        r = log(z) + 0.5*y+p/q;
        #= 8.0 ≤ x < 2^58 =#
    elseif ix < 0x43900000
        t = log(x)
        z = 1.0 / x
        y = z * z
        w = w0+z*(w1+y*(w2+y*(w3+y*(w4+y*(w5+y*w6)))));
	    r = (x-0.5)*(t-1.0)+w;
    else
        #= 2^58 ≤ x ≤ Inf =#
        r = x * (log(x) - 1.0)
    end
    if hx < 0
        r = nadj - r
    end
    return r
end
