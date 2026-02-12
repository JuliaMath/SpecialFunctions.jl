# Implementation of Wright's generalized Bessel function Φ, see
# https://dlmf.nist.gov/10.46.E1
#
# Based on SciPy's Cython implementation, licensed under BSD-3-Clause:
# Copyright: Christian Lorentzen
# Copyright (c) 2001-2002 Enthought, Inc. 2003, SciPy Developers.
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions
# are met:
#
# 1. Redistributions of source code must retain the above copyright
#    notice, this list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above
#    copyright notice, this list of conditions and the following
#    disclaimer in the documentation and/or other materials provided
#    with the distribution.
#
# 3. Neither the name of the copyright holder nor the names of its
#    contributors may be used to endorse or promote products derived
#    from this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
# A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
# OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
# LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
# DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
# THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

# Implementation Overview:
#
# First, different functions are implemented valid for certain domains of the
# three arguments.
# Finally they are put together in `wrightbessel``.

# rgamma_zero: smallest value x for which rgamma(x) == 0 as x gets large
const rgamma_zero = 178.47241115886637
# exp_inf: smallest value x for which exp(x) == inf
const exp_inf = 709.7827128933841

# Compute reciprocal gamma via loggamma
@inline rgamma(y::Real) = exp(-loggamma(y))

# Compute exp(x) / gamma(y) safely via loggamma
@inline function _exp_rgamma(x::Real, y::Real)
    # Compute exp(x) / Γ(y) = exp(x - loggamma(y)) to avoid overflow
    return exp(x - loggamma(y))
end


# 1. Taylor series expansion in x=0, for x <= 1
# Phi(a, b, x) = sum_k x^k / k! / Gamma(a*k+b)
# Note that every term, and therefore also Phi(a, b, x), is monotone
# decreasing with increasing a or b.
function _wb_series(a::Float64, b::Float64, x::Float64, nstart::Int, nstop::Int)
    xk_k = x^nstart * rgamma(nstart + 1)  # x^k/k!
    res = xk_k * rgamma(nstart * a + b)
    # term k=nstart+1 , +2 +3, ...
    if nstop > nstart && a > 0.0
        # series expansion until term k such that a*k+b <= rgamma_zero
        k_max = floor(Int, (rgamma_zero - b) / a)
        if nstop > k_max
            nstop = k_max
        end
        for k in (nstart + 1):(nstop - 1)
            xk_k *= x / k
            res += xk_k * rgamma(a * k + b)
        end
    end
    return res
end


# 2. Taylor series expansion in x=0, for large a.
# Phi(a, b, x) = sum_k x^k / k! / Gamma(a*k+b)
# Use Stirling formula to find k=k_max, the maximum term.
# Then use n terms of Taylor series around k_max.
function _wb_large_a(a::Float64, b::Float64, x::Float64, n::Int)
    k_max = floor(Int, (a^(-a) * x)^(1 / (1 + a)))
    n_start = k_max - (n ÷ 2)
    if n_start < 0
        n_start = 0
    end

    res = 0.0
    lnx = log(x)
    for k in n_start:(n_start + n - 1)
        res += exp(k * lnx - loggamma(k + 1.0) - loggamma(a * k + b))
    end
    return res
end


# 3. Taylor series in a=0 up to order 5, for tiny a and not too large x
#
# Phi(a, b, x) = exp(x)/Gamma(b)
#                * (1 - a*x * Psi(b) + a^2/2*x*(1+x) * (Psi(b)^2 - Psi'(b)
#                   + ... )
#                + O(a^6))
#
# where Psi is the digamma function.
#
# Parameter order takes effect only when b > 1e-3 and 2 <= order <= 5,
# otherwise it defaults to 2 or, if b <= 1e-3, to 5. The lower order is, the
# less polygamma functions have to be computed.
#
# For small b, i.e. b <= 1e-3, cancellation of poles of digamma(b)/Gamma(b)
# and polygamma needs to be carried out => series expansion in a=0 to order 5
# and in b=0 to order 4.
function _wb_small_a(a::Float64, b::Float64, x::Float64, order::Int)
    # A: coefficients of a^k  (1, -x * Psi(b), ...)
    # B: powers of b^k/k! or terms in polygamma functions
    # C: coefficients of a^k1 * b^k2
    # X: polynomials in x
    if b <= 1e-3
        # Series expansion of both a and b up to order 5:
        # γ = Euler Gamma aka Euler Mascheroni constant
        # ζ(3) = zeta(3)
        C1 = 1.0000000000000000
        # C2 = 2*γ
        C2 = 1.1544313298030657
        # C3 = 3*γ^2 - π^2/2
        C3 = -3.9352684291215233
        # C4 = 4*γ^3 - 2*γ*π^2 + 8*ζ(3)
        C4 = -1.0080632408182857
        # C5 = 5*γ^4 - 5*γ^2*π^2 + 40*γ*ζ(3) + π^4/12
        C5 = 19.984633365874979
        X1 = 1.0
        X2 = x
        X3 = x * (x + 1.0)
        X4 = x * (x * (x + 3.0) + 1.0)
        X5 = x * (x * (x * (x + 6.0) + 7.0) + 1.0)
        X6 = x * (x * (x * (x * (x + 10.0) + 25.0) + 15.0) + 1.0)
        B = cumprod(ntuple(i -> b / i, 4))
        A = (
            rgamma(b),
            X2         * (C1 + C2 * b + C3 * B[2] + C4 * B[3] + C5 * B[4]),
            X3 / 2.0   *      (C2     + C3 * b    + C4 * B[2] + C5 * B[3]),
            X4 / 6.0   *               (C3        + C4 * b    + C5 * B[2]),
            X5 / 24.0  *                           (C4        + C5 * b),
            X6 / 120.0 *                                       (C5)
        )
        res = exp(x) * evalpoly(a, A)
    else
        # Φ(a, b, x) = exp(x)/gamma(b) * sum(A[i] * X[i] * B[i], i=0..5)
        # A[n] = a^n/n!
        # But here, we repurpose A[n] = X[n] * B[n] / n!
        dg = digamma(b)
        # pg1 = polygamma(1, b)
        pg1 = zeta(2.0, b)
        if order <= 2
            res = 1.0 + a * x * (-dg + 0.5 * a * (1.0 + x) * (dg^2 - pg1))
        else
            order = min(order, 5)
            # pg2 = polygamma(2, b)
            pg2 = -2.0 * zeta(3.0, b)
            X1 = 1.0
            X2 = x
            X3 = x * (x + 1.0)
            X4 = x * (x * (x + 3.0) + 1.0)
            B1 = 1.0
            B2 = -dg
            B3 = dg^2 - pg1
            B4 = (-dg^2 + 3.0 * pg1) * dg - pg2
            A1 = 1.0
            A2 = X2 * B2
            A3 = A4 = A5 = A6 = NaN

            if order >= 2
                A3 = X3 * B3 / 2.0

                if order >= 3
                   A4 = X4 * B4 / 6.0

                    if order >= 4
                        # pg3 = polygamma(3, b)
                        pg3 =  6.0 * zeta(4.0, b)
                        X5 = x * (x * (x * (x + 6.0) + 7.0) + 1.0)
                        B5 = ((dg^2 - 6.0 * pg1) * dg + 4.0 * pg2) * dg + 3.0 * pg1^2 - pg3
                        A5 = X5 * B5 / 24.0

                        if order >= 5
                            # pg4 = polygamma(4, b)
                            pg4 = -24.0 * zeta(5.0, b)
                            X6 = x * (x * (x * (x * (x + 10.0) + 25.0) + 15.0) + 1.0)
                            B6 = ((((-dg^2 + 10.0 * pg1) * dg - 10.0 * pg2) * dg
                                    - 15.0 * pg1^2 + 5.0 * pg3) * dg
                                    + 10.0 * pg1 * pg2 - pg4)
                            A6 = X6 * B6 / 120.0
                        end
                    end
                end
            end
            res = evalpoly(a, (A1, A2, A3, A4, A5, A6)[1:(order+1)])
        end

        # res *= exp(x) * rgamma(b)
        res *= _exp_rgamma(x, b)
    end
    return res
end


# 4. Asymptotic expansion for large x up to order 8
#
# Phi(a, b, x) ~ Z^(1/2-b) * exp((1+a)/a * Z) * sum_k (-1)^k * C_k / Z^k
#
# with Z = (a*x)^(1/(1+a)).
function _wb_asymptotic(a::Float64, b::Float64, x::Float64)
    A = ntuple(i -> a^(i-1), 15)        # powers of a
    B = ntuple(i -> b^(i-1), 17)        # powers of b
    Ap1 = ntuple(i -> (1 + a)^(i-1), 9) # powers of (1+a)

    # C: coefficients of asymptotic series a_k

    C1 = 1.0 / sqrt(2.0 * π * Ap1[2])

    C2 = C1 / (24.0 * Ap1[2]) *
        ( (2.0*a + 1.0) * (2.0 + a) - 12.0 * b * (1.0 + a - b) )

    C3 = C1 / (1152.0 * Ap1[3]) * (
          144.0 * B[5]
        - 96.0 * B[4] * (5.0 * a + 1.0)
        + 24.0 * B[3] * (20.0 * A[3] + 5.0 * a - 4.0)
        - 24.0 * b * Ap1[2] * (6.0 * A[3] - 7.0 * a - 2.0)
        + (a + 2.0) * (2.0 * a + 1.0) * (2.0 * A[3] - 19.0 * a + 2.0)
    )

    C4 = C1 / (414720.0 * Ap1[4]) * (
          8640.0 * B[7]
        - 8640.0 * B[6] * (7.0 * a - 1.0)
        + 10800.0 * B[5] * (14.0 * A[3] - 7.0 * a - 2.0)
        - 1440.0 * B[4] * (112.0 * A[4] - 147.0 * A[3] - 63.0 * a + 8.0)
        + 180.0 * B[3] * (364.0 * A[5] - 1288.0 * A[4] - 567.0 * A[3] + 392.0 * a + 76.0)
        - 180.0 * b * Ap1[2] * (20.0 * A[5] - 516.0 * A[4] + 417.0 * A[3] + 172.0 * a - 12.0)
        - (a + 2.0) * (2.0 * a + 1.0) * (556.0 * A[5] + 1628.0 * A[4] - 9093.0 * A[3] + 1628.0 * a + 556.0)
    )

    C5 = C1 / (39813120.0 * Ap1[5]) * (
          103680.0 * B[9]
        - 414720.0 * B[8] * (3.0 * a - 1.0)
        + 725760.0 * B[7] * a * (8.0 * a - 7.0)
        - 48384.0 * B[6] * (274.0 * A[4] - 489.0 * A[3] + 39.0 * a + 26.0)
        + 30240.0 * B[5] * (500.0 * A[5] - 1740.0 * A[4] + 495.0 * A[3] + 340.0 * a - 12.0)
        - 2880.0 * B[4] * (
            2588.0 * A[6] - 19780.0 * A[5] + 14453.0 * A[4] + 9697.0 * A[3] - 1892.0 * a - 404.0
        )
        + 48.0 * B[3] * (
            11488.0 * A[7] - 547836.0 * A[6] + 1007484.0 * A[5] + 593353.0 * A[4]
            - 411276.0 * A[3] - 114396.0 * a + 4288.0
        )
        + 48.0 * b * Ap1[2] * (
            7784.0 * A[7] + 48180.0 * A[6] - 491202.0 * A[5] + 336347.0 * A[4]
            + 163734.0 * A[3] - 28908.0 * a - 5560.0
        )
        - (a + 2.0) * (2.0 * a + 1.0) * (
            4568.0 * A[7] - 226668.0 * A[6] - 465702.0 * A[5]
            + 2013479.0 * A[4] - 465702.0 * A[3] - 226668.0 * a + 4568.0
        )
    )

    C6 = C1 / (6688604160.0 * Ap1[6]) * (
          1741824.0 * B[11] - 2903040.0 * B[10] * (11.0 * a - 5.0)
        + 2177280.0 * B[9] * (110.0 * A[3] - 121.0 * a + 14.0)
        - 580608.0 * B[8] * (1628.0 * A[4] - 3333.0 * A[3] + 1023.0 * a + 52.0)
        + 169344.0 * B[7] * (12364.0 * A[5] - 43648.0 * A[4] + 26763.0 * A[3] + 1232.0 * a - 788.0)
        - 24192.0 * B[6] * (
            104852.0 * A[6] - 646624.0 * A[5] + 721391.0 * A[4] - 16841.0 * A[3]
            - 74096.0 * a + 148.0
        )
        + 2016.0 * B[5] * (
            710248.0 * A[7] - 8878716.0 * A[6] + 17928834.0 * A[5] - 3333407.0 * A[4]
            - 4339566.0 * A[3] + 287364.0 * a + 89128.0
        )
        - 1344.0 * B[4] * (
            87824.0 * A[8] - 7150220.0 * A[7] + 29202756.0 * A[6] - 15113527.0 * A[5]
            - 14223011.0 * A[4] + 3462492.0 * A[3] + 1137092.0 * a - 18896.0
        )
        - 84.0 * B[3] * (
            1690480.0 * A[9] + 14139136.0 * A[8] - 232575464.0 * A[7]
            + 296712592.0 * A[6] + 215856619.0 * A[5] - 152181392.0 * A[4]
            - 47718440.0 * A[3] + 5813632.0 * a + 943216.0
        )
        + 84.0 * b * Ap1[2] * (
            82224.0 * A[9] - 5628896.0 * A[8] - 26466520.0 * A[7]
            + 168779208.0 * A[6] - 104808005.0 * A[5] - 56259736.0 * A[4]
            + 15879912.0 * A[3] + 4020640.0 * a - 63952.0
        )
        + (a + 2.0) * (2.0 * a + 1.0) * (
            2622064.0 * A[9] + 12598624.0 * A[8]
            - 167685080.0 * A[7] - 302008904.0 * A[6]
            + 1.115235367e9 * A[5] - 302008904.0 * A[4]
            - 167685080.0 * A[3] + 12598624.0 * a + 2622064.0
        )
    )

    C7 = C1 / (4.8157949952e12 * Ap1[7]) * (
          104509440.0 * B[13] - 209018880.0 * B[12] * (13.0 * a - 7.0)
        + 574801920.0 * B[11] * (52.0 * A[3] - 65.0 * a + 12.0)
        - 63866880.0 * B[10] * (2834.0 * A[4] - 6279.0 * A[3] + 2769.0 * a - 134.0)
        + 23950080.0 * B[9] * (27404.0 * A[5] - 98228.0 * A[4] + 78663.0 * A[3] - 10868.0 * a - 1012.0)
        - 13685760.0 * B[8] * (
            105612.0 * A[6] - 599196.0 * A[5] + 791843.0 * A[4] - 224913.0 * A[3]
            - 27612.0 * a + 4540.0
        )
        + 2661120.0 * B[7] * (
            693680.0 * A[7] - 6473532.0 * A[6] + 13736424.0 * A[5]
            - 7047469.0 * A[4] - 723840.0 * A[3] + 471588.0 * a + 7376.0
        )
        - 2661120.0 * B[6] * (
            432536.0 * A[8] - 7850804.0 * A[7] + 27531114.0 * A[6]
            - 24234457.0 * A[5] - 703001.0 * A[4] + 3633474.0 * A[3]
            - 36244.0 * a - 45128.0
        )
        + 166320.0 * B[5] * (
            548912.0 * A[9] - 75660832.0 * A[8] + 502902712.0 * A[7]
            - 764807992.0 * A[6] + 91248287.0 * A[5] + 217811464.0 * A[4]
            - 20365384.0 * A[3] - 9776416.0 * a + 37936.0
        )
        + 10080.0 * B[4] * (
            18759728.0 * A[10] + 165932208.0 * A[9] - 4.71041844e9 * A[8]
            + 1.3686052536e10 * A[7] - 5.456818809e9 * A[6]
            - 6.834514245e9 * A[5] + 1.919299512e9 * A[4]
            + 7.52176152e8 * A[3] - 4.56612e7 * A[2] - 8.616848e6
        )
        - 360.0 * B[3] * (
            3.274336e7 * A[11] - 3.381871792e9 * A[10] - 2.1488827776e10 * A[9]
            + 2.00389923864e11 * A[8] - 1.9870800534e11 * A[7]
            - 1.71633799779e11 * A[6] + 1.23124874028e11 * A[5]
            + 4.0072774872e10 * A[4] - 9.13799328e9 * A[3]
            - 1.895843248e9 * a + 1.8929728e7
        )
        - 360.0 * b * Ap1[2] * (
            5.7685408e7 * A[11] + 4.06929456e8 * A[10] - 6.12537576e9 * A[9]
            - 2.709491892e10 * A[8] + 1.28752249410e11 * A[7]
            - 7.4866710561e10 * A[6] - 4.2917416470e10 * A[5]
            + 1.6256951352e10 * A[4] + 4.3752684e9 * A[3]
            - 3.16500688e8 * a - 4.7197152e7
        )
        + (a + 2.0) * (2.0 * a + 1.0) * (
            1.67898208e8 * A[11] - 2.2774946512e10 * A[10] - 8.8280004528e10 * A[9]
            + 6.11863976472e11 * A[8] + 1.041430242126e12 * A[7]
            - 3.446851131657e12 * A[6] + 1.041430242126e12 * A[5]
            + 6.11863976472e11 * A[4] - 8.8280004528e10 * A[3]
            - 2.2774946512e10 * a + 1.67898208e8
        )
    )

    C8 = C1 / (1.155790798848e14 * Ap1[8]) * (
          1.79159040e8 * B[15] - 1.254113280e9 * B[14] * (5.0 * a - 3.0)
        + 1.358622720e9 * B[13] * (70.0 * A[3] - 95.0 * a + 22.0)
        - 9.05748480e8 * B[12] * (904.0 * A[4] - 2109.0 * A[3] + 1119.0 * a - 112.0)
        + 1.245404160e9 * B[11] * (
            3532.0 * A[5] - 12824.0 * A[4] + 11829.0 * A[3] - 2824.0 * a + 44.0
        )
        - 5.9304960e7 * B[10] * (
            256820.0 * A[6] - 1397680.0 * A[5] + 2025545.0 * A[4]
            - 869495.0 * A[3] + 52000.0 * a + 8788.0
        )
        + 1.4826240e7 * B[9] * (
            2274536.0 * A[7] - 18601572.0 * A[6] + 40698318.0 * A[5]
            - 28230079.0 * A[4] + 3916398.0 * A[3] + 832668.0 * a - 65176.0
        )
        - 5.9304960e7 * B[8] * (
            760224.0 * A[8] - 9849164.0 * A[7] + 32495784.0 * A[6]
            - 34813869.0 * A[5] + 9175207.0 * A[4] + 1898688.0 * A[3]
            - 469788.0 * a - 13184.0
        )
        + 2.5945920e7 * B[7] * (
            1.167504e6 * A[9] - 2.8779840e7 * A[8] + 1.49752856e8 * A[7]
            - 2.46026112e8 * A[6] + 1.11944073e8 * A[5] + 1.83416e7 * A[4]
            - 1.2131496e7 * A[3] - 2.74368e5 * a + 1.028e5
        )
        - 1.57248e5 * B[6] * (
            1.2341872e7 * A[10] - 3.122991216e9 * A[9]
            + 2.9900054232e10 * A[8] - 7.8024816720e10 * A[7]
            + 5.8914656739e10 * A[6] + 4.637150811e9 * A[5]
            - 1.152340248e10 * A[4] + 2.36218968e8 * A[3]
            + 3.37923216e8 * a + 1.592048e6
        )
        - 2.8080e4 * B[5] * (
            2.65154912e8 * A[11] + 2.276098704e9 * A[10]
            - 1.05569461008e11 * A[9] + 4.96560666360e11 * A[8]
            - 6.27891462858e11 * A[7] + 4.1935358025e10 * A[6]
            + 2.03913875814e11 * A[5] - 2.3984801544e10 * A[4]
            - 1.3869306e10 * A[3] + 3.72786832e8 * a + 1.03532640e8
        )
        + 1.4400e3 * B[4] * (
            3.10292864e8 * A[12] - 5.5169117872e10 * A[11]
            - 3.58957020112e11 * A[10] + 5.714152556088e12 * A[9]
            - 1.3241597459352e13 * A[8] + 4.220720097141e12 * A[7]
            + 6.845418090249e12 * A[6] - 2.129559215808e12 * A[5]
            - 9.09225098472e11 * A[4] + 1.07518582576e11 * A[3]
            + 2.5619444368e10 * a - 1.13832704e8
        )
        + 12.0 * B[3] * (
            1.35319651136e11 * A[13] + 1.119107842176e12 * A[12]
            - 2.219351817432e13 * A[11] - 1.33421793595520e14 * A[10]
            + 8.60103051087996e14 * A[9] - 7.03353374803080e14 * A[8]
            - 7.04240127687381e14 * A[7] + 5.13111704637960e14 * A[6]
            + 1.66909061348316e14 * A[5] - 5.7671564069120e13 * A[4]
            - 1.2453426246e13 * A[3] + 6.95901207936e11 * a + 9.3786157376e10
        )
        - 12.0 * b * Ap1[2] * (
            4.365353408e9 * A[13] - 7.20248637504e11 * A[12]
            - 4.222331152560e12 * A[11] + 2.9413934270560e13 * A[10]
            + 1.32123980710980e14 * A[9] - 5.11247376962820e14 * A[8]
            + 2.83403639131779e14 * A[7] + 1.70415792320940e14 * A[6]
            - 7.9274388426588e13 * A[5] - 2.1009953050400e13 * A[4]
            + 3.284035340880e12 * A[3] + 5.89294339776e11 * a
            - 3.693760576e9
        )
        - (a + 2.0) * (2.0 * a + 1.0) * (
            3.4221025984e10 * A[13] + 2.26022948160e11 * A[12]
            - 5.067505612464e12 * A[11] - 1.8868361443936e13 * A[10]
            + 8.6215425028308e13 * A[9] + 1.43500920544692e14 * A[8]
            - 4.37682618704613e14 * A[7] + 1.43500920544692e14 * A[6]
            + 8.6215425028308e13 * A[5] - 1.8868361443936e13 * A[4]
            - 5.067505612464e12 * A[3] + 2.26022948160e11 * a
            + 3.4221025984e10
        )
    )

    C9 = C1 / (22191183337881600.0 * Ap1[9]) * (
          2149908480.0 * B[17] - 5733089280.0 * B[16]*(17*a - 11)
        + 7166361600.0 * B[15] * (272.0 * A[3] - 391*a + 104)
        - 3344302080.0 * B[14] * (6766.0 * A[4] - 16371.0 * A[3] + 9741.0 * a - 1306.0)
        + 1811496960.0 * B[13] * (
                             93092.0 * A[5] - 341564.0 * A[4] + 344199.0 * A[3] - 104924.0 * a
                             + 6308.0
                            )
        - 517570560.0 * B[12] * (
                           1626220.0 * A[6] - 8641508.0 * A[5] + 13274773.0 * A[4]
                           - 6952303.0 * A[3] + 1007420.0 * a + 5564.0
                          )
        + 284663808.0 * B[11] * (
                           9979136.0 * A[7] - 75766892.0 * A[6] + 169256148.0 * A[5]
                           - 136824959.0 * A[4] + 35714348.0 * A[3] - 463692.0 * a
                           - 293664.0
                          )
        - 1423319040.0 * B[10] * (
                            4466648.0 * A[8] - 49231116.0 * A[7] + 157507414.0 * A[6]
                            - 187114257.0 * A[5] + 78372295.0 * A[4]
                            - 4470082.0 * A[3] - 1913996.0 *a + 82424.0
                           )
        + 266872320.0 * B[9] * (
                          33133136.0 * A[9] - 564264544.0 * A[8]
                          + 2618606424.0 * A[7] - 4491310104.0 * A[6]
                          + 2853943765.0 * A[5] - 374694552.0 * A[4]
                          - 135365288.0 * A[3] + 17623968.0 * a + 696912.0
                         )
        - 2156544.0 * B[8] * (
                        2914256144.0 * A[10] - 93491712432.0 * A[9]
                        + 664876176984.0 * A[8] - 1661362937880.0 * A[7]
                        + 1563719627313.0 * A[6] - 382840842843.0 * A[5]
                        - 115399415640.0 * A[4] + 34565562936.0 * A[3]
                        + 1609337232.0 * a - 217321904.0
                       )
        + 179712.0 * B[7] * (
                       1266018560.0 * A[11] - 789261834512.0 * A[10]
                       + 10186841596896.0 * A[9] - 38877799073352.0 * A[8]
                       + 54334425968952.0 * A[7] - 22529574889533.0 * A[6]
                       - 5132942328000.0 * A[5] + 3438377465592.0 * A[4]
                       + 84287641248.0 * A[3] - 72493479440.0 * a - 807415936.0
                      )
        + 13824.0 * B[6] * (
                      156356794976.0 * A[12] + 1180898077328.0 * A[11]
                      - 90615270907936.0 * A[10] + 609258947056248.0 * A[9]
                      - 1312655191366722.0 * A[8] + 885900509321745.0 * A[7]
                      + 112162151855265.0 * A[6] - 212803071513258.0 * A[5]
                      + 6805217831352.0 * A[4] + 10051742651296.0 * A[3]
                      - 55035924848.0 * a - 52946379296.0
                     )
        - 576.0 * B[5] * (
                    143943926464.0 * A[13] - 60115486481856.0 * A[12]
                    - 376366989757200.0 * A[11] + 9534223075576160.0 * A[10]
                    - 35603777465262396.0 * A[9] + 39375990156664980.0 * A[8]
                    - 868175004137259.0 * A[7] - 14279180718355020.0 * A[6]
                    + 1985747535239364.0 * A[5] + 1264001337603680.0 * A[4]
                    - 75972792514320.0 * A[3] - 23855850572736.0 * a
                    - 4996648256.0
                   )
        - 384.0 * B[4] * (
                    2038525473856.0 * A[14] + 16057322146112.0 * A[13]
                    - 502133360559024.0 * A[12] - 2985686417468080.0 * A[11]
                    + 32418922182093292.0 * A[10] - 63665380623022452.0 * A[9]
                    + 16481208821092575.0 * A[8] + 34161547357596099.0 * A[7]
                    - 11490298497454932.0 * A[6] - 5117272758337156.0 * A[5]
                    + 933703210750480.0 * A[4] + 234855186762000.0 * A[3]
                    - 7860524600000.0 * a - 1226607567040.0
                   )
        + 96.0 * B[3] * (
                   324439754752.0 * A[15] - 77231415197120.0 * A[14]
                   - 539102931841856.0 * A[13] + 4618258299956336.0 * A[12]
                   + 28588485529469792.0 * A[11] - 141383982651179428.0 * A[10]
                   + 98783147840417772.0 * A[9] + 112831723492305801.0 * A[8]
                   - 83329761150975036.0 * A[7] - 26553582937192900.0 * A[6]
                   + 12469117738765952.0 * A[5] + 2587165396642160.0 * A[4]
                   - 340406368038080.0 * A[3] - 53659641606080.0 * a
                   + 219671272960.0
                  )
        + 96.0 * b * Ap1[2] * (
                       1026630779520.0 * A[15] + 8781958472768.0 * A[14]
                       - 210659786204384.0 * A[13] - 1222283505284208.0 * A[12]
                       + 5064251967491416.0 * A[11] + 24013052207628140.0 * A[10]
                       - 79710880160087370.0 * A[9] + 42596558293213227.0 * A[8]
                       + 26570293386695790.0 * A[7] - 14407831324576884.0 * A[6]
                       - 3617322833922440.0 * A[5] + 950664948554384.0 * A[4]
                       + 172358006894496.0 * A[3] - 7430887938496.0 * a
                       - 889746675584.0
                      )
        - (a + 2.0) * (2.0 * a + 1) * (
                             573840801152.0 * A[15]
                             - 156998277198784.0 * A[14]
                             - 898376974770592.0 * A[13]
                             + 8622589006459984.0 * A[12]
                             + 32874204024803560.0 * A[11]
                             - 111492707520083828.0 * A[10]
                             - 184768503480287646.0 * A[9]
                             + 528612016938984183.0 * A[8]
                             - 184768503480287646.0 * A[7]
                             - 111492707520083828.0 * A[6]
                             + 32874204024803560.0 * A[5]
                             + 8622589006459984.0 * A[4]
                             - 898376974770592.0 * A[3]
                             - 156998277198784.0 * a
                             + 573840801152.0
                            )
            )

    Z = (a * x)^(1.0 / Ap1[2])
    Zp = 1.0
    res = C1
    C = (C1, C2, C3, C4, C5, C6, C7, C8, C9)
    for k in 2:length(C)
        Zp /= Z
        res += (-1)^(k + 1) * C[k] * Zp
    end
    res *= Z^(0.5 - b) * exp(Ap1[2] / a * Z)
    return res
end


# Compute integrand Kmod(eps, a, b, x, r) for Gauss-Laguerre quadrature.
# K(a, b, x, r+eps) = exp(-r-eps) * Kmod(eps, a, b, x, r)
@inline function _Kmod(eps::Float64, a::Float64, b::Float64, x::Float64, r::Float64)
    x_r_a = x * (r + eps)^(-a)
    return exp(x_r_a * cospi(a)) * (r + eps)^(-b) * sin(x_r_a * sinpi(a) + π * b)
end


# Compute integrand P for Gauss-Legendre quadrature.
# P(eps, a, b, x, ϕ) =
#                      * exp(eps * cos(ϕ) + x * eps^(-a) * cos(a*ϕ))
#                      * cos(eps * sin(ϕ) - x * eps^(-a) * sin(a*ϕ)
#                            + (1-b)*ϕ)
@inline function _P(eps::Float64, a::Float64, b::Float64, x::Float64, ϕ::Float64)
    x_eps_a = x * eps^(-a)
    return exp(eps * cos(ϕ) + x_eps_a * cos(a * ϕ)) *
           cos(eps * sin(ϕ) - x_eps_a * sin(a * ϕ) + (1.0 - b) * ϕ)
end


# 5. Integral representation
#
# K(a, b, x, r) = exp(-r + x * r^(-a) * cos(π*a)) * r^(-b)
#               * sin(x * r^(-a) * sin(π*a) + π * b)
# P(eps, a, b, x, ϕ) =
#                      * exp(eps * cos(ϕ) + x * eps^(-a) * cos(a*ϕ))
#                      * cos(eps * sin(ϕ) - x * eps^(-a) * sin(a*ϕ)
#                            + (1-b)*ϕ)
#
# Φ(a, b, x) = 1/π * int_eps^inf K(a, b, x, r) * dr
#              + eps^(1-b)/π * int_0^π P(eps, a, b, x, ϕ) * dϕ
#
# for any eps > 0.
#
# Note that P has a misprint in Luchko (2008). Eq. 9, the cos(ϕ(beta-1)) at
# the end of the first line should be removed and the −sin(ϕ(beta−1)) at
# the end of the second line should read +(1-b)*ϕ.
# This integral representation introduced the free parameter eps (from the
# radius of complex contour integration). We try to choose eps such that
# the integrand behaves smoothly.
#
# As K has a leading exp(-r), we factor this out and apply Gauss-Laguerre
# quadrature rule:
#
# int_0^inf K(a, b, x, r+eps) dr = exp(-eps) int_0^inf exp(-r) Kmod(.., r) dr
#
# Note the shift r -> r+eps to have integration from 0 to infinity.
# The integral over P is done via a Gauss-Legendre quadrature rule.
#
# Note: Hardest argument range is large z, large b and small eps.

# roots of laguerre polynomial of order 50
# FastGaussQuadrature.gausslaguerre(50)[1]
const x_laguerre = [
    0.02863051833937908, 0.1508829356769337, 0.3709487815348964,
    0.6890906998810479, 1.105625023539913, 1.620961751102501,
    2.23561037591518, 2.950183366641835, 3.765399774405782,
    4.682089387559285, 5.70119757478489, 6.823790909794551,
    8.051063669390792, 9.384345308258407, 10.82510903154915,
    12.37498160875746, 14.03575459982991, 15.80939719784467,
    17.69807093335025, 19.70414653546156, 21.83022330657825,
    24.0791514444115, 26.45405784125298, 28.95837601193738,
    31.59588095662286, 34.37072996309045, 37.28751061055049,
    40.35129757358607, 43.56772026999502, 46.94304399160304,
    50.48426796312992, 54.19924488016862, 58.09682801724853,
    62.18705417568891, 66.48137387844482, 70.99294482661949,
    75.73701154772731, 80.73140480247769, 85.99721113646323,
    91.55969041253388, 97.44956561485056, 103.7048912366923,
    110.3738588076403, 117.5191982031112, 125.2254701334734,
    133.6120279227287, 142.8583254892541, 153.2603719726036,
    165.3856433166825, 180.6983437092145
]
# weights for laguerre polynomial of order 50
# FastGaussQuadrature.gausslaguerre(50)[2]
const w_laguerre = [
    0.07140472613518988, 0.1471486069645884, 0.1856716275748313,
    0.1843853825273539, 0.1542011686063556, 0.1116853699022688,
    0.07105288549019586, 0.04002027691150833, 0.02005062308007171,
    0.008960851203646281, 0.00357811241531566, 0.00127761715678905,
    0.0004080302449837189, 0.0001165288322309724, 2.974170493694165e-5,
    6.777842526542028e-6, 1.37747950317136e-6, 2.492886181720092e-7,
    4.010354350427827e-8, 5.723331748141425e-9, 7.229434249182665e-10,
    8.061710142281779e-11, 7.913393099943723e-12, 6.81573661767678e-13,
    5.13242671658949e-14, 3.365624762437814e-15, 1.913476326965035e-16,
    9.385589781827253e-18, 3.950069964503411e-19, 1.417749517827512e-20,
    4.309970276292175e-22, 1.101257519845548e-23, 2.344617755608987e-25,
    4.11854415463823e-27, 5.902246763596448e-29, 6.812008916553065e-31,
    6.237449498812102e-33, 4.452440579683377e-35, 2.426862352250487e-37,
    9.852971481049686e-40, 2.891078872318428e-42, 5.906162708112361e-45,
    8.01287459750397e-48, 6.789575424396417e-51, 3.308173010849252e-54,
    8.250964876440456e-58, 8.848728128298018e-62, 3.064894889844417e-66,
    1.988708229330752e-71, 6.049567152238783e-78
    ]
# roots of legendre polynomial of order 50
# FastGaussQuadrature.gausslegendre(50)[1]
const x_legendre = [
    -0.998866404420071, -0.9940319694320907, -0.9853540840480059,
    -0.9728643851066921, -0.9566109552428079, -0.9366566189448779,
    -0.9130785566557919, -0.885967979523613, -0.8554297694299461,
    -0.8215820708593359, -0.7845558329003993, -0.7444943022260685,
    -0.7015524687068223, -0.6558964656854394, -0.6077029271849502,
    -0.5571583045146501, -0.5044581449074642, -0.4498063349740388,
    -0.3934143118975651, -0.3355002454194374, -0.276288193779532,
    -0.2160072368760418, -0.1548905899981459, -0.09317470156008614,
    -0.03109833832718888, 0.03109833832718888, 0.09317470156008614,
    0.1548905899981459, 0.2160072368760418, 0.276288193779532,
    0.3355002454194374, 0.3934143118975651, 0.4498063349740388,
    0.5044581449074642, 0.5571583045146501, 0.6077029271849502,
    0.6558964656854394, 0.7015524687068223, 0.7444943022260685,
    0.7845558329003993, 0.8215820708593359, 0.8554297694299461,
    0.885967979523613, 0.9130785566557919, 0.9366566189448779,
    0.9566109552428079, 0.9728643851066921, 0.9853540840480059,
    0.9940319694320907, 0.998866404420071
    ]
# weights for legendre polynomial of order 50
# FastGaussQuadrature.gausslegendre(50)[2]
const w_legendre = [
    0.002908622553155141, 0.006759799195745401, 0.01059054838365097,
    0.01438082276148557, 0.01811556071348939, 0.02178024317012479,
    0.02536067357001239, 0.0288429935805352, 0.03221372822357802,
    0.03545983561514615, 0.03856875661258768, 0.0415284630901477,
    0.04432750433880328, 0.04695505130394843, 0.04940093844946632,
    0.05165570306958114, 0.05371062188899625, 0.05555774480621252,
    0.05718992564772838, 0.05860084981322245, 0.05978505870426546,
    0.06073797084177022, 0.06145589959031666, 0.06193606742068324,
    0.06217661665534726, 0.06217661665534726, 0.06193606742068324,
    0.06145589959031666, 0.06073797084177022, 0.05978505870426546,
    0.05860084981322245, 0.05718992564772838, 0.05555774480621252,
    0.05371062188899625, 0.05165570306958114, 0.04940093844946632,
    0.04695505130394843, 0.04432750433880328, 0.0415284630901477,
    0.03856875661258768, 0.03545983561514615, 0.03221372822357802,
    0.0288429935805352, 0.02536067357001239, 0.02178024317012479,
    0.01811556071348939, 0.01438082276148557, 0.01059054838365097,
    0.006759799195745401, 0.002908622553155141
    ]

function _wb_integral(a::Float64, b::Float64, x::Float64)
    # Fitted parameters for optimal choice of eps
    A = (0.41037, 0.30833, 6.9952, 18.382, -2.8566, 2.1122)

    # We use the free choice of eps to make the integral better behaved.
    # 1. Concern is oscillatory behaviour of P. Therefore, we'd like to
    #    make the change in the argument of cosine small, i.e. make arc length
    #    int_0^ϕ sqrt(1 + f'(ϕ)^2) dϕ small, with
    #    f(ϕ) = eps * sin(ϕ) - x * eps^(-a) * sin(a*ϕ) + (1-b)*ϕ
    #    Proxy, make |f'(ϕ)| small.
    # 2. Concern is int_0 K ~ int_0 (r+eps)^(-b) .. dr
    #    This is difficult as r -> 0  for large b. It behaves better for larger
    #    values of eps.

    # Minimize oscillatory behaviour of P
    eps = (A[1] * b * exp(-0.5 * a) +
           exp(A[2] + 1.0 / (1.0 + a) * log(x) - A[3] * exp(-A[4] * a) +
               A[5] / (1.0 + exp(A[6] * a))))
    if a >= 4.0 && x >= 100.0
        eps += 1.0  # This part is hard to fit.
    end

    # Large b
    if b >= 8.0
        # Make P small compared to K by setting eps large enough.
        # int K ~ exp(-eps) and int P ~ eps^(1-b)
        eps = max(eps, b^(-b / (1.0 - b)) + 0.1 * b)
    end

    # safeguard, higher better for larger a, lower better for tiny a.
    eps = min(eps, 150.0)
    eps = max(eps, 3.0)  # Note: 3 seems to be a pretty good choice in general.

    res1 = 0.0
    res2 = 0.0
    for k in eachindex(x_laguerre, w_laguerre, x_legendre)
        res1 += w_laguerre[k] * _Kmod(eps, a, b, x, x_laguerre[k])
        # y = (b-a)*(x+1)/2 + a  for integration from a=0 to b=π
        y = π * (x_legendre[k] + 1.0) / 2.0
        res2 += w_legendre[k] * _P(eps, a, b, x, y)
    end
    res1 *= exp(-eps)
    res2 *= π / 2.0
    res2 *= eps^(1.0 - b)

    return (res1 + res2) / π
end


"""
    wrightbessel(a::Float64, b::Float64, x::Float64)

Compute Wright's generalized Bessel function, defined as:

``\\Phi(a, b; x) = \\sum_{k=0}^\\infty \\frac{x^k}{k! \\Gamma(a k + b)}``

So far, only non-negative values of `a`, `b` and `x` (also noted as `ρ`, `β` and `z`)
are implemented.

External links:
[DLMF 10.46.E1](https://dlmf.nist.gov/10.46.E1), 
[Wikipedia](https://en.wikipedia.org/wiki/Bessel%E2%80%93Maitland_function)

# Implementation

One of 5 different computation methods is used depending on the ranges of the arguments:

1. Taylor series expansion in `x=0` ([DLMF 10.46.E1](https://dlmf.nist.gov/10.46.E1)), for `x <= 1`.
   Involves gamma functions in each term.
2. Taylor series expansion in `x=0` [(Dunn and Smyth 2005)](@cite dunn_smyth_2005), for large `a`.
3. Taylor series in `a=0`, for tiny `a` and not too large `x`.
4. Asymptotic expansion for large `x` ([Wright 1935](@cite wright_1935), [Paris 2017](@cite paris_2017)).
   Suitable for large `x` while still small `a` and `b`.
5. Integral representation [(Luchko 2008)](@cite luchko_2008), in principle for all arguments.

Accuracy is generally higher than 1e-11, though for some parameter values it can
be as low as 1e-8.
"""
function wrightbessel(a::Float64, b::Float64, x::Float64)
    if isnan(a) || isnan(b) || isnan(x)
        throw(ArgumentError("arguments cannot be NaN (got a=$a, b=$b, x=$x)"))
    elseif isinf(a) || isinf(b)
        throw(ArgumentError("`a` and `b` must be finite (got a=$a, b=$b)"))
    elseif a < 0.0 || b < 0.0 || x < 0.0
        throw(ArgumentError("arguments must be positive (got a=$a, b=$b, x=$x)"))
    elseif isinf(x)
        return Inf
    elseif a >= rgamma_zero || b >= rgamma_zero
        return NaN
    elseif x == 0.0
        return rgamma(b)
    elseif a == 0.0
        # return exp(x) * rgamma(b)
        return _exp_rgamma(x, b)
    elseif (a <= 1e-3 && b <= 50.0 && x <= 9.0) ||
           (a <= 1e-4 && b <= 70.0 && x <= 100.0) ||
           (a <= 1e-5 && b <= 170.0 && x < exp_inf)
        # Taylor Series expansion in a=0 to order=order => precision <= 1e-11
        # If beta is also small => precision <= 1e-11.
        # max order = 5
        if a <= 1e-5
            if x <= 1.0
                order = 2
            elseif x <= 10.0
                order = 3
            elseif x <= 100.0
                order = 4
            else # x < exp_inf:
                order = 5
            end
        elseif a <= 1e-4
            if x <= 1e-2
                order = 2
            elseif x <= 1.0
                order = 3
            elseif x <= 10.0
                order = 4
            else # x <= 100
                order = 5
            end
        else # a <= 1e-3
            if x <= 1e-5
                order = 2
            elseif x <= 1e-1
                order = 3
            elseif x <= 1.0
                order = 4
            else # x <= 9
                order = 5
            end
        end
        return _wb_small_a(a, b, x, order)
    elseif x <= 1.0
        # 18 term Taylor Series => error mostly smaller 5e-14
        return _wb_series(a, b, x, 0, 18)
    elseif x <= 2.0
        # 20 term Taylor Series => error mostly smaller 1e-12 to 1e-13
        return _wb_series(a, b, x, 0, 20)
    elseif a >= 5.0
        # Taylor series around the approximate maximum term.
        # Set number of terms=order.
        if a >= 10.0
            if x <= 1e11
                order = 6
            else
                order = min(floor(Int, log10(x) - 5 + b / 10), 30)
            end
        else
            if x <= 1e4
                order = 6
            elseif x <= 1e8
                order = floor(Int, 2 * log10(x))
            elseif x <= 1e10
                order = floor(Int, 4 * log10(x) - 16)
            else
                order = min(floor(Int, 6 * log10(x) - 36), 100)
            end
        end
        return _wb_large_a(a, b, x, order)
    elseif (0.5 <= a) && (a <= 1.8) && (100.0 <= b) && (1e5 <= x)
        return NaN
    elseif (a * x)^(1 / (1 + a)) >= 14 + b * b / (2 * (1 + a))
        # Asymptotic expansion in Z = (a*x)^(1/(1+a)) up to 8th term 1/Z^8.
        # For 1/Z^k, the highest term in b is b^(2*k) * a0 / (2^k k! (1+a)^k).
        # As a0 is a common factor to all orders, this explains a bit the
        # domain of good convergence set above.
        # => precision ~ 1e-11 but can go down to ~1e-8 or 1e-7
        # Note: We ensured a <= 5 as this is a bad approximation for large a.
        return _wb_asymptotic(a, b, x)
    else
        return _wb_integral(a, b, x)
    end
end