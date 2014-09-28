function ZABS(Z::Complex128)
    __ZABS__::Float64 = 0
    Q::Float64 = 0
    S::Float64 = 0
    U::Float64 = 0
    V::Float64 = 0
    ZI::Float64 = 0
    ZR::Float64 = 0
    ZR = real(Z)
    ZI = imag(Z)
    U = DABS(ZR)
    V = DABS(ZI)
    S = U + V
    S = S * 1.0
    if S == 0.0
        @goto line20
    end
    if U > V
        @goto line10
    end
    Q = U / V
    __ZABS__ = V * DSQRT(1.0 + Q * Q)
    return __ZABS__
    @label line10
    Q = V / U
    __ZABS__ = U * DSQRT(1.0 + Q * Q)
    return __ZABS__
    @label line20
    __ZABS__ = 0.0
    return __ZABS__
end
