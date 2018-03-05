# Create test set for elliptic functions

import mpmath

mpmath.mp.prec = 256

def isplusinf(x): return x == mpmath.mp.inf or (x > 0 and x.exp + x.bc >= 2**63)
def isminusinf(x): return isplusinf(-x)
def isinf(x): return isplusinf(x) or isminusinf(x)

def parse_int8(x):
    return x.to_bytes(1, "big", signed=True)

def parse_int64(x):
    return x.to_bytes(8, "big", signed=True)

def reverse_bits(x,nbits):
    y = 0
    for i in range(nbits):
        y <<= 1
        y += (x&1)
        x >>= 1
    return y

def parse_mantissa(x):
    if x == 0:
        return bytes(32)
    m = reverse_bits(x.man, x.man.bit_length())
    b = m.to_bytes(32,"big",signed=False)
    return bytes(map(lambda b: reverse_bits(b,8), b))

def parse_exponent(x):
    return parse_int64(x.exp+x.bc)

def parse_mpf(x):
    s = bytes()
    if x < 0 and not isminusinf(x):
        s = parse_int8(-1)
    if x == 0:
        s = parse_int8(0)
    if x > 0 and not isplusinf(x):
        s = parse_int8(1)
    if isminusinf(x):
        s = parse_int8(2)
        x = mpmath.mpf(0)
    if isplusinf(x):
        s = parse_int8(3)
        x = mpmath.mpf(0)
    if mpmath.isnan(x):
        s = parse_int8(4)
        x = mpmath.mpf(0)
    return s + parse_mantissa(x) + parse_exponent(x)

def parse_mpc(x):
    return parse_mpf(x.real) + parse_mpf(x.imag)


def test_io():
    with open("io.bin","wb") as f:
        f.write(parse_int64(0))
        f.write(parse_int64(1))
        f.write(parse_int64(-1))
        f.write(parse_mpf(mpmath.mpf( 0)))
        f.write(parse_mpf(mpmath.mpf( 1)))
        f.write(parse_mpf(mpmath.mpf( 2)))
        f.write(parse_mpf(mpmath.mpf(-1)))
        f.write(parse_mpc(mpmath.mpc(1j)))
        f.write(parse_mpf(mpmath.sqrt(mpmath.mpf(2))))
        f.write(parse_mpf(mpmath.mp.inf))
        f.write(parse_mpf(-mpmath.mp.inf))
        f.write(parse_mpf(mpmath.mp.nan))
test_io()


def dump_testset(name,u,m):
    with open(name,"wb") as f:
        f.write(parse_int64(len(u)*len(m)))
        for ui in u:
            for mi in m:
                sn,cn,dn = ellipj(ui,mi)
                f.write(parse_mpc(ui))
                f.write(parse_mpc(mi))
                f.write(parse_mpc(sn))
                f.write(parse_mpc(cn))
                f.write(parse_mpc(dn))

def ellipj(u,m):
    sn = mpmath.ellipfun("sn",u,m=m)
    cn = mpmath.ellipfun("cn",u,m=m)
    dn = mpmath.ellipfun("dn",u,m=m)
    if u == 0: sn = 0*sn # https://github.com/fredrik-johansson/mpmath/issues/386
    return sn,cn,dn

def test_ellipj():
    s = [
        mpmath.mpc(1,0),
        mpmath.sqrt(mpmath.mpc(0,1)),
        mpmath.mpc(0,1),
        mpmath.mpc(-1,0),
    ]
    r = [
        mpmath.mpf(mpmath.mp.eps),
        mpmath.sqrt(mpmath.mp.eps),
        1/mpmath.mp.e,
        mpmath.mpf("1"),
        mpmath.mp.e,
    ]
    x = [mpmath.mpc(0)] + [si*ri for si in s for ri in r]
    dump_testset("ellipj.bin",x, x+[1+mi for mi in x])
test_ellipj()


def dump_testset(name,m):
    with open(name,"wb") as f:
        f.write(parse_int64(len(m)))
        for mi in m:
            K = mpmath.ellipk(mi)
            f.write(parse_mpc(mi))
            f.write(parse_mpc(K))

def test_K():
    ns = 16
    s = [
        mpmath.exp(2*mpmath.mp.pi*1j*i/ns)
        for i in range(ns)
    ]
    r = [
        mpmath.mpf(mpmath.mp.eps),
        mpmath.sqrt(mpmath.mp.eps),
        1/mpmath.mp.pi,
        1/mpmath.mp.e,
        mpmath.mpf("0.5"),
        mpmath.mpf("1"),
        mpmath.mpf("2"),
        mpmath.mp.e,
        mpmath.mp.pi,
    ]
    m = [mpmath.mpc(0)] + [si*ri for si in s for ri in r]
    dump_testset("K.bin",m + [1+mi for mi in m])
test_K()
