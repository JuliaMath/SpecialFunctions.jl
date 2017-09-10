module Amos

const D1MACH1 = 2.23e-308
const D1MACH2 = 1.79e-308
const D1MACH3 = 1.11e-16
const D1MACH4 = 2.22e-16
const D1MACH5 = 0.3010299956639812

const I1MACH1 = Int32(5)
const I1MACH2 = Int32(6)
const I1MACH3 = Int32(0)
const I1MACH4 = Int32(0)
const I1MACH5 = Int32(32)
const I1MACH6 = Int32(4)
const I1MACH7 = Int32(2)
const I1MACH8 = Int32(31)
const I1MACH9 = Int32(2147483647)
const I1MACH10 = Int32(2)
const I1MACH11 = Int32(24)
const I1MACH12 = Int32(-125)
const I1MACH13 = Int32(127)
const I1MACH14 = Int32(53)
const I1MACH15 = Int32(-1021)
const I1MACH16 = Int32(1023)

const SNGL = Float32
const DBLE = Float64
const FLOAT = Float32
const COMPLEX = Complex128

DSIGN(x::Float64, y::Float64)::Float64 = copysign(x, y)
DSQRT(x::Float64)::Float64 = sqrt(x)
DABS(x::Float64)::Float64 = abs(x)
DMAX1(x::Float64, y::Float64)::Float64 = max(x, y)
DMIN1(x::Float64, y::Float64)::Float64 = min(x, y)
DSIN(x::Float64)::Float64 = sin(x)
DCOS(x::Float64)::Float64 = cos(x)
DSINH(x::Float64)::Float64 = sinh(x)
DCOSH(x::Float64)::Float64 = cosh(x)
DATAN(x::Float64)::Float64 = atan(x)
DEXP(x::Float64)::Float64 = exp(x)
DLOG(x::Float64)::Float64 = log(x)

INT(x::Float32)::Int32 = trunc(Int32, x)
IABS(x::Int32)::Int32 = abs(x)
MOD(x::Int32,  y::Int32)::Int32 = mod(x, y)
MAX0(x::Int32, y::Int32)::Int32 = max(x, y)
MIN0(x::Int32, y::Int32)::Int32 = min(x, y)

include("dgamln.jl")
include("zabs.jl")
include("zacai.jl")
include("zacon.jl")
include("zairy.jl")
include("zasyi.jl")
include("zbesh.jl")
include("zbesi.jl")
include("zbesj.jl")
include("zbesk.jl")
include("zbesy.jl")
include("zbinu.jl")
include("zbiry.jl")
include("zbknu.jl")
include("zbuni.jl")
include("zbunk.jl")
include("zdiv.jl")
include("zexp.jl")
include("zkscl.jl")
include("zlog.jl")
include("zmlri.jl")
include("zmlt.jl")
include("zrati.jl")
include("zs1s2.jl")
include("zseri.jl")
include("zshch.jl")
include("zsqrt.jl")
include("zuchk.jl")
include("zunhj.jl")
include("zuni1.jl")
include("zuni2.jl")
include("zunik.jl")
include("zunk1.jl")
include("zunk2.jl")
include("zuoik.jl")
include("zwrsk.jl")

end # module
