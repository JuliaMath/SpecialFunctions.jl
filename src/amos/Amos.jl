module Amos

 int = Int
 int32 = Int32
 float32 = Float32
 float64 = Float64
 complex64 = ComplexF32
 complex128 = ComplexF64
itrunc(x) = trunc(Int,x)

 D1MACH1 = 2.23e-308
 D1MACH2 = 1.79e-308
 D1MACH3 = 1.11e-16
 D1MACH4 = 2.22e-16
 D1MACH5 = 0.3010299956639812

 I1MACH1 = int32(5)
 I1MACH2 = int32(6)
 I1MACH3 = int32(0)
 I1MACH4 = int32(0)
 I1MACH5 = int32(32)
 I1MACH6 = int32(4)
 I1MACH7 = int32(2)
 I1MACH8 = int32(31)
 I1MACH9 = int32(2147483647)
 I1MACH10 = int32(2)
 I1MACH11 = int32(24)
 I1MACH12 = int32(-125)
 I1MACH13 = int32(127)
 I1MACH14 = int32(53)
 I1MACH15 = int32(-1021)
 I1MACH16 = int32(1023)

# fortran intrinsics
 DABS = abs
 DATAN = atan
 DBLE = float64
 DCOS = cos
 DCOSH = cosh
 DEXP = exp
 DLOG = log
 DMAX1 = max
 DMIN1 = min
 DSIGN = copysign
 DSIN = sin
 DSINH = sinh
 DSQRT = sqrt
 COMPLEX = complex128
 FLOAT = float32
 IABS = abs
 INT = itrunc
 MAX0 = max
 MIN0 = min
 MOD = mod
 SNGL = float32

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
