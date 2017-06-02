using SpecialFunctions
using Base.Test

const SF = SpecialFunctions

# test functions, see runtests.jl
relerr(z, x) = z == x ? 0.0 : abs(z - x) / abs(x)
relerrc(z, x) = max(relerr(real(z),real(x)), relerr(imag(z),imag(x)))
≅(a,b) = relerrc(a,b) ≤ 1e-9

# test functions, similar to runtests.jl, but can't do relative error when some values are zero
abserr(z, x) = z == x ? 0.0 : abs(z - x)
abserrc(z, x) = max(abserr(real(z),real(x)), abserr(imag(z),imag(x)))
≒(a,b) = abserrc(a,b) ≤ 1e-10

S = Array(Complex, (0,))
Z = Array(Complex, (0,))
errors = Array(Complex, (0,))
rel_errors = Array(Complex, (0,))
L = 0
L2 = 0
diff = 0
rel_diff = 0

file = "polylog_testdata_alpha.dat"
f = open(file, "r")
for line in eachline(f)
    if ismatch(r"#", line)
    else
        v = split(line, ',')

        s = parse(Float64, v[1]) + im*parse(Float64, v[2])
        z = parse(Float64, v[3]) + im*parse(Float64, v[4])
        L = parse(Float64, v[5]) + im*parse(Float64, v[6])
        L2 = SF.polylog(s,z,1.0e-18)
        diff = L - L2
        rel_diff = complex(real(diff)/real(L),  imag(diff)/imag(L))
        
        push!(S, s)
        push!(Z, z)
        push!(errors, diff)
        push!(rel_errors, rel_diff)

        # println("  L = $L,  L2 =$L2,  diff = $diff,  \t\t rel diff = $rel_diff")
        println("  s=$s, z=$z, |error|=$(ceil(log10(abs(diff))))")
    end
end

# test particular values from Wolfram alpha
# @testset "polylog tests from Wolfram alpha" begin

#     @test SF.polylog(-5.000001, 0.2) ≅ 6.904305364987799724230921116283822405784472717903797757356
#     @test SF.polylog(-2.000001, 0.2) ≅ 0.468750244180354668389455477288861311399345621141103816518
#     @test SF.polylog(-1.000001, 0.2) ≅ 0.312500094184052926057717828023167764362038234067425674984
#     @test SF.polylog( 0.000000, 0.2) ≅ 0.25
#     @test SF.polylog( 0.000001, 0.2) ≅ 0.249999960605813559100141040899469061220130457763090931203
#     @test SF.polylog( 1.000001, 0.2) ≅ 0.223143533840628659834227863703211050928449673084530824304
#     @test SF.polylog( 2.000001, 0.2) ≅ 0.211003767368667938482167666994614717744311348708200162427
#     @test SF.polylog( 3.000001, 0.2) ≅ 0.205324191902667822690683746355210013921766645655171501890
#     @test SF.polylog( 5.000001, 0.2) ≅ 0.201284594885756384926736867628954299128217148725340685679
    
#     @test SF.polylog(-5.000001, 0.99) ≅ 1.16439554606503726735945485177408748716681200000005e14
#     @test SF.polylog(-2.000001, 0.99) ≅ 1.97011088076187708022430311310006082286580942497213e6
#     @test SF.polylog(-1.000001, 0.99) ≅ 9900.04972775403451817385172007212640686495346816117
#     @test SF.polylog( 0.000000, 0.99) ≅ 99
#     @test SF.polylog( 0.000001, 0.99) ≅ 98.9995988050883517986218058007486858307489619231956
#     @test SF.polylog( 1.000001, 0.99) ≅ 4.605161353580737753746071585981515640203854412492030454580 # error ~ 1.0e11
#     @test SF.polylog( 2.000001, 0.99) ≅ 1.588624649826159899085393338960947692670814871114246941712 # error ~ 1.0e11
#     @test SF.polylog( 3.000001, 0.99) ≅ 1.185832744102043627496839259642313182039620569301974825714
#     @test SF.polylog( 5.000001, 0.99) ≅ 1.026110449210262375841139988837721699091019919425544651696

#     @test SF.polylog( 1.000001, 1.20) ≅ 1.60944129470676213105637136369425072009730867713452098010 - 3.14158912002728351076393334268714060723117930694360052182im
#     @test SF.polylog( 1.000001, 1.50) ≅ 0.69315092620532843831758654220946253504678190003195418383 - 3.1415916309839162783984857611820275451978840493033368692im
#     @test SF.polylog( 1.000001, 1.70) ≅ 0.35667861585275033749153552050383961466011713941159572502 - 3.1415924761565638528178063690254522808240730779853094090im
#     @test SF.polylog( 1.000001, 1.999999) ≒ 
#     @test SF.polylog( 1.000001, 2.00) ≅ 
#     @test SF.polylog( 1.000001, 2.10) ≅ -0.0953067628118213777754205736836756047381793027303437193 - 3.14159352922832327133509242559303053858610422117044987im

# end    
