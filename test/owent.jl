# Owen's T function tests

using Test
using IrrationalConstants
using LinearAlgebra
using SpecialFunctions

# test values for accurate and precise calculation
hvec = [0.0625, 6.5, 7.0, 4.78125, 2.0, 1.0, 0.0625, 1, 1, 1, 1, 0.5, 0.5, 0.5, 0.5, 0.25, 0.25, 0.25, 0.25, 0.125, 0.125, 0.125, 0.125, 0.0078125
, 0.0078125, 0.0078125, 0.0078125, 0.0078125, 0.0078125, 0.0625, 0.5, 0.9, 2.5, 7.33, 0.6, 1.6, 2.33, 2.33]
avec = [0.25, 0.4375, 0.96875, 0.0625, 0.5, 0.9999975, 0.999999125, 0.5, 1, 2, 3, 0.5, 1, 2, 3, 0.5, 1, 2, 3, 0.5, 1, 2, 3, 0.5, 1, 2, 3, 10, 100
, 0.999999999999999, 0.999999999999999, 0.999999999999999, 0.999999999999999, 0.999999999999999, 0.999999999999999, 0.999999999999999, 0.999999999999999
, 0.99999]
cvec = [big"0.0389119302347013668966224771378", big"2.00057730485083154100907167685e-11", big"6.399062719389853083219914429e-13"
, big"1.06329748046874638058307112826e-7", big"0.00862507798552150713113488319155", big"0.0667418089782285927715589822405"
, big"0.1246894855262192"
, big"0.04306469112078537", big"0.06674188216570097", big"0.0784681869930841", big"0.0792995047488726", big"0.06448860284750375", big"0.1066710629614485"
, big"0.1415806036539784", big"0.1510840430760184", big"0.07134663382271778", big"0.1201285306350883", big"0.1666128410939293", big"0.1847501847929859"
, big"0.07317273327500386", big"0.1237630544953746", big"0.1737438887583106", big"0.1951190307092811", big"0.07378938035365545"
, big"0.1249951430754052", big"0.1761984774738108", big"0.1987772386442824", big"0.2340886964802671", big"0.2479460829231492"
, big"0.1246895548850743676554299881345328280176736760739893903915691894"
, big"0.1066710629614484543187382775527753945264849005582264731161129477"
, big"0.0750909978020473015080760056431386286348318447478899039422181015"
, big"0.0030855526911589942124216949767707430484322201889568086810922629"
, big"5.7538182971139187466647478665179637676531179007295252996453e-14", big"0.0995191725772188724714794470740785702586033387786949658229016920"
, big"0.0258981646643923680014142514442989928165349517076730515952020227"
, big"0.0049025023268168675126146823752680242063832053551244071400100690"
, big"0.0049024988349089450612896251009169062698683918433614542387524648"]

# test values for type stability checking
ht = [0.0625, 0.0, 0.5, 0.5, 0.5]
at = [0.025, 0.5, 0.0, 1.0, +Inf]

#=
# Interactive tests, displayin error and error magnitude
mvbck = "\033[1A\033[58G"
mva = "\033[72G"
mve = "\033[90G"

warmup = owent(0.0625,0.025)
warmup = owent(0.0625,0.999999125)

for i in 1:size(hvec,1)
    if i == 1
        println("\t\tExecution Time","\033[58G","h",mva,"a",mve,"error\t\t","log10 error")
    end
    h = hvec[i]
    a = avec[i]
    c = cvec[i]
    print(i,"\t")
    @time tx = owent(h,a)
    err = Float64(tx-c)
    logerr = Float64(round(log10(abs(tx-c)),sigdigits=3))
    println(mvbck,h,mva,a,mve,err,"\t",logerr)
end
=#

@testset "Owen's T" begin

    # check that error for calc is within specification
    @testset "Owen's T value checks" begin

        for i in 1:size(hvec,1)
            h = hvec[i] # h test value
            a = avec[i] # a test value
            c = cvec[i] # the "correct" answer
            t = owent(h,a)

            err = round(log10(abs(t-c)),digits=0)

            @test err ≤ -16.0
        end

    end

    @testset "Owen's T Shortcut Evaluations" begin
        @test owent(-0.0625,0.025) == owent(0.0625,0.025) # if h<0 equals owent(abs(h),a)
        @test owent(0.0,0.025) == atan(0.025)*inv2π  # if h=0, t = atan(a)*inv2π
        @test owent(0.0625,0.0) == 0.0  # when a=0, owent=0
        @test owent(0.0625,-0.025) == -owent(0.0625,0.025) # if a<0, t = -owent(h,abs(a))
        @test owent(0.0625,1.0) == 0.125*erfc(-0.0625*invsqrt2)*erfc(0.0625*invsqrt2) # if a=1, t = (formula)
        @test owent(0.0625,+Inf) == 0.25*erfc(sqrt(0.0625^2)*invsqrt2) # if a=∞, t = (formula)
        @test owent(0.0625,-Inf) == -0.25*erfc(sqrt(0.0625^2)*invsqrt2) # should also work, due to a<0 condition
    end

    @testset "Owen's T type stability" begin

        for T1 in [Float16, Float32, Float64, BigFloat]
            for T2 in [BigFloat, Float64, Float32, Float16]
                for i in 1:size(ht,1)
                    h=T1(ht[i])
                    a=T2(at[i])
                    (p1, p2) = promote(h,a)
                    t=owent(h,a)
                    @test typeof(t) == typeof(p1)
                    #println("T1: ",T1," ",typeof(p1),"\tT2: ",T2," ",typeof(p2),"\t",i,"\th ",h,"\ta ",a,"\tt ",t,"\t",typeof(t)," ",typeof(t)==typeof(p1))
                end
            end
        end

    end # test type stability

end
# Owen's T Function Tests
