function ZBUNK(ZR::Float64,ZI::Float64,FNU::Float64,KODE::Int32,MR::Int32,N::Int32,YR::AbstractArray{Float64},YI::AbstractArray{Float64},NZ::Int32,TOL::Float64,ELIM::Float64,ALIM::Float64)
    AX::Float64 = zero(Float64)
    AY::Float64 = zero(Float64)
    NZ = Int32(0)
    AX = abs(ZR) * 1.7321
    AY = abs(ZI)
    if AY > AX
        @goto line10
    end
    (NZ,) = ZUNK1(ZR,ZI,FNU,KODE,MR,N,YR,YI,NZ,TOL,ELIM,ALIM)
    @goto line20
    @label line10
    (NZ,) = ZUNK2(ZR,ZI,FNU,KODE,MR,N,YR,YI,NZ,TOL,ELIM,ALIM)
    @label line20
    return NZ
end
