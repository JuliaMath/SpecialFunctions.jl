function ZUCHK(YR::Float64,YI::Float64,NZ::Int32,ASCLE::Float64,TOL::Float64)
    SS::Float64 = zero(Float64)
    ST::Float64 = zero(Float64)
    WI::Float64 = zero(Float64)
    WR::Float64 = zero(Float64)
    NZ = Int32(0)
    WR = abs(YR)
    WI = abs(YI)
    ST = min(WR,WI)
    if ST > ASCLE
        return NZ
    end
    SS = max(WR,WI)
    ST = ST / TOL
    if SS < ST
        NZ = Int32(1)
    end
    return NZ
end
