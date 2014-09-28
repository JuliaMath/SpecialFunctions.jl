function ZUCHK(YR::Float64,YI::Float64,NZ::Integer,ASCLE::Float64,TOL::Float64)
    SS::Float64 = 0
    ST::Float64 = 0
    WI::Float64 = 0
    WR::Float64 = 0
    NZ = 0
    WR = DABS(YR)
    WI = DABS(YI)
    ST = DMIN1(WR,WI)
    if ST > ASCLE
        return NZ
    end
    SS = DMAX1(WR,WI)
    ST = ST / TOL
    if SS < ST
        NZ = 1
    end
    return NZ
end
