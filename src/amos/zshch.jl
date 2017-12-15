function ZSHCH(ZR::Float64,ZI::Float64,CSHR::Float64,CSHI::Float64,CCHR::Float64,CCHI::Float64)
    CH::Float64 = zero(Float64)
    CN::Float64 = zero(Float64)
    SH::Float64 = zero(Float64)
    SN::Float64 = zero(Float64)
    SH = sinh(ZR)
    CH = cosh(ZR)
    SN = sin(ZI)
    CN = cos(ZI)
    CSHR = SH * CN
    CSHI = CH * SN
    CCHR = CH * CN
    CCHI = SH * SN
    return (CSHR,CSHI,CCHR,CCHI)
end
