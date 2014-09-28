function ZSHCH(ZR::Float64,ZI::Float64,CSHR::Float64,CSHI::Float64,CCHR::Float64,CCHI::Float64)
    CH::Float64 = 0
    CN::Float64 = 0
    SH::Float64 = 0
    SN::Float64 = 0
    SH = DSINH(ZR)
    CH = DCOSH(ZR)
    SN = DSIN(ZI)
    CN = DCOS(ZI)
    CSHR = SH * CN
    CSHI = CH * SN
    CCHR = CH * CN
    CCHI = SH * SN
    return (CSHR,CSHI,CCHR,CCHI)
end
