function ZEXP(AR::Float64,AI::Float64,BR::Float64,BI::Float64)
    CA::Float64 = zero(Float64)
    CB::Float64 = zero(Float64)
    ZM::Float64 = zero(Float64)
    ZM = DEXP(AR)
    CA = ZM * DCOS(AI)
    CB = ZM * DSIN(AI)
    BR = CA
    BI = CB
    return (BR,BI)
end
