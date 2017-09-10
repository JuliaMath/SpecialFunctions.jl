function ZMLT(AR::Float64,AI::Float64,BR::Float64,BI::Float64,CR::Float64,CI::Float64)
    CA::Float64 = zero(Float64)
    CB::Float64 = zero(Float64)
    CA = AR * BR - AI * BI
    CB = AR * BI + AI * BR
    CR = CA
    CI = CB
    return (CR,CI)
end
