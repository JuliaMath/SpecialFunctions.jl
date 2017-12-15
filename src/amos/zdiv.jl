function ZDIV(AR::Float64,AI::Float64,BR::Float64,BI::Float64,CR::Float64,CI::Float64)
    BM::Float64 = zero(Float64)
    CA::Float64 = zero(Float64)
    CB::Float64 = zero(Float64)
    CC::Float64 = zero(Float64)
    CD::Float64 = zero(Float64)
    BM = 1.0 / abs(COMPLEX(BR,BI))
    CC = BR * BM
    CD = BI * BM
    CA = (AR * CC + AI * CD) * BM
    CB = (AI * CC - AR * CD) * BM
    CR = CA
    CI = CB
    return (CR,CI)
end
