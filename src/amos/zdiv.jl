function ZDIV(AR::Float64,AI::Float64,BR::Float64,BI::Float64,CR::Float64,CI::Float64)
    BM::Float64 = 0
    CA::Float64 = 0
    CB::Float64 = 0
    CC::Float64 = 0
    CD::Float64 = 0
    BM = 1.0 / ZABS(COMPLEX(BR,BI))
    CC = BR * BM
    CD = BI * BM
    CA = (AR * CC + AI * CD) * BM
    CB = (AI * CC - AR * CD) * BM
    CR = CA
    CI = CB
    return (CR,CI)
end
