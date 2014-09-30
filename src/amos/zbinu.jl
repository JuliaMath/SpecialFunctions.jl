const _ZBINU_CWR = Array(Float64,2)
const _ZBINU_CWI = Array(Float64,2)
function ZBINU(ZR::Float64,ZI::Float64,FNU::Float64,KODE::Integer,N::Integer,CYR::AbstractArray{Float64},CYI::AbstractArray{Float64},NZ::Integer,RL::Float64,FNUL::Float64,TOL::Float64,ELIM::Float64,ALIM::Float64)
    AZ::Float64 = 0
    const CWI = _ZBINU_CWI
    const CWR = _ZBINU_CWR
    DFNU::Float64 = 0
    I::Int32 = 0
    INW::Int32 = 0
    NLAST::Int32 = 0
    NN::Int32 = 0
    NUI::Int32 = 0
    NW::Int32 = 0
    ZEROI::Float64 = 0
    ZEROR::Float64 = 0
    begin 
        ZEROR = 0.0
        ZEROI = 0.0
    end
    NZ = 0
    AZ = ZABS(COMPLEX(ZR,ZI))
    NN = N
    DFNU = FNU + DBLE(FLOAT(N - 1))
    if AZ <= 2.0
        @goto line10
    end
    if AZ * AZ * 0.25 > DFNU + 1.0
        @goto line20
    end
    @label line10
    (NW,) = ZSERI(ZR,ZI,FNU,KODE,NN,CYR,CYI,NW,TOL,ELIM,ALIM)
    INW = IABS(NW)
    NZ = NZ + INW
    NN = NN - INW
    if NN == 0
        return NZ
    end
    if NW >= 0
        @goto line120
    end
    DFNU = FNU + DBLE(FLOAT(NN - 1))
    @label line20
    if AZ < RL
        @goto line40
    end
    if DFNU <= 1.0
        @goto line30
    end
    if AZ + AZ < DFNU * DFNU
        @goto line50
    end
    @label line30
    (NW,) = ZASYI(ZR,ZI,FNU,KODE,NN,CYR,CYI,NW,RL,TOL,ELIM,ALIM)
    if NW < 0
        @goto line130
    end
    @goto line120
    @label line40
    if DFNU <= 1.0
        @goto line70
    end
    @label line50
    (NW,) = ZUOIK(ZR,ZI,FNU,KODE,1,NN,CYR,CYI,NW,TOL,ELIM,ALIM)
    if NW < 0
        @goto line130
    end
    NZ = NZ + NW
    NN = NN - NW
    if NN == 0
        return NZ
    end
    DFNU = FNU + DBLE(FLOAT(NN - 1))
    if DFNU > FNUL
        @goto line110
    end
    if AZ > FNUL
        @goto line110
    end
    @label line60
    if AZ > RL
        @goto line80
    end
    @label line70
    (NW,) = ZMLRI(ZR,ZI,FNU,KODE,NN,CYR,CYI,NW,TOL)
    if NW < 0
        @goto line130
    end
    @goto line120
    @label line80
    (NW,) = ZUOIK(ZR,ZI,FNU,KODE,2,2,CWR,CWI,NW,TOL,ELIM,ALIM)
    if NW >= 0
        @goto line100
    end
    NZ = NN
    for I = 1:NN
        CYR[I] = ZEROR
        CYI[I] = ZEROI
        @label line90
    end
    return NZ
    @label line100
    if NW > 0
        @goto line130
    end
    (NW,) = ZWRSK(ZR,ZI,FNU,KODE,NN,CYR,CYI,NW,CWR,CWI,TOL,ELIM,ALIM)
    if NW < 0
        @goto line130
    end
    @goto line120
    @label line110
    NUI = INT(SNGL(FNUL - DFNU)) + 1
    NUI = MAX0(NUI,0)
    (NW,NLAST) = ZBUNI(ZR,ZI,FNU,KODE,NN,CYR,CYI,NW,NUI,NLAST,FNUL,TOL,ELIM,ALIM)
    if NW < 0
        @goto line130
    end
    NZ = NZ + NW
    if NLAST == 0
        @goto line120
    end
    NN = NLAST
    @goto line60
    @label line120
    return NZ
    @label line130
    NZ = -1
    if NW == -2
        NZ = -2
    end
    return NZ
end
