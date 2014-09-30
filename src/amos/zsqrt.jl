function ZSQRT(AR::Float64,AI::Float64,BR::Float64,BI::Float64)
    DPI::Float64 = zero(Float64)
    DRT::Float64 = zero(Float64)
    DTHETA::Float64 = zero(Float64)
    ZM::Float64 = zero(Float64)
    begin 
        DRT = 0.7071067811865476
        DPI = 3.141592653589793
    end
    ZM = ZABS(COMPLEX(AR,AI))
    ZM = DSQRT(ZM)
    if AR == 0.0
        @goto line10
    end
    if AI == 0.0
        @goto line20
    end
    DTHETA = DATAN(AI / AR)
    if DTHETA <= 0.0
        @goto line40
    end
    if AR < 0.0
        DTHETA = DTHETA - DPI
    end
    @goto line50
    @label line10
    if AI > 0.0
        @goto line60
    end
    if AI < 0.0
        @goto line70
    end
    BR = 0.0
    BI = 0.0
    return (BR,BI)
    @label line20
    if AR > 0.0
        @goto line30
    end
    BR = 0.0
    BI = DSQRT(DABS(AR))
    return (BR,BI)
    @label line30
    BR = DSQRT(AR)
    BI = 0.0
    return (BR,BI)
    @label line40
    if AR < 0.0
        DTHETA = DTHETA + DPI
    end
    @label line50
    DTHETA = DTHETA * 0.5
    BR = ZM * DCOS(DTHETA)
    BI = ZM * DSIN(DTHETA)
    return (BR,BI)
    @label line60
    BR = ZM * DRT
    BI = ZM * DRT
    return (BR,BI)
    @label line70
    BR = ZM * DRT
    BI = -ZM * DRT
    return (BR,BI)
end
