function ZLOG(AR::Float64,AI::Float64,BR::Float64,BI::Float64,IERR::Int32)
    DHPI::Float64 = zero(Float64)
    DPI::Float64 = zero(Float64)
    DTHETA::Float64 = zero(Float64)
    ZM::Float64 = zero(Float64)
    begin 
        DPI = 3.141592653589793
        DHPI = 1.5707963267948966
    end
    IERR = Int32(0)
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
    if AI == 0.0
        @goto line60
    end
    BI = DHPI
    BR = DLOG(DABS(AI))
    if AI < 0.0
        BI = -BI
    end
    return (BR,BI,IERR)
    @label line20
    if AR > 0.0
        @goto line30
    end
    BR = DLOG(DABS(AR))
    BI = DPI
    return (BR,BI,IERR)
    @label line30
    BR = DLOG(AR)
    BI = 0.0
    return (BR,BI,IERR)
    @label line40
    if AR < 0.0
        DTHETA = DTHETA + DPI
    end
    @label line50
    ZM = abs(COMPLEX(AR,AI))
    BR = DLOG(ZM)
    BI = DTHETA
    return (BR,BI,IERR)
    @label line60
    IERR = Int32(1)
    return (BR,BI,IERR)
end
