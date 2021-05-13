

LoopController = function()
    local pause = {}
    local Pcisph_Iters = {}
    local Pcisph_MaxIters = {}
    local simIters = {}
    local simMaxIters = {}

    getAttr("RENDERER","CURRENT","Pause",0,pause)
    getAttr("RENDERER","CURRENT","Pcisph_Iters",   0,Pcisph_Iters)
    getAttr("RENDERER","CURRENT","Pcisph_MaxIters",0,Pcisph_MaxIters)
    getAttr("RENDERER","CURRENT","Sim_Iters",   0,simIters)
    getAttr("RENDERER","CURRENT","Sim_MaxIters",0,simMaxIters)

    if simIters[1] >= simMaxIters[1] and simMaxIters[1] > 0 then
        pause[1] = 1
        setAttr("RENDERER","CURRENT","Pause",0,pause)
    else
        if pause[1] == 0 then
            simIters[1] = simIters[1]+1
            
            setAttr("RENDERER","CURRENT","Sim_Iters",   0,simIters) 
        end
    end

    if Pcisph_Iters[1] >= Pcisph_MaxIters[1] and Pcisph_MaxIters[1] > 0 then
        -- poe iters a 0 e devolve falso para sair da pipeline. Da proxima vez que entrar na pipeline ele vai correr outras 7 vezes
        Pcisph_Iters[1] = 0
        setAttr("RENDERER","CURRENT","Pcisph_Iters",0,Pcisph_Iters)

        return false;
    else
        if pause[1] == 0 then
            Pcisph_Iters[1] = Pcisph_Iters[1]+1
            setAttr("RENDERER","CURRENT","Pcisph_Iters",   0,Pcisph_Iters) 
        end
    end
    return true
end

Print = function()
    print("Print no lua")
    return true
end


Iter = 0

SmoothPre = function ()
    Iter = 0    
end

SmoothPipeline = function()
    local numIter = {}
    getAttr("RENDERER","CURRENT","Blur_PingPong",0,numIter)

    if Iter < numIter[1] then
        Iter = Iter + 1
        return true
    end
    return false
end