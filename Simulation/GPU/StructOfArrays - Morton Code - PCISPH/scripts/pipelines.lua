

SimController = function()
    local pause = {}
    local simIters = {}
    local simMaxIters = {}
    getAttr("RENDERER","CURRENT","Pause",0,pause)
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