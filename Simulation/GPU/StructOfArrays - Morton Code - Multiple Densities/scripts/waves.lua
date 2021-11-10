

waves = function()

    local WAVE_NUM_DIR = {}
    local WAVE_PERIOD = {}
    local WAVE_STRENGTH = {}
    local TIMESTEP = {}

    local WAVE_CLEAR = {}
    getAttr("RENDERER","CURRENT","WAVE_CLEAR",0,WAVE_CLEAR)
    
    if WAVE_CLEAR[1] == 1.0 then
        local value = {}
        getAttr("RENDERER","CURRENT","GRAVITY",0,value)
        value[1] = 0 
        value[3] = 0
        setAttr("RENDERER","CURRENT","GRAVITY",0,value)

        getAttr("RENDERER","CURRENT","WAVE_NUM_DIR",0,WAVE_NUM_DIR)
        WAVE_NUM_DIR[1] = 0
        setAttr("RENDERER","CURRENT","WAVE_NUM_DIR",0,WAVE_NUM_DIR)
        WAVE_CLEAR[1] = 0.0
        setAttr("RENDERER","CURRENT","WAVE_CLEAR",0,WAVE_CLEAR)
    end

    getAttr("RENDERER","CURRENT","WAVE_NUM_DIR",0,WAVE_NUM_DIR)
    getAttr("RENDERER","CURRENT","WAVE_PERIOD",0,WAVE_PERIOD)
    getAttr("RENDERER","CURRENT","WAVE_STRENGTH",0,WAVE_STRENGTH)
    getAttr("RENDERER","CURRENT","TIMESTEP",0,TIMESTEP)

    TIMESTEP[1] = TIMESTEP[1]/0.01
    
    if WAVE_NUM_DIR[1] > 0.0 then
        local timer = {}
        local value = {}
        getAttr("RENDERER","CURRENT","GRAVITY",0,value)
        getAttr("RENDERER","CURRENT","Sim_Iters",0,timer)

        timer[1] = timer[1]*TIMESTEP[1]

        value[3] = WAVE_STRENGTH[1]*math.sin(timer[1]/math.max(WAVE_PERIOD[1],1))
        if WAVE_NUM_DIR[1] > 1.0 then
            value[1] = WAVE_STRENGTH[1]*math.cos(timer[1]/math.max(WAVE_PERIOD[1],1))
        else 
            value[1] = 0
        end
        setAttr("RENDERER","CURRENT","GRAVITY",0,value)
    end
end

 