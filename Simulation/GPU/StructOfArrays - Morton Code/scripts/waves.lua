

waves = function()

    local WAVE_NUM_DIR = {}
    local WAVE_PERIOD = {}
    local WAVE_STRENGTH = {}
    local TIMESTEP = {}

    local WAVE_CLEAR = {}
    getAttr("RENDERER","CURRENT","WAVE_CLEAR",0,WAVE_CLEAR)
    
    if WAVE_CLEAR[1] == 1.0 then
        local SAVEXMIN={}
        local SAVEZMIN={}
        getAttr("RENDERER","CURRENT","SAVEXMIN",0,SAVEXMIN)
        getAttr("RENDERER","CURRENT","SAVEZMIN",0,SAVEZMIN)
        setAttr("RENDERER","CURRENT","XMIN",0,SAVEXMIN)
        setAttr("RENDERER","CURRENT","ZMIN",0,SAVEZMIN)

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
        local XMIN = {}
        local ZMIN = {}
        local SAVEXMIN={}
        local SAVEZMIN={}
        getAttr("RENDERER","CURRENT","SAVEXMIN",0,SAVEXMIN)
        getAttr("RENDERER","CURRENT","SAVEZMIN",0,SAVEZMIN)
    
        getAttr("RENDERER","CURRENT","TIMER",0,timer)

        timer[1] = timer[1]*TIMESTEP[1]

        value[1] = WAVE_STRENGTH[1]*math.sin(timer[1]/math.max(WAVE_PERIOD[1],1))
        if WAVE_NUM_DIR[1] > 1.0 then
            value[3] = WAVE_STRENGTH[1]*math.cos(timer[1]/math.max(WAVE_PERIOD[1],1))
            ZMIN[1] = SAVEZMIN[1] +  math.abs(value[3]);
            setAttr("RENDERER","CURRENT","ZMIN",0,ZMIN)
        else 
            value[3] = 0
        end
        
        XMIN[1] = SAVEXMIN[1] + math.abs(value[1]);
        setAttr("RENDERER","CURRENT","XMIN",0,XMIN)
    end
end

 