#version 460
#define GLSLIFY 1

#pragma glslify: import ../collisionSDF.glsl

const vec3 inf = vec3(1e20, 1e20, 1e20);

#define LIGHTINGFRONT 1
#define FRONTANDBACK  2
#define THICKNESS     3

uniform float       fov;
uniform float       Particle_Thickness;
uniform float       densityThreshold;
uniform float       eta;
uniform float       stepmultiplier;
uniform float       foamThreshold;
uniform int         DebugParticles;
uniform int         maxrefrat;
uniform int         SkyboxType;
uniform int         shlickpower;
uniform int         specular_power;
uniform mat4        m_pv;
uniform sampler2D   ground;
uniform sampler2D   height;
uniform sampler2D   foam;
uniform sampler3D   voxelRep;
uniform samplerCube skybox;
uniform vec2        windowSize;
uniform vec3        lightDir;
uniform vec3        right;
uniform vec3        up;
uniform vec3        view;
uniform vec3        simmax;
uniform vec3        simmin;
uniform vec3        voxdimvec;
uniform vec3        voxmax;
uniform vec3        voxmin;
uniform vec4        position;
uniform vec4        Particle_Color;
uniform vec4        Thin_Color;
uniform vec4        foamcolor;

in vec2 uv;
out vec4 color;

vec2 getVoxelVal(vec3 index){
    float tol = 0.01;
    if(index.x > simmax.x - tol || index.y > simmax.y - tol || index.z > simmax.z - tol 
    || index.x < simmin.x + tol || index.y < simmin.y + tol || index.z < simmin.z + tol)
        return vec2(0.0);
    // if(abs(index.x) > simmin.x-tol || abs(index.y) > simmin.y-tol || abs(index.z) > simmin.z-tol )
    //     return 0.0;
    //index = (index-voxmin)/(voxmax-voxmin) ;// clamp(, vec3(0.0), vec3(1.0));
    index = clamp((index-voxmin)/(voxmax-voxmin),vec3(0.0),vec3(1.0)) ;// clamp(, vec3(0.0), vec3(1.0));

    return texture(voxelRep,index   ).xy;
    // return voxels[index.x + voxDim*(index.y + voxDim*index.z)];
}

vec4 getFoamNormal(vec3 index, float off){
    vec3 offset = vec3(off);
    vec3 normal = vec3(0);
    normal.x = 1.0 * getVoxelVal(index-vec3(offset.x,0,0)).y- getVoxelVal(index+vec3(offset.x,0,0)).y;
    normal.y = 1.0 * getVoxelVal(index-vec3(0,offset.y,0)).y- getVoxelVal(index+vec3(0,offset.y,0)).y;    
    normal.z = 1.0 * getVoxelVal(index-vec3(0,0,offset.z)).y- getVoxelVal(index+vec3(0,0,offset.z)).y;
    
    
    if (normal.x == 0 && normal.y == 0 && normal.z == 0 ){
        return vec4(0);
    }
    return vec4(normalize(normal),0); 
}

vec4 getNormal(vec3 index, float off){
    vec3 offset =  vec3(off);
    vec3 normal = vec3(0);
    normal.x = 1.0 * getVoxelVal(index-vec3(offset.x,0,0)).x- getVoxelVal(index+vec3(offset.x,0,0)).x;
    normal.y = 1.0 * getVoxelVal(index-vec3(0,offset.y,0)).x- getVoxelVal(index+vec3(0,offset.y,0)).x;    
    normal.z = 1.0 * getVoxelVal(index-vec3(0,0,offset.z)).x- getVoxelVal(index+vec3(0,0,offset.z)).x;
    
    
    if (normal.x == 0 && normal.y == 0 && normal.z == 0 ){
        return vec4(0);
    }
    return vec4(normalize(normal),0); 
}

vec3 rayDir()
{
    float aspectRatio = windowSize.x / windowSize.y;

    vec2 screenCoords = uv * 2.0 - 1.0;
    float height = tan(radians(fov) * 0.5);

    //vec3 right  = normalize(cross(view, vec3(0.0, 1.0, 0.0)));
    //vec3 up     = normalize(cross(right, view));
    vec3 rayDir = view + (right * height * aspectRatio * screenCoords.x) + (up * height * screenCoords.y);
    //rayDir = vec3(0, 0, -1) + (vec3(1, 0, 0) * height * aspectRatio * screenCoords.x) + (vec3(0, 1, 0) * height * screenCoords.y);
    return normalize(rayDir);
}

/*
 * Tutorial from https://www.scratchapixel.com/lessons/3d-basic-rendering/minimal-ray-tracer-rendering-simple-shapes/ray-box-intersection
 */
vec3 rayAABBIntersection(vec3 rayPos, vec3 rayDir)
{
    float tmin, tmax, tymin, tymax, tzmin, tzmax;
    vec3 invDir = 1.0 / rayDir;
    vec3 bounds[2] = { simmin, simmax };
    int signs[3] = { invDir.x < 0 ? 1 : 0, invDir.y < 0 ? 1 : 0, invDir.z < 0 ? 1 : 0 };
 
    tmin  = (bounds[  signs[0]].x - rayPos.x) * invDir.x;
    tmax  = (bounds[1-signs[0]].x - rayPos.x) * invDir.x;
    tymin = (bounds[  signs[1]].y - rayPos.y) * invDir.y;
    tymax = (bounds[1-signs[1]].y - rayPos.y) * invDir.y;
 
    if ((tmin > tymax) || (tymin > tmax))
        return inf;
    if (tymin > tmin)
        tmin = tymin;
    if (tymax < tmax)
        tmax = tymax;
 
    tzmin = (bounds[  signs[2]].z - rayPos.z) * invDir.z;
    tzmax = (bounds[1-signs[2]].z - rayPos.z) * invDir.z;
 
    if ((tmin > tzmax) || (tzmin > tmax))
        return inf;
    if (tzmin > tmin)
        tmin = tzmin;
    if (tzmax < tmax)
        tmax = tzmax;

    // if tmin is negative, it intersected behind the ray
    return tmin > 0.0 ? rayPos + rayDir * tmin : inf;
}

bool isInsideBoundingBox(vec3 pos,float ctol)
{
    float tol = ctol;
    return pos.x > simmin.x -tol && pos.x < simmax.x + tol && 
           pos.z > simmin.z -tol && pos.z < simmax.z + tol && 
           pos.y > simmin.y -tol && pos.y < simmax.y + tol;
}

bool isInsideBoundingBox(vec3 pos)
{
    float tol = 0.001;
    return pos.x > simmin.x -tol && pos.x < simmax.x + tol && 
           pos.z > simmin.z -tol && pos.z < simmax.z + tol && 
           pos.y > simmin.y -tol && pos.y < simmax.y + tol;
}

uniform float HeightMultiplier;
uniform int   Render;
uniform vec4  AUXVEC1;
uniform vec4  AUXVEC2;
uniform vec4  AUXVEC3;

vec4 getSkybox(vec3 raypos,vec3 rayDir){

    //if(rayDir.y < 0 && raypos.y >= -2){
    //    vec3 p = raypos + rayDir*(abs((raypos.y+1.97))/-rayDir.y);
    //    return texture(ground,p.xz/4.0);
    //    //return color
    //}
    if( SkyboxType == 0)
        return texture(skybox, normalize(-1*rayDir));
    else 
    if( SkyboxType == 1)
        return vec4(1);
    return vec4(vec3(0),1);
}

float getheight(vec2 pos,vec2 rep){
    float heightmult = HeightMultiplier+0.05;
    //float y = (pnoise(pos,rep));
    //y =  y* smoothstep(2.1,1.6,abs(pos.x))*smoothstep(2.9,2.0,abs(pos.y));
    return heightmult*
        texture(height,clamp(((pos.xy-voxmin.xz)/(voxmax.xz-voxmin.xz)).xy,vec2(0.0),vec2(1.0))).x
        ;//*smoothstep(0,1 ,abs(y) );
}

    

bool colisionCheck( vec3 current,vec3 currentStep,inout vec3 color)
{
    
    float tol = 0.05;
    vec2 rep = vec2(10);

    vec3 next = (current + currentStep).xyz;
    vec2 pos = next.xz;
    float off = 0.01;

    float y = getheight(pos,rep);
    //    return false;
    
    if ( ( Render == 0 || Render ==2 ) && mapBoder(next)<-0.0000){
        color = (0.5+0.5*dot(-calcNormalBorder(current),lightDir.xyz))*vec3(1.0);
        return true;
    }

    if ( ( Render == 0 || Render ==2 ) && mapSDF(next)<-0.0000){
        color = (0.5+0.5*dot(-calcNormal(current),lightDir.xyz))*vec3(0.0,0.5,0.5);
        return true;
    }

     if ( (Render < 2) && ((next).y < -2.0 + HeightMultiplier) &&  ((next).y < -1.99+y )){
        vec3 t  = vec3(+off,(getheight(pos+vec2(off,0),rep)),0.0);
        vec3 b  = vec3(-off,(getheight(pos-vec2(off,0),rep)),0.0);
        vec3 l  = vec3(0,   (getheight(pos+vec2(0,off),rep)),off);
        vec3 r  = vec3(0,   (getheight(pos-vec2(0,off),rep)),-off);

        
	    vec3 normal = cross(normalize(t-b),normalize(l-r));
        float intensity = 0.3+0.7*dot(normal,lightDir);
        color = intensity*texture(ground,(clamp(((pos.xy-voxmin.xz)/(voxmax.xz-voxmin.xz)).xy,vec2(0.0),vec2(1.0)))).xyz;//vec3(0.6);
        // color = (current/2)/2+0.5;
        // color *= 2;
        return true;
    }
    return false;
}

vec3 worldpos2clip(vec3 current){
    current.z -=2.9;
    vec4 clip_space_pos = ((m_pv*vec4(current.xyz,1)));
    clip_space_pos = clip_space_pos/clip_space_pos.w;
    clip_space_pos.z  = (clip_space_pos.z +1.0)/2.0;
    return vec3(clip_space_pos.xyz);
}

vec4 rayMarchTranslucent(vec3 rayPos, vec3 rayDir, int numRefracs)
{

    vec4 cout = vec4(0.0);
    float weight = 1.0;

    vec3 accum = vec3(0.0);
    // float world_step = texel_step * (sceneLf - sceneRt);
    float marchingStep = 5*stepmultiplier* length(((voxmax-voxmin)/voxdimvec));
    bool inFluid = false;
    float lastH = 0;
    float lastZ = 0;
    vec3 firstnormal = vec3(0.0);
    float threshold = densityThreshold;

    vec3 current = rayPos;
    vec3 firstPos;
    vec3 lastPos = current;
    vec4 totspec = vec4(0);
                float etaInverse = 1.0/eta;
    

    if(getVoxelVal(current).x > threshold){
        inFluid = true;
        firstnormal = normalize( 1-rayDir);
        firstPos = current;

        
        
    }
    int maxrefracts = maxrefrat;
    if (DebugParticles == FRONTANDBACK )
        maxrefracts = 1;
    vec4 fog;
    bool colided = false;
    bool colidedFoam = false;
    for (int i = 1000; i < 2000; i++)
    {
        if (!isInsideBoundingBox(current)) 
            break;

        vec3 color; 
        
        //TODO calculate step of ray based on texel size and ray direction
        float st = (marchingStep + 2*exp(i/2000)*marchingStep);
        vec3 currentStep = rayDir * st;

        vec3 last = current;
        current += currentStep;
        float currentLength = length(current - position.xyz);

        if (colisionCheck(current,currentStep,color)) {
            if(inFluid)
                 accum += 0.05*Particle_Thickness*(1.0-Particle_Color.xyz);
                 
            cout.xyz += color.xyz -(accum);
            colided = true;
            break;

        
             //break;
        }
        vec2 voxelval = getVoxelVal(current);
        // CHECK IF VOXELw
        if ( voxelval.x > threshold){
            
            lastPos = current;
            accum += voxelval.x*Particle_Thickness*st*(1.0-Particle_Color.xyz);
            

            //if ( !colidedFog && voxelval.y > 0.25)
            //    fog +=   10.0*marchingStep*max(vec4(0.0),(voxelval.y) );a
            //    else
            
            if (voxelval.y > foamThreshold && !colidedFoam){
                vec4 normal = getFoamNormal(current,st);
                vec4 col;  
                vec2 aux = normalize(normal.xy);
                //col = mix(mix(
                //      texture(foam,((current.xz)/0.5)), 
                //      texture(foam,((current.yz)/0.5)),abs(normal.x) ), 
                //      texture(foam,((current.xy)/0.5)),abs(normal.z) );
                col = mix(mix(
                      texture(foam,(2.0*(current.xz+0.05*(normal.xz*current.y))))  , 
                      texture(foam,(2.0*(current.yz+0.05*(normal.yz*current.x))))  ,abs(normal.x) ), 
                      texture(foam,(2.0*(current.xy+0.05*(normal.xy*current.z))))  ,abs(normal.z) );

                //col = texture(foam,((((current.xz))+0.1*sin((current.y)*(normal.xz)))));
               
                if(col.a > 0.95){

                    //return col;
                    col = weight*(col)-vec4(accum,0.0)+cout;
                    colidedFoam = true;
                    return mix(
                               ((foamcolor*col) )
                        ,Thin_Color*((foamcolor*col) ),Thin_Color.a
                    );
                }
            }

        }
        

        if (numRefracs < maxrefracts  && (
                ( voxelval.x > threshold && !inFluid)
                || (voxelval.x <= threshold && inFluid)
            ) 
        ){
           
            //vec3 currentStep2 = rayDir * (i-1) * marchingStep;
            vec3 current2 = current - currentStep;
            vec4 normal = getNormal(current2,st);// + vec4(0.01,0.0,0.0,0.0);
            if(inFluid)
                normal = -1*normal;

            float shlick;
            if (inFluid){
                float f = pow(((1 - etaInverse) * (1 - etaInverse)) / ((1 + etaInverse) * (1 + etaInverse)),2);
		        shlick = f + (1 - f) * pow(1 - dot(-1*rayDir, normal.xyz), shlickpower);
         
            }
            else {
                float f = pow(((1 - eta) * (1 - eta)) / ((1 + eta) * (1 + eta)),2);
		        shlick = f + (1 - f) * pow(1 - dot(-1*rayDir, normal.xyz), shlickpower);
            }
            shlick = clamp(shlick,0.0 , 1.0);    
            //if(shlick < 0.7){

            vec3 h = normalize(lightDir - rayDir);
		    // compute the specular intensity
		    float intSpec = max(dot(h, -1*normal.xyz), 0.0);
		    // compute the specular term into spec
		    vec4 spec = vec4(1) * pow(intSpec, specular_power);
            vec3 tol = vec3(0.01,0.01,0.01);
            
            
            cout += weight*(shlick)*getSkybox(current,reflect(rayDir,normal.xyz)) -vec4(accum,0.0);
            
            totspec += weight*spec-vec4(accum,0.0);

            weight *= (1-shlick);
           
            numRefracs++;
            if(inFluid){
                rayDir = normalize(refract(rayDir, normal.xyz, etaInverse));
            }
            else
                rayDir = normalize(refract(rayDir, normal.xyz, eta));
            inFluid = !inFluid;
        
            
            if( DebugParticles == LIGHTINGFRONT)
                     return Particle_Color*vec4(0.2+0.8*clamp(dot(-1*normal.xyz,lightDir), 0.0,1.0));
                
            if ((firstnormal.x == 0 && firstnormal.y == 0 && firstnormal.z == 0 )){
                firstnormal = normal.xyz;
                firstPos = current2;
            }
            
            //return 0.5*(normal+1);
            // return refract(rayDir,normal ,0.7 );
            
            
            
            //}
        }

        //    //return(vec4((current),1));w
        if (rayDir == vec3(0.0))
            return getSkybox(current,rayDir);
        if(accum.x > 1.0 && accum.y > 1.0 && accum.z > 1.0){
            return cout+weight*(1-vec4(accum,0.0));
        }

    }

    
    if ( length(accum) == 0 && !colided)
        return getSkybox(current, rayDir);

    //return Particle_Color*vec4(0.2+0.8*clamp(dot(getNormal(lastPos).xyz,lightDir), 0.0,1.0));
    if (DebugParticles == THICKNESS)
        return vec4(1-vec4(accum,0.0) );

    if (DebugParticles == FRONTANDBACK){
        vec3 normal = -1*getNormal(lastPos,marchingStep).xyz;
        float f = pow(((1 - etaInverse) * (1 - etaInverse)) / ((1 + etaInverse) * (1 + etaInverse)),2);
		float shlick = f + (1 - f) * pow(1 - dot(-1*rayDir, normal.xyz), shlickpower);
        if(!colided)

        cout = cout + clamp(

        
            mix(
                (getSkybox(lastPos,refract(rayDir, normal.xyz, etaInverse))-vec4(accum,0.0)), 
                (getSkybox(lastPos,reflect(rayDir, normal.xyz))-vec4(accum,0.0)), 
                shlick
            ), 
            vec4(0.0),vec4(1.0)
        )
        
        ;
    }
    else
        if(!colided)
            cout = cout + clamp(    weight*(getSkybox(current,rayDir)-vec4(accum,0.0))
                                    ,vec4(0.0),vec4(1.0)
            );
    // if (abs(lastPos.x) < 1.9 && abs(lastPos.y) < 1.9 && abs(lastPos.z) < 1.9)
    //     return mix(vec4(0.8),cout,smoothstep(0.0, 0.1*Particle_Thickness, accum ));
    //if(!colided)
    if(accum != vec3(0))
        cout.xyz = mix(cout.xyz, Thin_Color.xyz * cout.xyz, smoothstep(0, 1, (Thin_Color.a)));
    
    
    //cout += fog;
    cout = mix(cout,fog,fog.a);
    if (fog.a<0.2){

        cout += totspec;
    }
    
    //cout.xyz += fog.xyz;
    //return vec4(accum * Particle_Color);
    return cout;
    
}

vec4 rayMarch(vec3 rayPos, vec3 rayDir)
{

    vec4 cout = vec4(0.0);
    float weight = 1.0;

    vec3 accum = vec3(0.0);
    // float world_step = texel_step * (sceneLf - sceneRt);
    float marchingStep = stepmultiplier* length(((voxmax-voxmin)/voxdimvec));
    bool inFluid = false;
    float lastH = 0;
    float lastZ = 0;
    vec3 firstnormal = vec3(0.0);
    float threshold = densityThreshold;

    vec3 current = rayPos;
    vec3 firstPos;
    vec3 lastPos = current;
    vec4 totspec = vec4(0);
                float etaInverse = 1.0/eta;
    int numRefracs = 0;

    if(getVoxelVal(current).x > threshold){
        inFluid = true;
        firstnormal = normalize( 1-rayDir);
        firstPos = current;

        if (DebugParticles == FRONTANDBACK )
            numRefracs = 1;
        
    }
    int maxrefracts = maxrefrat;
    if (DebugParticles == FRONTANDBACK )
        maxrefracts = 1;
    vec4 fog;
    bool colided = false;
    bool colidedFoam = false;
    for (int i = 0; i < 2000; i++)
    {
        if (!isInsideBoundingBox(current)) 
            break;

        vec3 color; 
        
        //TODO calculate step of ray based on texel size and ray direction
        float st = (marchingStep + 2*exp(i/2000)*marchingStep);
        vec3 currentStep = rayDir * st;
        //vec3 currentStep = rayDir * (marchingStep + i*0.01*marchingStep);
        vec3 last = current;
        current += currentStep;
        float currentLength = length(current - position.xyz);

        if (colisionCheck(current,currentStep,color)) {
            if(inFluid)
                 accum += 0.05*Particle_Thickness*(1.0-Particle_Color.xyz);
                 
            cout.xyz += color.xyz -(accum);
            colided = true;
            break;

        
             //break;
        }
        vec2 voxelval = getVoxelVal(current);
        // CHECK IF VOXELw
        if ( voxelval.x > threshold){
            
            lastPos = current;
            accum += voxelval.x*Particle_Thickness*st*(1.0-Particle_Color.xyz);
            

            //if ( !colidedFog && voxelval.y > 0.25)
            //    fog +=   10.0*marchingStep*max(vec4(0.0),(voxelval.y) );a
            //    else
            
            if (voxelval.y > foamThreshold && !colidedFoam){
                vec4 normal = getFoamNormal(current,st);
                vec4 col;  
                vec2 aux = normalize(normal.xy);
                //col = mix(mix(
                //      texture(foam,((current.xz)/0.5)), 
                //      texture(foam,((current.yz)/0.5)),abs(normal.x) ), 
                //      texture(foam,((current.xy)/0.5)),abs(normal.z) );
                col = mix(mix(
                      texture(foam,(2.0*(current.xz+0.05*(normal.xz*current.y))))  , 
                      texture(foam,(2.0*(current.yz+0.05*(normal.yz*current.x))))  ,abs(normal.x) ), 
                      texture(foam,(2.0*(current.xy+0.05*(normal.xy*current.z))))  ,abs(normal.z) )
                      ;

                //col = texture(foam,((((current.xz))+0.1*sin((current.y)*(normal.xz)))));
               
                if(col.a > 0.95){

                    //return col;
                    colidedFoam = true;
                    col = weight*(col)-vec4(accum,0.0)+cout;
                    return mix(
                               ((foamcolor*col) )
                        ,Thin_Color*((foamcolor*col) ),Thin_Color.a
                    );
                }
            }

        }
        

        if (numRefracs < maxrefracts  && (
                ( voxelval.x > threshold && !inFluid)
                || (voxelval.x <= threshold && inFluid)
            ) 
        ){
           
            //vec3 currentStep2 = rayDir * (i-1) * marchingStep;
            vec3 current2 = current - currentStep;
            vec4 normal = getNormal(current2,st);// + vec4(0.01,0.0,0.0,0.0);
            if(inFluid)
                normal = -1*normal;

            float shlick;
            if (inFluid){
                float f = pow(((1 - etaInverse) * (1 - etaInverse)) / ((1 + etaInverse) * (1 + etaInverse)),2);
		        shlick = f + (1 - f) * pow(1 - dot(-1*rayDir, normal.xyz), shlickpower);
         
            }
            else {
                float f = pow(((1 - eta) * (1 - eta)) / ((1 + eta) * (1 + eta)),2);
		        shlick = f + (1 - f) * pow(1 - dot(-1*rayDir, normal.xyz), shlickpower);
            }
            shlick = clamp(shlick,0.0 , 1.0);    
            //if(shlick < 0.7){

            vec3 h = normalize(lightDir - rayDir);
		    // compute the specular intensity
		    float intSpec = max(dot(h, -1*normal.xyz), 0.0);
		    // compute the specular term into spec
		    vec4 spec = vec4(1) * pow(intSpec, specular_power);
            vec3 tol = vec3(0.01,0.01,0.01);
            
            
            if( numRefracs < 1)
                cout += clamp (weight*(shlick)*rayMarchTranslucent(current2,reflect(rayDir,normal.xyz),numRefracs+2)  -vec4(accum,0.0),
                        vec4(0.0),
                        vec4(1.0)
                );
            else
                cout += weight*(shlick)*getSkybox(current,reflect(rayDir,normal.xyz)) -vec4(accum,0.0);
            
            totspec += weight*spec-vec4(accum,0.0);

            weight *= (1-shlick);
           
            numRefracs++;
            if(inFluid){
                rayDir = normalize(refract(rayDir, normal.xyz, etaInverse));
            }
            else
                rayDir = normalize(refract(rayDir, normal.xyz, eta));
            inFluid = !inFluid;
        
            
            if( DebugParticles == LIGHTINGFRONT)
                     return Particle_Color*vec4(0.2+0.8*clamp(dot(-1*normal.xyz,lightDir), 0.0,1.0));
                
            if ((firstnormal.x == 0 && firstnormal.y == 0 && firstnormal.z == 0 )){
                firstnormal = normal.xyz;
                firstPos = current2;
            }
            
            //return 0.5*(normal+1);
            // return refract(rayDir,normal ,0.7 );
            
            
            
            //}
        }

        //    //return(vec4((current),1));w
        if (rayDir == vec3(0.0))
            return getSkybox(current,rayDir);
        if(accum.x > 1.0 && accum.y > 1.0 && accum.z > 1.0){
            return cout+weight*(1-vec4(accum,0.0));
        }

    }

    
    if ( length(accum) == 0 && !colided)
        return getSkybox(current, rayDir);

    //return Particle_Color*vec4(0.2+0.8*clamp(dot(getNormal(lastPos).xyz,lightDir), 0.0,1.0));
    if (DebugParticles == THICKNESS)
        return vec4(1-vec4(accum,0.0) );

    if (DebugParticles == FRONTANDBACK){
        vec3 normal = -1*getNormal(lastPos,marchingStep).xyz;
        float f = pow(((1 - etaInverse) * (1 - etaInverse)) / ((1 + etaInverse) * (1 + etaInverse)),2);
		float shlick = f + (1 - f) * pow(1 - dot(-1*rayDir, normal.xyz), shlickpower);
        if(!colided)

        cout = cout + clamp(

        
            mix(
                (getSkybox(lastPos,refract(rayDir, normal.xyz, etaInverse))-vec4(accum,0.0)), 
                (getSkybox(lastPos,reflect(rayDir, normal.xyz))-vec4(accum,0.0)), 
                shlick
            ), 
            vec4(0.0),vec4(1.0)
        )
        
        ;
    }
    else
        if(!colided)
            cout = cout + clamp(    weight*(getSkybox(current,rayDir)-vec4(accum,0.0))
                                    ,vec4(0.0),vec4(1.0)
            );
    // if (abs(lastPos.x) < 1.9 && abs(lastPos.y) < 1.9 && abs(lastPos.z) < 1.9)
    //     return mix(vec4(0.8),cout,smoothstep(0.0, 0.1*Particle_Thickness, accum ));
    //if(!colided)
    if(accum != vec3(0))
        cout.xyz = mix(cout.xyz, Thin_Color.xyz * cout.xyz, smoothstep(0, 1, (Thin_Color.a)));
    
    
    //cout += fog;
    cout = mix(cout,fog,fog.a);
    if (fog.a<0.2){

        cout += totspec;
    }
    
    //cout.xyz += fog.xyz;
    //return vec4(accum * Particle_Color);
    return cout;
    
}

void main()
{
    const vec3 rayDir = normalize(rayDir());
    vec3 rayPos = position.xyz + rayDir*0.02;
    color = getSkybox(rayPos,rayDir);
    
    
    
    // Ray eye = Ray( rayPos, normalize(rayDir) );
    // AABB aabb = AABB(vec3(-2.0), vec3(+2.0));
    //if(abs(rayPos.x) > 5 || abs(rayPos.y) > 5 || abs(rayPos.z) > 5)
    //    discard;
    if (!isInsideBoundingBox(rayPos))
    {
        
        rayPos = rayAABBIntersection(rayPos, rayDir);
    }
  
    vec4 particleColor = vec4(0.0);
    if (rayPos != inf)
    {
        color = ((rayMarch(rayPos, rayDir)));
        //particleColor = vec4(0, 0.7, 0, 0.5);
    }
    
        
    //else 
    //    color = vec4(0, 0.0, 0.8, 0.5);

    color = color;//mix(color, particleColor, particleColor.a);
}