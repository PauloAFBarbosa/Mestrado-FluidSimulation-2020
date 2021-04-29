#version 430


uniform sampler2D backgroundTex;
uniform sampler2D fluidColor;

uniform sampler2D front_depth;
uniform sampler2D back_depth;
uniform sampler2D depthMap;

uniform samplerCube texUnit; 

uniform vec3 lightDir;
uniform vec4 Particle_Color;
uniform vec4 Thin_Color;
uniform float eta;
uniform int DebugParticles;
uniform vec4 m_campos;
uniform mat4 m_inv_proj;
uniform mat4 m_inv_view;
uniform vec2 resolution;
uniform int thicktype;

in vec2 uv;

out vec4 color;
out float gl_FragDepth ;

uniform float Particle_alpha;




vec4 getWorldPosBack(vec2 offset){
  vec2 coords = uv+(offset/resolution);
  vec4 clip_space_pos = vec4(coords * 2.0 - 1.0, (texture(back_depth,coords).z), 1.0);
  vec4 viewspace = m_inv_proj*clip_space_pos;
  vec4 worldspacepos = m_inv_view*((viewspace/viewspace.w));

  return worldspacepos; 
}



vec4 getWorldPosFront(vec2 offset){
  vec2 coords = uv+(offset/resolution);
  vec4 clip_space_pos = vec4(coords * 2.0 - 1.0, 1-texture(front_depth,coords).z, 1.0);
  vec4 viewspace = m_inv_proj*clip_space_pos;
  vec4 worldspacepos = m_inv_view*(viewspace/viewspace.w);

  return worldspacepos; 
}



vec4 getWorldPos(vec2 offset){

  //return getWorldPosBack(offset);
  //return getWorldPosFront(offset);
  vec2 coords = uv+(offset/resolution);
  vec4 clip_space_pos = vec4(coords * 2.0 - 1.0, 1-texture(depthMap,coords).z, 1.0);
  vec4 viewspace = m_inv_proj*clip_space_pos;
  vec4 worldspacepos = m_inv_view*(viewspace/viewspace.w);

  return worldspacepos; 
}


vec4 normalFromFrontDepth(){
  vec4 c = getWorldPosFront(vec2(0,0)) ;  
  vec4 t = getWorldPosFront(vec2(0,2)) ;  
  vec4 l = getWorldPosFront(vec2(2,0)) ;

  return -1*vec4(normalize(cross((t-c).xyz,(l-c).xyz)),0);
}
vec4 normalFromBackDepth(){
  vec4 c = getWorldPosBack(vec2(0,0)) ;  
  vec4 t = getWorldPosBack(vec2(0,2)) ;  
  vec4 l = getWorldPosBack(vec2(2,0)) ;

  return -1*vec4(normalize(cross((t-c).xyz,(l-c).xyz)),0);
}

vec4 normalFromDepth(){
  //float dzdx = (-0.5*(textureOffset(depthMap,uv,ivec2(1,0)) - textureOffset(depthMap,uv,ivec2(-1,0)))).z;
  //float dzdy = (-0.5*(textureOffset(depthMap,uv,ivec2(0,1)) - textureOffset(depthMap,uv,ivec2(0,-1)))).z;
  //vec4 viewspace = m_inv_proj*(vec4((normalize(vec3(dzdx,dzdy,1))),0));
  //return  m_inv_view*(viewspace/viewspace.w);
  //vec3 t = vec3(uv.x,uv.y,textureOffset());
  //vec3 l = vec3(uv.x,uv.y,textureOffset());
  vec4 c = getWorldPos(vec2(0,0)) ;  
  vec4 t = getWorldPos(vec2(0,2)) ;  
  vec4 l = getWorldPos(vec2(2,0)) ;

  return -1*vec4(normalize(cross((t-c).xyz,(l-c).xyz)),0);
}


vec3 iterate(){
  float j = 0;
  vec4 t = 1-texture(front_depth,uv + vec2( 1,0)/resolution);
  vec4 b = 1-texture(front_depth,uv + vec2(-1,0)/resolution);
  vec4 r = 1-texture(front_depth,uv + vec2(0, 1)/resolution);
  vec4 l = 1-texture(front_depth,uv + vec2(0,-1)/resolution);
  vec3 normal = vec3(normalize(cross((t-b).xyz,(r-l).xyz)));

  vec3 raydir = -1*normalize(refract(normal,vec3(0,0,1),eta));
  float stepLength = 0.01;
  vec3 rayStep = stepLength*raydir;
  vec3 rayPos=vec3(uv.xy,1-texture(front_depth,uv).x);
 //# if rayStep.z == 0
  for(  ; 
        rayPos.z < 1 && rayPos.z < texture(back_depth,rayPos.xy).x; 
        rayPos+= rayStep){
    //j+= texture(front_depth,uv + ivec2(i,0)/resolution).z;
  }
  return vec3(rayPos);
}



vec3 refratBackColor(vec3 lookAtVec, vec3 ldir){

  vec3 n = normalize(normalFromBackDepth().xyz); 
	vec3 e = normalize(lookAtVec.xyz);
	vec3 t = normalize(reflect(e, n));
	vec3 ref = texture(texUnit, (-1*t)).rgb;
  
  vec3 r = refract(e, n, eta);
  //float k = 1.0 - eta * eta * (1.0 - dot(n, e) * dot(n, e));
  vec3 refrat;
  refrat = vec3(texture(texUnit, (-1*r)));

  float fresnel = abs(dot(-normalize(e), normalize(n)));
  vec3 refRefac = vec3((fresnel) * refrat + (1-fresnel) * (ref));
  vec3 h = normalize(ldir - e);	
      // compute the specular intensity
  float intSpec = max(dot(h,n), 0.0);
      // compute the specular term into spec
  vec4 spec = vec4(1) * pow(intSpec,128);

  return refRefac + spec.xyz;
}


uniform int Surface_Num;

void main()
{
   // color = texture(spherePosTex,uv);  return;
    float depth = 1-texture(front_depth,uv).z;
    color = texture(backgroundTex, uv);
    //vec2 resolution = vec2(1920,1080);
    vec4 thickness;
    //thickness.xyz = iterate();
    if(thicktype == 1){
      thickness = vec4(distance(getWorldPosFront(vec2(0)).xyz,getWorldPosBack(vec2(0)).xyz));  
    }else if(thicktype == 2){
      thickness =  texture(fluidColor,uv)*Particle_alpha;   
    }else if(thicktype == 3){
      thickness = vec4(
        min(vec4(distance(getWorldPosFront(vec2(0)).xyz,getWorldPosBack(vec2(0)).xyz)),
              texture(fluidColor,uv).x*Particle_alpha)
        );

    }
    vec3 normal = normalFromDepth().xyz;//normalize((smooth_Custom(normalMapSphere,uv,resolution)));

    if(depth == 0)
      normal = normalFromBackDepth().xyz;
    vec4 worldpos = getWorldPos(vec2(0)); 

    vec4 lookAtVec = normalize(worldpos-m_campos);

      vec4 particleColor =  texture( fluidColor, uv);
      vec3 n = normalize(normal.xyz);
      vec3 ldir   = normalize(-lightDir);

      float diff = (0.2+0.8*max(dot(n, ldir), 0.0));
    if(DebugParticles != 0 ){
      if(depth != 1 ){
        if(DebugParticles == 1)
          color = vec4(vec3(diff)*Particle_Color.xyz,1);
        else if(DebugParticles == 2)
          color = vec4(vec3((worldpos)),1);
        else if(DebugParticles == 3)
          color = vec4(vec3((normal+vec3(1))/2),1);
        else if(DebugParticles == 4)
          color = vec4((thickness.xyz)*Particle_Color.xyz,1);
      } 
      else 
        color = vec4(0);
      return;
    }

    if( depth <  1){
      vec3 n = normalize(normal.xyz);
	    vec3 e = normalize(lookAtVec.xyz);
	    vec3 t = normalize(reflect(e, n));
	    vec3 ref = texture(texUnit, (-1*t)).rgb;
  
      //vec3 refratnormal = normalize(normalFromBackDepth().xyz); 
	   
	    //float k = 1.0 - eta * eta * (1.0 - dot(n, e) * dot(n, e));
      vec3 ldir   = normalize(-lightDir);

	    vec3 refrat;
      if (Surface_Num == 1){
        vec3 r = refract(e, n, eta);
	      refrat = vec3(texture(texUnit, (-1*r)));
      }
      else if(Surface_Num == 2){
	      refrat = refratBackColor(lookAtVec.xyz,ldir);// vec3(texture(texUnit, (-1*r)));
      }
      
      
      
      vec3 h = normalize(ldir - e);	
      // compute the specular intensity
      float intSpec = max(dot(h,n), 0.0);
      // compute the specular term into spec
      vec4 spec = vec4(1) * pow(intSpec,128);

      float fresnel = abs(dot(-normalize(e), normalize(n)));
      refrat = mix(refrat,refrat*(Particle_Color.xyz),smoothstep(0,1,(thickness.x)));
      vec3 refRefac = vec3((fresnel) * refrat + (1-fresnel) * (ref));
      refRefac = mix(refRefac,Thin_Color.xyz*refRefac,smoothstep(0,1,(Thin_Color.a)));
      //if(depth == 0)
      //  refRefac = refrat;
      color = vec4(
        ( refRefac + spec.xyz)
        , 1.0);
    }
}