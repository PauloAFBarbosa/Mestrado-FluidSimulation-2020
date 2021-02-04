#version 430


uniform sampler2D backgroundTex;
uniform sampler2D normalMapSphere;
uniform sampler2D fluidColor;
uniform sampler2D spherePosTex;
uniform samplerCube texUnit; 

uniform vec3 lightDir;
uniform vec4 Particle_Color;
uniform float eta;
uniform float Alpha_Threshold;
uniform int Blur;
uniform int WhiteParticles;
in vec2 uv;

out vec4 color;
out float gl_FragDepth ;

uniform float Blur_Size;
uniform float Blur_Quality;



vec4 blur(sampler2D image, vec2 uv, vec2 resolution) {

  float Pi = 6.28318530718; // Pi*2
    
  // GAUSSIAN BLUR SETTINGS {{{
  float Directions = 16.0; // BLUR DIRECTIONS (Default 16.0 - More is better but slower)
  float Quality = Blur_Quality; // BLUR QUALITY (Default 4.0 - More is better but slower)
  float Size = Blur_Size; // BLUR SIZE (Radius)
  // GAUSSIAN BLUR SETTINGS }}}

  vec2 Radius = Size/resolution.xy;
    
  // Pixel colour
  vec4 Color = texture(image, uv);
    
  float k = 0;
  for( float d=0.0; d<Pi; d+=Pi/Directions)
    {
		for(float i=1.0/Quality; i<=1.0; i+=1.0/Quality)
        {
          vec2 coords = uv+vec2(cos(d),sin(d))*Radius*i;
          if (coords.x >= 0 && coords.x <= 1 && coords.y >= 0 && coords.y <= 1){
            Color += texture( image, coords);		
            k+=1.0;
          }
        }
    }

  return Color/k;
}

vec4 smooth_Custom(sampler2D image,vec2 uv, vec2 resolution){
  if(Blur == 1)
    return blur (image,uv,resolution);
  return texture(image,uv);
}

void main()
{
    color = texture(backgroundTex, uv);
    vec2 resolution = vec2(1920,1080);
    vec4 particleColorBlured =  smooth_Custom(fluidColor,uv,resolution );   
    vec4 normal = normalize((smooth_Custom(normalMapSphere,uv,resolution)));
    vec4 lookAtVec = normalize((smooth_Custom(spherePosTex,uv,resolution)));

    if(WhiteParticles == 1){
      vec4 particleColor =  texture( fluidColor, uv);
      vec3 n = normalize(-1*normal.xyz);
      vec3 ldir   = normalize(-lightDir);

      float diff = max(dot(n, ldir), 0.0);
      if(particleColor.x > 0)
        color = vec4(vec3(diff),1);
      return;
    }

    if( particleColorBlured.x >  Alpha_Threshold){
      vec3 n = normalize(-1*normal.xyz);
	    vec3 e = normalize(lookAtVec.xyz);
	    vec3 t = normalize(reflect(e, n));
	    vec3 ref = texture(texUnit, (-1*t)).rgb;
  
  
	    vec3 r = refract(e, n, eta);
	    float k = 1.0 - eta * eta * (1.0 - dot(n, e) * dot(n, e));
	    vec3 refrat;
      refrat = vec3(texture(texUnit, (-1*r)));
      
      vec3 ldir   = normalize(-lightDir);
      
      
      vec3 h = normalize(ldir - e);	
      // compute the specular intensity
      float intSpec = max(dot(h,n), 0.0);
      // compute the specular term into spec
      vec4 spec = vec4(1) * pow(intSpec,128);

      float fresnel = abs(dot(-normalize(e), normalize(n)));

      vec3 refRefac = vec3((k) * refrat + (1-k) * (ref));
      refRefac = mix(refRefac,refRefac*(Particle_Color.xyz),smoothstep(0,1,(particleColorBlured.x)*(Particle_Color.a)));
      color = vec4(
        ( refRefac)
        , 1.0);
    }
}