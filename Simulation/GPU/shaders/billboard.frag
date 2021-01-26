
#version 330


uniform vec3 lightDir;
uniform sampler2D sphereNormal;
uniform vec4 m_nup;
uniform vec4 m_nright;
uniform vec4 m_campos;
uniform float Particle_alpha;


in vec2 iuv;
in vec4 ipos;
out vec4 outColor;

void main()
{

    vec4 color = vec4(1.0,1.0,1.0,1);

    
    if (distance(iuv,vec2(0.5,0.5))<=0.47){
        
        outColor = vec4(color.xyz, Particle_alpha*(1-distance(iuv,vec2(0.5,0.5)))/0.47);
        
    }
    else
        discard;
}