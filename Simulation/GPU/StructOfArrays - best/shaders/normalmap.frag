
#version 330


uniform vec3 lightDir;
uniform sampler2D sphereNormal;
uniform vec4 m_nup;
uniform vec4 m_nright;
uniform vec4 m_campos;
uniform float Particle_alpha;
//in vec3 normal;

in vec2 iuv;
in vec4 ipos;
in vec4 center;
in vec4 icenter;
in mat3 rotation;

out vec4 outColor;
void main()
{
    vec4 color = vec4(0.14,1,1.0,1);

    vec3 ldir   = normalize(-lightDir);
    vec3 ambient = color.xyz * 0.2; // ambient
    vec3 diffuse = color.xyz * 0.8; // diffuse
    
    if (distance(iuv,vec2(0.5,0.5))<=0.47){
        vec3 normal = -1*texture(sphereNormal,vec2(1,1)*iuv.yx).xyz;
        
        normal = (rotation*normal);
        
        outColor = vec4(normalize(normal.xyz),Particle_alpha*(1-distance(iuv,vec2(0.5,0.5)))/0.47);
    }
    else
        discard;
}