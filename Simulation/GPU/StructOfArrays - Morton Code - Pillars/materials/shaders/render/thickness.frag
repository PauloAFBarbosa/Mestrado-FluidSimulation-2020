#version 430



uniform vec2  view_size;

uniform mat4 m_model;
uniform vec4 m_nup;
uniform vec4 m_nright;
uniform mat4 m_pv;
uniform mat3 normalMatrix;
uniform vec4 m_campos;
uniform vec3 Particle_Size;
uniform float Particle_alpha;

in vec2 iuv;
in vec4 ipos;
in vec4 icenter;
in mat3 rotation;


out vec4 outColor;
#define PI 3.1415926538

float linearize_depth(float d,float zNear,float zFar)
{
    return zNear * zFar / (zFar + d * (zNear - zFar));
}

void main()
{
    vec4 color = vec4(1.0,1.0,1.0,1);

    
    if (distance(iuv,vec2(0.5,0.5))<=0.5){
        float a = PI*iuv.x;
        float b = PI*iuv.y;
        float z = sin(a)*sin(b);
        outColor = vec4(color.xyz, 2*z*Particle_Size.z);
    }
    else
        discard;
}