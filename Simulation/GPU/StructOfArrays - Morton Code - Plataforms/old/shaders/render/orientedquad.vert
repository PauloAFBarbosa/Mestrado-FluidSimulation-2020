#version 430

uniform	mat4 PROJECTION_VIEW;

layout(std430, binding = 0) buffer Position{
    vec4 testpos[];
};

uniform mat4 m_model;
uniform vec4 m_nup;
uniform vec4 m_nright;
uniform mat4 m_pv;
uniform mat3 normalMatrix;
uniform vec4 m_campos;
uniform vec3 Particle_Size;

in vec4 position;

out vec4 oposition;

out vec2 iuv;
out vec4 ipos;
out vec4 icenter;
out mat3 rotation;

void main()
{
    vec4 pos = vec4(0);
    int particleId= gl_VertexID/6;
    if(gl_VertexID%6 == 0){
        iuv = vec2(0,0);
        pos = vec4(-1,-1,0,1);
    }
    if(gl_VertexID%6 == 1){
        iuv = vec2(0,1);
        pos = vec4(1,-1,0,1);
    }
    if(gl_VertexID%6 == 2){
        iuv = vec2(1,1);
        pos = vec4(1,1,0,1);
    }
    if(gl_VertexID%6 == 3){
        iuv = vec2(1,0);
        pos = vec4(-1,1,0,1);
    }
    if(gl_VertexID%6 == 4){
        iuv = vec2(0,0);
        pos = vec4(-1,-1,0,1);
    }
    if(gl_VertexID%6 == 5){
        iuv = vec2(1,1);
        pos = vec4(1,1,0,1);
    }
    icenter = testpos[particleId];

    //vec3 forward3 = normalize(m_campos.xyz-icenter.xyz);
    //vec3 right3 =   ((normalize(cross(m_nup.xyz,forward3))));
    //vec3 up3 =      normalize(cross(forward3,m_nright.xyz));
    
    vec3 forward3 =  normalize(cross(m_nright.xyz,m_nup.xyz));
    vec3 right3 =   m_nright.xyz;// ((normalize(cross(m_nup.xyz,forward3))));
    vec3 up3 =      m_nup.xyz;//normalize(cross(forward3,m_nright.xyz));
    mat3 w;

    w =  (mat3(right3,up3,forward3)) * mat3( vec3(Particle_Size.x,0,0),vec3(0,Particle_Size.y,0),vec3(0,0,Particle_Size.z)) ;
    
    rotation = w;

    pos.xyz = w*pos.xyz;
    pos.xyz += icenter.xyz;
    ipos = pos;
    //pos.x += (trianglenum/1000)/100.0;
    //pos.z -= (trianglenum%1000)/100.0;
    gl_Position = (PROJECTION_VIEW * pos).xyzw;
    oposition = pos;
}