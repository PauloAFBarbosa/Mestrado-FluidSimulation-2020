#version 430



uniform vec2  view_size;

uniform mat4 m_model;
uniform vec4 m_nup;
uniform vec4 m_nright;
uniform mat4 m_pv;
uniform mat3 normalMatrix;
uniform vec4 m_campos;
uniform vec3 Particle_Size;

in vec2 iuv;
in vec4 ipos;
in vec4 icenter;
in mat3 rotation;

layout (location = 0) out vec4 depth;
//out vec4 outColor;
out float gl_FragDepth;
#define PI 3.1415926538

float linearize_depth(float d,float zNear,float zFar)
{
    return zNear * zFar / (zFar + d * (zNear - zFar));
}

void main()
{
    if (distance(iuv,vec2(0.5,0.5))<=0.5){
        float a = PI*iuv.x;
        float b = PI*iuv.y;
        //float y = sin(a)*cos(b);//sen a cos b
        float z = sin(a)*sin(b);//sen a sen b
        //float x = cos(a);//cos a
        
        vec3 n = -1*vec3(0,0,z); //vec3(iuv.x-0.5,iuv.y-0.5,0.5-distance(iuv,vec2(0.5,0.5)));//texture(sphereNormal,vec2(1,1)*iuv).xyz;
        

        vec4 normal;
        normal.xyz = ((rotation*n)) ;
        normal.a   = 0;

        //*Front-----------------
        vec4 pos;
        pos.xyz = ipos.xyz - normal.xyz;// icenter.xyz-(normal.xyz);
        pos.a = 1;

        vec4 clip_space_pos = ((m_pv*pos));
        float ndc_depth =clamp(clip_space_pos.z / clip_space_pos.w,0,1);
    
        pos.xyz = icenter.xyz+(normal.xyz);
        pos.a = 1;
        clip_space_pos = ((m_pv*pos));
        float far_ndc_depth  = clip_space_pos.z / clip_space_pos.w;


        
        depth   = vec4(vec2(1-ndc_depth,far_ndc_depth),1,1);
    }
    else
        discard;
}