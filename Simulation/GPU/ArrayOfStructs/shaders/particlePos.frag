
#version 330


uniform vec3 lightDir;
uniform sampler2D sphereNormal;
uniform vec4 m_nup;
uniform vec4 m_nright;
uniform vec4 m_campos;
//in vec3 normal;


in vec2 iuv;
in vec4 ipos;
in vec4 icenter;
in mat3 rotation;
out vec4 outColor;

void main()
{
    vec4 color = vec4(0.5,0.5,1.0,1);

    if (distance(iuv,vec2(0.5,0.5))<=0.47){

        vec3 normal = -1*texture(sphereNormal,vec2(1,1)*iuv.yx).xyz;
        normal = (rotation*normal);
        //vec3 pos = icenter.xyz+normal*0.01;
        outColor.xyz = normalize(ipos.xyz-m_campos.xyz);
        outColor.a = 1;
    }
    else
        discard;
}