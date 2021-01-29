#version 440
layout(points) in;
layout(triangle_strip, max_vertices=300) out;

uniform mat4 m_model;
uniform vec4 m_nup;
uniform vec4 m_nright;
uniform mat4 m_pv;
uniform mat3 normalMatrix;

uniform vec4 m_campos;

in data{
    vec4 pos;
    uint id;
    
}i[1];


out vec2 iuv;
out vec4 ipos;
out vec4 icenter;
out mat3 rotation;
uniform vec2 Particle_Size;

void GenerateVertex(vec4 p, vec2 uv)
{
  iuv = uv;
  gl_Position = m_pv * p;
  ipos = p;
  EmitVertex();
}

void GenerateSquare(vec4 circleY)
{
    vec4 center = circleY;
    vec4 centerInScreen = m_pv * center;
    icenter = center;
    vec3 forward3 = normalize(m_campos.xyz-icenter.xyz);
    vec3 right3 =   ((normalize(cross(m_nup.xyz,forward3))));
    vec3 up3 =      normalize(cross(forward3,m_nright.xyz));

    
    mat3 w =  mat3(vec3(Particle_Size.y/Particle_Size.x,0,0),vec3(0,Particle_Size.x/Particle_Size.y,0),vec3(0,0,1))*(mat3(right3,up3,forward3)) ;
    rotation = w;
    vec4 up    = Particle_Size.y*vec4(up3,0);
    vec4 right = Particle_Size.x*vec4(right3,0);  
    vec4 vert1 = vec4(center.x,center.y,center.z,1) - up - right;
    vec4 vert2 = vec4(center.x,center.y,center.z,1) + up - right;
    vec4 vert3 = vec4(center.x,center.y,center.z,1) - up + right;
    vec4 vert4 = vec4(center.x,center.y,center.z,1) + up + right;

    GenerateVertex(vert3,vec2(0,1));
    GenerateVertex(vert4,vec2(1,1));
    GenerateVertex(vert1,vec2(0,0));
    GenerateVertex(vert2,vec2(1,0));

    EndPrimitive();
}

vec3 gNormals[14] =
	{
		vec3(0.0f, 0.0f, -1.0f),
		vec3(0.0f, 0.0f, -1.0f),
		vec3(0.0f, 0.0f, -1.0f),
		vec3(0.0f, 0.0f, -1.0f),

		vec3(1.0f, 0.0f,  0.0f),
		vec3(1.0f, 0.0f,  0.0f),
		vec3(1.0f, 0.0f,  0.0f),

		vec3(-1.0f, 0.0f,  0.0f),
		vec3(-1.0f, 0.0f,  0.0f),
		vec3(-1.0f, 0.0f,  0.0f),
		vec3(-1.0f, 0.0f,  0.0f),

		vec3(0.0f, 0.0f,  1.0f),
		vec3(0.0f, 0.0f,  1.0f),
		vec3(0.0f, 0.0f,  1.0f)
	};








void main()
{ 
    vec4 p = i[0].pos;

    GenerateSquare(p);
}