#version 440


uniform mat4 m_model;

in vec3 position;


out data{
    vec4 pos;
    uint id;
}o;

void main()
{
    o.id = gl_VertexID;
    uint layer = o.id%100;
    o.pos = vec4(position,1);
}