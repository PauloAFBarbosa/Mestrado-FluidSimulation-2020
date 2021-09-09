#version 330

in vec2 uv;

uniform sampler2D depth;

out vec4 color;

void main() {
    color = texture(depth,uv);
    if(color.z == 1)
        color.x = 1-color.x;
    else{
        color.x = 1;
        color.y = 1;
        color.z = 0;
    } 
    
}