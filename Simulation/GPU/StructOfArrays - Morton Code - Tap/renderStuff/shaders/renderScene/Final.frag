#version 130

uniform sampler2D color;
uniform sampler2D depth;

in vec2 texCoordV;

out vec4 outColor;

void main() {


	outColor = texture(depth, texCoordV);

}