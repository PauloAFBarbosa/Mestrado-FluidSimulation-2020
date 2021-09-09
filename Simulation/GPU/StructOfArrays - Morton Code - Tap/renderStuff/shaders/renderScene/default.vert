#version 430


// uniform mat3 NormalMatrix;

uniform mat4 PVM, VM;

uniform vec4 lightDirection;

in vec4 position;
in vec4 normal;
in vec4 texCoord0;


out vec3 normalV;
out vec2 texCoordV;
out vec3 eyeV;
out vec3 lightDirV;
out vec4 pos;

void main () {
	
	texCoordV = vec2(texCoord0);

	normalV = normalize(vec3(normal));
	lightDirV = normalize(-vec3(lightDirection));

	pos = position;
	vec3 p = vec3(VM * position);
	eyeV = normalize(-p);

	gl_Position = PVM * position;	
	
}