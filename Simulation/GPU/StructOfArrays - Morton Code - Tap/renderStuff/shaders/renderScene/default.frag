#version 430


uniform vec4 diffuse;
uniform vec4 specular, emission;
uniform float shininess;
uniform int texCount;
uniform sampler2D texUnit;


in vec3 normalV;
in vec3 eyeV;
in vec3 lightDirV;
in vec2 texCoordV;
in vec4 pos;

layout (location = 0) out vec4 colorOut;
layout (location = 1) out vec4 depthOut;

void main() {

	vec4 dif;
	vec4 spec;

	float intensity = max(dot(normalize(normalV),lightDirV), 0.0);


	vec3 h = normalize(lightDirV + normalize(eyeV));	
	float intSpec = max(dot(h,normalize(normalV)), 0.0);
	spec = specular * pow(intSpec,shininess);
	dif = diffuse;

	colorOut = vec4(intensity + 0.3);
	depthOut = pos;
	//if (texCount == 0) {
	//	colorOut = (intensity + 0.3) * dif + emission + spec ;
	//}
	//else {
	//	colorOut = ((intensity + 0.3) * dif + emission) * texture(texUnit, texCoordV) + spec;
	//}
	
	//colorOut = (intensity + 0.2) * dif + spec ;
	
}