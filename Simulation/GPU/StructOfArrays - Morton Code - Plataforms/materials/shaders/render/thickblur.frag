#version 430

uniform sampler2D depthMap;
uniform sampler2D thickness;

uniform int Blur;
uniform float Blur_Size;
uniform float Blur_Quality;
uniform vec2 resolution;
in vec2 uv;


layout(location = 0) out vec4 blur_thick;

vec4 blur(sampler2D image, vec2 uv, vec2 resolution)
{

	float Pi = 6.28318530718; // Pi*2

	// GAUSSIAN BLUR SETTINGS {{{
	float Directions = 16.0;	  // BLUR DIRECTIONS (Default 16.0 - More is better but slower)
	float Quality = Blur_Quality; // BLUR QUALITY (Default 4.0 - More is better but slower)
	float Size = Blur_Size;		  // BLUR SIZE (Radius)
	// GAUSSIAN BLUR SETTINGS }}}

	vec2 Radius = Size / resolution.xy;

	// Pixel colour
	vec4 Color = texture(image, uv);

	float k = 0;
	for (float d = 0.0; d < Pi; d += Pi / Directions)
	{
		for (float i = 1.0 / Quality; i <= 1.0; i += 1.0 / Quality)
		{
			vec2 coords = uv + vec2(cos(d), sin(d)) * Radius * i;
			if (coords.x >= 0 && coords.x <= 1 && coords.y >= 0 && coords.y <= 1)
			{
				Color += texture(image, coords);
				k += 1.0;
			}
		}
	}

	return Color / k;
}

vec4 blur5(sampler2D image, vec2 direction)
{
	vec4 color = vec4(0.0);
	vec2 off1 = vec2(1.3333333333333333) * direction;
	color += texture(image, uv) * 0.29411764705882354;
	color += texture(image, uv + (off1 / resolution)) * 0.35294117647058826;
	color += texture(image, uv - (off1 / resolution)) * 0.35294117647058826;
	return color;
}

vec4 blur13(sampler2D image, vec2 direction)
{
	vec4 color = vec4(0.0);
	vec2 off1 = vec2(1.411764705882353) * direction;
	vec2 off2 = vec2(3.2941176470588234) * direction;
	vec2 off3 = vec2(5.176470588235294) * direction;
	color += texture(image, uv) * 0.1964825501511404;
	color += texture(image, uv + (off1 / resolution)) * 0.2969069646728344;
	color += texture(image, uv - (off1 / resolution)) * 0.2969069646728344;
	color += texture(image, uv + (off2 / resolution)) * 0.09447039785044732;
	color += texture(image, uv - (off2 / resolution)) * 0.09447039785044732;
	color += texture(image, uv + (off3 / resolution)) * 0.010381362401148057;
	color += texture(image, uv - (off3 / resolution)) * 0.010381362401148057;
	return color;
}

vec4 smooth_Custom(sampler2D image)
{
	//if(Blur == 1)
	//  return (blur13(image,vec2(1,0)) + blur13(image,vec2(0,1)) + blur13(image,vec2(1,-1)) + blur13(image,vec2(1,1)))/4;
	return blur(image, uv, resolution);
	//return texture(image,uv);
}



uniform int filterRadius;

float gaussian_new(sampler2D img)
{

	vec2 texCoord = uv;
	float depth = textureOffset(img, texCoord, ivec2(0)).x;

	float sum = 0;
	float wsum = 0;

	for (int y = -filterRadius; y <= filterRadius; y += 1)
	{
		for (int x = -filterRadius; x <= filterRadius; x += 1)
		{

			vec2 p = vec2(float(x), float(y)) / resolution + texCoord;

			vec4 value = textureOffset(img, texCoord, ivec2(x, y));

			vec2 d = (vec2(float(x), float(y))) / filterRadius;
			float sqrDist = dot(d, d);
			float w = (-sqrDist / (0.5f * 0.5f));

			float weight = exp(w);
			sum += (value.x * weight);
			wsum += weight;
		}
	}

	if (wsum > 0.f)
	{
		sum/= wsum;
	}
	

	return sum;
}

void main()
{
	if (texture(depthMap, uv).z == 0)
	{
		discard;
	}
	if (Blur > 0)
	{
		blur_thick = vec4(vec3(gaussian_new(thickness)), 1);
	}
	else
	{
		blur_thick = texture(thickness, uv);
	}
}