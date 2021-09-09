#version 430

uniform sampler2D depthMap;

uniform int Blur;
uniform float Blur_Size;
uniform float Blur_Quality;
uniform vec2 resolution;
in vec2 uv;

layout(location = 0) out vec4 blur_depth;

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

#define EPS 1e-5

float linearize_depth(float d, float zNear, float zFar)
{
	return 0.5 + 0.5 * zNear * zFar / (zFar + d * (zNear - zFar));
}

float lum(in vec4 color)
{
	return linearize_depth(color.z, 0.02, 100.0);
}

/*

vec3 centreColour = texture2D(image,texCoord);
vec3 result = centreColour;
float normalization = 1;
for (int i = 0; i < NUM_SAMPLES; ++1)
{
   vec3 sample = texture2D(image, texCoord + offsets[i]);
   float gaussianCoeff = computeGaussianCoeff(offsets[i]);

   //depends on your implementation, this is quick
   float closeness = 1.0f - distance(sample,centreColour) / length(vec3(1,1,1));

   float sampleWeight = closeness * gaussian;

   result += sample * sampleWeight;
   noralization += sampleWeight;
}
vec3 bilateral = result / normalization;


*/

uniform int filterRadius;
uniform float sphereRadius;
uniform float bilateralParameter;

vec2 linearize_depth(vec2 d)
{
	float zNear = 0.02;
	float zFar = 100;
	return zNear * zFar / (zFar + d * (zNear - zFar));
}

float linearize_depth(float d)
{
	float zNear = 0.02;
	float zFar = 100;
	return zNear * zFar / (zFar + d * (zNear - zFar));
}

vec2 unlinearize_depth(vec2 d)
{
	float zNear = 0.02;
	float zFar = 100;
	return (-zFar + (zNear * zFar / d)) / (zNear - zFar);
}

float unlinearize_depth(float d)
{
	float zNear = 0.02;
	float zFar = 100;
	return (-zFar + (zNear * zFar / d)) / (zNear - zFar);
}

vec2 bilateral_new(sampler2D img)
{

	vec2 texCoord = uv;
	vec2 depth = linearize_depth(textureOffset(img, texCoord, ivec2(0)).xy);

	vec2 sum = vec2(0.f);
	vec2 wsum = vec2(0.f);

	for (int y = -filterRadius; y <= filterRadius; y += 1)
	{
		for (int x = -filterRadius; x <= filterRadius; x += 1)
		{

			vec2 p = vec2(float(x), float(y)) / resolution + texCoord;
			vec2 value = linearize_depth(textureOffset(img, texCoord, ivec2(x, y)).xy);

			vec2 d = (vec2(float(x), float(y))) / filterRadius;
			float sqrDist = dot(d, d);
			float w = (-sqrDist / (0.5f * 0.5f));
			vec2 r;
			//if (depth > 0.1)
			r = (value - depth) / (sphereRadius * bilateralParameter);

			//else
			// r = (value - depth) / (0.01*sphereRadius * bilateralParameter);

			vec2 g = -(r * r) / (0.5f);
			vec2 weight = exp(w + g);
			sum += (value * weight);
			wsum += weight;
		}
	}

	if (wsum.x > 0.f)
	{
		sum.x /= wsum.x;
	}
	if (wsum.y > 0.f)
	{
		sum.y /= wsum.y;
	}

	return unlinearize_depth(sum);
}

vec2 gaussian_new(sampler2D img)
{

	vec2 texCoord = uv;
	vec2 depth = linearize_depth(textureOffset(img, texCoord, ivec2(0)).xy);

	vec2 sum = vec2(0.f);
	vec2 wsum = vec2(0.f);

	for (int y = -filterRadius; y <= filterRadius; y += 1)
	{
		for (int x = -filterRadius; x <= filterRadius; x += 1)
		{

			vec2 p = vec2(float(x), float(y)) / resolution + texCoord;
			vec2 value = linearize_depth(textureOffset(img, texCoord, ivec2(x, y)).xy);

			vec2 d = (vec2(float(x), float(y))) / filterRadius;
			float sqrDist = dot(d, d);
			float w = (-sqrDist / (0.5f * 0.5f));

			vec2 weight = exp(vec2(w));
			sum += (value * weight);
			wsum += weight;
		}
	}

	if (wsum.x > 0.f)
	{
		sum.x /= wsum.x;
	}
	if (wsum.y > 0.f)
	{
		sum.y /= wsum.y;
	}

	return unlinearize_depth(sum);
}

vec2 box(sampler2D img)
{

	vec2 texCoord = uv;
	vec2 depth = linearize_depth(textureOffset(img, texCoord, ivec2(0)).xy);

	vec2 sum = vec2(0.f);
	vec2 wsum = vec2(0.f);

	for (int y = -filterRadius; y <= filterRadius; y += 1)
	{
		for (int x = -filterRadius; x <= filterRadius; x += 1)
		{

			vec2 p = vec2(float(x), float(y)) / resolution + texCoord;
			vec2 value = linearize_depth(textureOffset(img, texCoord, ivec2(x, y)).xy);

			sum += (value );
			wsum += 1;
		}
	}

	if (wsum.x > 0.f)
	{
		sum.x /= wsum.x;
	}
	if (wsum.y > 0.f)
	{
		sum.y /= wsum.y;
	}

	return unlinearize_depth(sum);
}

uniform float narrow_threshold;
uniform float narrow_threshold_far;
uniform float narrow_threshold_close;

vec2 narrow_gaussian(sampler2D img)
{

	vec2 texCoord = uv;
	vec2 depth = linearize_depth(textureOffset(img, texCoord, ivec2(0)).xy);

	vec2 sum = vec2(0.f);
	vec2 wsum = vec2(0.f);
	float n = 1;
	vec2 upper = depth.xy + vec2(mix(narrow_threshold_close, narrow_threshold_far, min(n, depth.x) / n),
								 mix(narrow_threshold_close, narrow_threshold_far, min(n, depth.y) / n));
	vec2 lower = depth.xy - vec2(mix(narrow_threshold_close, narrow_threshold_far, min(n, depth.x) / n),
								 mix(narrow_threshold_close, narrow_threshold_far, min(n, depth.y) / n));

	for (int y = -filterRadius; y <= filterRadius; y += 1)
	{
		for (int x = -filterRadius; x <= filterRadius; x += 1)
		{

			vec2 p = vec2(float(x), float(y)) / resolution + texCoord;
			vec2 value = linearize_depth(textureOffset(img, texCoord, ivec2(x, y)).xy);

			//value.xy= clamp(value, lower, upper );

			vec2 d = (vec2(float(x), float(y))) / filterRadius;
			float sqrDist = dot(d, d);
			float w = (-sqrDist / (0.5f * 0.5f));
			//vec2 r;
			//if (depth > 0.1)
			//r = (value - depth) / (sphereRadius * bilateralParameter);

			//else
			// r = (value - depth) / (0.01*sphereRadius * bilateralParameter);

			//vec2 g = -(r * r) / (0.5f);
			vec2 weight = exp(vec2(w) );
			if (value.x > lower.x && value.x < upper.x)
			{
				sum.x += (value.x * weight.x);
				wsum.x += weight.x;
			}
			if (value.y > lower.y && value.y < upper.y)
			{
				sum.y += (value.y * weight.y);
				wsum.y += weight.y;
			}
		}
	}

	if (wsum.x > 0.f)
	{
		sum.x /= wsum.x;
	}
	if (wsum.y > 0.f)
	{
		sum.y /= wsum.y;
	}

	return unlinearize_depth(sum);
}

uniform mat4 PROJECTION_MATRIX;
uniform mat4 PROJECTION_INVERSE_MATRIX;

uniform float dt;
uniform float dzt;

// Mean curvature. From "Screen Space Fluid Rendering with Curvature Flow"
vec2 meanCurvature(sampler2D image)
{
	vec2 pos = uv;
	// Width of one pixel
	vec2 dx = vec2(1.0f / resolution.x, 0.0f);
	vec2 dy = vec2(0.0f, 1.0f / resolution.y);

	// Central z value
	vec2 zc = texture(image, pos).xy;
	vec2 lzc = linearize_depth(zc);
	// Take finite differences
	// Central differences give better results than one-sided here.
	// TODO better boundary conditions, possibly.
	// Remark: This is not easy, get to choose between bad oblique view smoothing
	// or merging of unrelated particles
	vec2 zdxp = texture(image, pos + dx).xy;
	vec2 zdxn = texture(image, pos - dx).xy;
	vec2 zdx = 0.5f * (zdxp - zdxn);
	zdx.x = (zdxp.x == 0.0f || zdxn.x == 0.0f) ? 0.0f : zdx.x;
	zdx.y = (zdxp.y == 0.0f || zdxn.y == 0.0f) ? 0.0f : zdx.y;

	vec2 zdyp = texture(image, pos + dy).xy;
	vec2 zdyn = texture(image, pos - dy).xy;
	vec2 zdy = 0.5f * (zdyp - zdyn);

	zdy.x = (zdyp.x == 0.0f || zdyn.x == 0.0f) ? 0.0f : zdy.x;
	zdy.y = (zdyp.y == 0.0f || zdyn.y == 0.0f) ? 0.0f : zdy.y;
	
		// Take second order finite differences
	vec2 zdx2 = zdxp + zdxn - 2.0f * zc;
	vec2 zdy2 = zdyp + zdyn - 2.0f * zc;

	// Second order finite differences, alternating variables
	vec2 zdxpyp = texture(image, pos + dx + dy).xy;
	vec2 zdxnyn = texture(image, pos - dx - dy).xy;
	vec2 zdxpyn = texture(image, pos + dx - dy).xy;
	vec2 zdxnyp = texture(image, pos - dx + dy).xy;
	vec2 zdxy = (zdxpyp + zdxnyn - zdxpyn - zdxnyp) / 4.0f;

	// Projection transform inversion terms
	float cx = 2.0f / (resolution.x * -PROJECTION_MATRIX[0][0]);
	float cy = 2.0f / (resolution.y * -PROJECTION_MATRIX[1][1]);

	// Normalization term
	vec2 d = cy * cy * zdx * zdx + cx * cx * zdy * zdy + cx * cx * cy * cy * zc * zc;

	// Derivatives of said term
	vec2 ddx = cy * cy * 2.0f * zdx * zdx2 + cx * cx * 2.0f * zdy * zdxy + cx * cx * cy * cy * 2.0f * zc * zdx;
	vec2 ddy = cy * cy * 2.0f * zdx * zdxy + cx * cx * 2.0f * zdy * zdy2 + cx * cx * cy * cy * 2.0f * zc * zdy;

	// Temporary variables to calculate mean curvature
	vec2 ex = 0.5f * zdx * ddx - zdx2 * d;
	vec2 ey = 0.5f * zdy * ddy - zdy2 * d;

	// Finally, mean curvature
	vec2 h;
	h.x = 0.5f * ((cy * ex.x + cx * ey.x) / pow(d.x, 1.5f));
	h.y = 0.5f * ((cy * ex.y + cx * ey.y) / pow(d.y, 1.5f));

	//vec3 dxyz = vec3(zdx, zdy, h);

	// Vary contribution with absolute depth differential - trick from pySPH
	vec2 outDepth = (zc) + (h * dt * (1.0f + (abs(zdx) + abs(zdy)) * dzt));

	return clamp(outDepth,0,1);

}

// Mean curvature. From "Screen Space Fluid Rendering with Curvature Flow"
vec2 narrowMeanCurvature(sampler2D image)
{
	vec2 pos = uv;
	// Width of one pixel
	vec2 dx = vec2(1.0f / resolution.x, 0.0f);
	vec2 dy = vec2(0.0f, 1.0f / resolution.y);

	// Central z value
	vec2 zc = texture(image, pos).xy;
	vec2 lzc = linearize_depth(zc);

	
	// Take finite differences
	// Central differences give better results than one-sided here.
	// TODO better boundary conditions, possibly.
	// Remark: This is not easy, get to choose between bad oblique view smoothing
	// or merging of unrelated particles
	vec2 zdxp = texture(image, pos + dx).xy;
	vec2 zdxn = texture(image, pos - dx).xy;
	vec2 zdx = 0.5f * (zdxp - zdxn);
	if (    abs(linearize_depth(zdxp.x)-lzc.x) > narrow_threshold ||  abs(linearize_depth(zdxn.x)-lzc.x) > narrow_threshold ||
			abs(linearize_depth(zdxp.y)-lzc.y) > narrow_threshold ||  abs(linearize_depth(zdxn.y)-lzc.y) > narrow_threshold )
				return zc;
	zdx.x = (zdxp.x == 0.0f || zdxn.x == 0.0f) ? 0.0f : zdx.x;
	zdx.y = (zdxp.y == 0.0f || zdxn.y == 0.0f) ? 0.0f : zdx.y;

	vec2 zdyp = texture(image, pos + dy).xy;
	vec2 zdyn = texture(image, pos - dy).xy;
	vec2 zdy = 0.5f * (zdyp - zdyn);

	zdy.x = (zdyp.x == 0.0f || zdyn.x == 0.0f) ? 0.0f : zdy.x;
	zdy.y = (zdyp.y == 0.0f || zdyn.y == 0.0f) ? 0.0f : zdy.y;
	
	if  (  abs(linearize_depth(zdyp.x)-lzc.x) > narrow_threshold ||  abs(linearize_depth(zdyn.x)-lzc.x) > narrow_threshold 
		|| abs(linearize_depth(zdyp.y)-lzc.y) > narrow_threshold ||  abs(linearize_depth(zdyn.y)-lzc.y) > narrow_threshold)
			return zc;
	

	// Take second order finite differences
	vec2 zdx2 = zdxp + zdxn - 2.0f * zc;
	vec2 zdy2 = zdyp + zdyn - 2.0f * zc;

	// Second order finite differences, alternating variables
	vec2 zdxpyp = texture(image, pos + dx + dy).xy;
	vec2 zdxnyn = texture(image, pos - dx - dy).xy;
	vec2 zdxpyn = texture(image, pos + dx - dy).xy;
	vec2 zdxnyp = texture(image, pos - dx + dy).xy;
	vec2 zdxy = (zdxpyp + zdxnyn - zdxpyn - zdxnyp) / 4.0f;
	if (abs(linearize_depth(zdxpyp.x)-lzc.x) > narrow_threshold || abs(linearize_depth(zdxpyp.y)-lzc.y) > narrow_threshold ||
		abs(linearize_depth(zdxnyn.x)-lzc.x) > narrow_threshold || abs(linearize_depth(zdxnyn.y)-lzc.y) > narrow_threshold ||
		abs(linearize_depth(zdxpyn.x)-lzc.x) > narrow_threshold || abs(linearize_depth(zdxpyn.y)-lzc.y) > narrow_threshold ||
		abs(linearize_depth(zdxnyp.x)-lzc.x) > narrow_threshold || abs(linearize_depth(zdxnyp.y)-lzc.y) > narrow_threshold )
		return zc;


	// Projection transform inversion terms
	float cx = 2.0f / (resolution.x * -PROJECTION_MATRIX[0][0]);
	float cy = 2.0f / (resolution.y * -PROJECTION_MATRIX[1][1]);

	// Normalization term
	vec2 d = cy * cy * zdx * zdx + cx * cx * zdy * zdy + cx * cx * cy * cy * zc * zc;

	// Derivatives of said term
	vec2 ddx = cy * cy * 2.0f * zdx * zdx2 + cx * cx * 2.0f * zdy * zdxy + cx * cx * cy * cy * 2.0f * zc * zdx;
	vec2 ddy = cy * cy * 2.0f * zdx * zdxy + cx * cx * 2.0f * zdy * zdy2 + cx * cx * cy * cy * 2.0f * zc * zdy;

	// Temporary variables to calculate mean curvature
	vec2 ex = 0.5f * zdx * ddx - zdx2 * d;
	vec2 ey = 0.5f * zdy * ddy - zdy2 * d;

	// Finally, mean curvature
	vec2 h;
	h.x = 0.5f * ((cy * ex.x + cx * ey.x) / pow(d.x, 1.5f));
	h.y = 0.5f * ((cy * ex.y + cx * ey.y) / pow(d.y, 1.5f));

	//vec3 dxyz = vec3(zdx, zdy, h);

	// Vary contribution with absolute depth differential - trick from pySPH
	vec2 outDepth = (zc) + (h * dt * (1.0f + (abs(zdx) + abs(zdy)) * dzt));

	return clamp(outDepth,0,1);
}
void main()
{
	if (texture(depthMap, uv).z == 0)
	{
		discard;
	}
	if (Blur > 0)
	{
		//if(Blur == 1){
		//
		//  blur_depth_front = smooth_Custom(frontdepthMap);
		//  blur_depth_front.y = blur_depth_back.x;
		//  blur_depth_front.z = 1;
		//
		//
		//}
		//else
		if (Blur == 1)
		{
			blur_depth.xy = gaussian_new(depthMap).xy;
			blur_depth.za = vec2(1.0);
		}
		else if (Blur == 2)
		{
			blur_depth.xy = bilateral_new(depthMap).xy;
			blur_depth.za = vec2(1.0);
		}
		else if (Blur == 3)
		{
			blur_depth.xy = (narrow_gaussian(depthMap));
			blur_depth.za = vec2(1.0);
		}
		else if (Blur == 4)
		{
			blur_depth.xy = (meanCurvature(depthMap));
			blur_depth.za = vec2(1.0);
		}
		else if (Blur == 5)
		{
			blur_depth.xy = (narrowMeanCurvature(depthMap));
			blur_depth.za = vec2(1.0);
		}else if (Blur == 6){
			blur_depth.xy = box(depthMap).xy;
			blur_depth.za = vec2(1.0);
		}

	}
	else
	{
		blur_depth = texture(depthMap, uv); //blur_depth_back = texture(backdepthMap,uv);
		blur_depth.za = vec2(1.0);
	}
}