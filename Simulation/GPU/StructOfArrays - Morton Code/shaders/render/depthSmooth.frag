#version 430


uniform sampler2D depthMap;
uniform sampler2D frontdepthMap;
uniform sampler2D backdepthMap;
uniform sampler2D thickness;

uniform int Blur;
uniform float Blur_Size;
uniform float Blur_Quality;
uniform vec2 resolution;
in vec2 uv;

layout (location = 0) out vec4 blur_depth_front;
layout (location = 1) out vec4 blur_depth_back;
layout (location = 2) out vec4 blur_thick;




vec4 blur(sampler2D image, vec2 uv, vec2 resolution) {

  float Pi = 6.28318530718; // Pi*2
    
  // GAUSSIAN BLUR SETTINGS {{{
  float Directions = 16.0; // BLUR DIRECTIONS (Default 16.0 - More is better but slower)
  float Quality = Blur_Quality; // BLUR QUALITY (Default 4.0 - More is better but slower)
  float Size = Blur_Size; // BLUR SIZE (Radius)
  // GAUSSIAN BLUR SETTINGS }}}

  vec2 Radius = Size/resolution.xy;
    
  // Pixel colour
  vec4 Color = texture(image, uv);
    
  float k = 0;
  for( float d=0.0; d<Pi; d+=Pi/Directions)
    {
		for(float i=1.0/Quality; i<=1.0; i+=1.0/Quality)
        {
          vec2 coords = uv+vec2(cos(d),sin(d))*Radius*i;
          if (coords.x >= 0 && coords.x <= 1 && coords.y >= 0 && coords.y <= 1){
            Color += texture( image, coords);		
            k+=1.0;
          }
        }
    }

  return Color/k;
}


vec4 blur5(sampler2D image, vec2 direction) {
  vec4 color = vec4(0.0);
  vec2 off1 = vec2(1.3333333333333333) * direction;
  color += texture(image, uv) * 0.29411764705882354;
  color += texture(image, uv + (off1 / resolution)) * 0.35294117647058826;
  color += texture(image, uv - (off1 / resolution)) * 0.35294117647058826;
  return color; 
}

vec4 blur13(sampler2D image, vec2 direction) {
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

vec4 smooth_Custom(sampler2D image){
  //if(Blur == 1)
    //  return (blur13(image,vec2(1,0)) + blur13(image,vec2(0,1)) + blur13(image,vec2(1,-1)) + blur13(image,vec2(1,1)))/4;
    return blur (image,uv,resolution);
  //return texture(image,uv);
}

#define EPS 1e-5

float linearize_depth(float d,float zNear,float zFar)
{
    return 0.5+0.5*zNear * zFar / (zFar + d * (zNear - zFar));
}

float lum(in vec4 color) {
    return linearize_depth(color.z,0.02,100.0);
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

uniform int filterRadius ;
uniform float sphereRadius ;
uniform float bilateralParameter ;

float linearize_depth(float d)
{
    float zNear = 0.02;
    float zFar = 100;
    return zNear * zFar / (zFar + d * (zNear - zFar));
}

float unlinearize_depth(float d)
{
    float zNear = 0.02;
    float zFar = 100;
    return (-zFar + (zNear * zFar / d) )/(zNear - zFar);
}

vec4 bilateral_new(sampler2D img){

	vec2 texCoord = uv;
	float depth = textureOffset(img, texCoord,ivec2(0)).x;

	float sum = (0.f);
	float wsum = (0.f);

    for(int y = -filterRadius; y <= filterRadius; y += 1) {
		for(int x = -filterRadius; x <= filterRadius; x += 1) {

			vec2 p = vec2(float(x), float(y))/resolution + texCoord;
			float value = textureOffset(img, texCoord,ivec2(x,y)).x;

			vec2 d = (p - texCoord) / filterRadius;
			float sqrDist = dot(d,d);
			float w = - sqrDist / (0.5f*0.5f);
      float r ;
      //if (depth > 0.1)
			r = (value - depth) / (sphereRadius*bilateralParameter);

      //else
			// r = (value - depth) / (0.01*sphereRadius * bilateralParameter);
			
			float g = -(r * r) / (0.5f);
			float weight = exp(w + g);
			sum += (value.x * weight.x);
			wsum += weight;

		}
	}

	if (wsum > 0.f) {
		sum /= wsum;
	}
	
	return vec4(vec3((sum)), 1.0);
}



uniform mat4 PROJECTION_MATRIX;
uniform mat4 PROJECTION_INVERSE_MATRIX;

// Mean curvature. From "Screen Space Fluid Rendering with Curvature Flow"
vec4 meanCurvature(sampler2D image) {
  vec2 pos = uv;
	// Width of one pixel
	vec2 dx = vec2(1.0f / resolution.x, 0.0f);
	vec2 dy = vec2(0.0f, 1.0f / resolution.y);

	// Central z value
	float zc =  texture(image, pos).x;

	// Take finite differences
	// Central differences give better results than one-sided here.
	// TODO better boundary conditions, possibly.
	// Remark: This is not easy, get to choose between bad oblique view smoothing
	// or merging of unrelated particles
	float zdxp = texture(image, pos + dx).x;
	float zdxn = texture(image, pos - dx).x;
	float zdx = 0.5f * (zdxp - zdxn);
	zdx = (zdxp == 0.0f || zdxn == 0.0f) ? 0.0f : zdx;

	float zdyp = texture(image, pos + dy).x;
	float zdyn = texture(image, pos - dy).x;
	float zdy = 0.5f * (zdyp - zdyn);
	zdy = (zdyp == 0.0f || zdyn == 0.0f) ? 0.0f : zdy;

	// Take second order finite differences
	float zdx2 = zdxp + zdxn - 2.0f * zc;
	float zdy2 = zdyp + zdyn - 2.0f * zc;

	// Second order finite differences, alternating variables
	float zdxpyp = texture(image, pos + dx + dy).z;
	float zdxnyn = texture(image, pos - dx - dy).z;
	float zdxpyn = texture(image, pos + dx - dy).z;
	float zdxnyp = texture(image, pos - dx + dy).z;
	float zdxy = (zdxpyp + zdxnyn - zdxpyn - zdxnyp) / 4.0f;

	// Projection transform inversion terms
	float cx = 2.0f / (resolution.x * -PROJECTION_MATRIX[0][0]);
	float cy = 2.0f / (resolution.y * -PROJECTION_MATRIX[1][1]);

	// Normalization term
	float d = cy * cy * zdx * zdx + cx * cx * zdy * zdy + cx * cx * cy * cy * zc * zc;
	
	// Derivatives of said term
	float ddx = cy * cy * 2.0f * zdx * zdx2 + cx * cx * 2.0f * zdy * zdxy + cx * cx * cy * cy * 2.0f * zc * zdx;
	float ddy = cy * cy * 2.0f * zdx * zdxy + cx * cx * 2.0f * zdy * zdy2 + cx * cx * cy * cy * 2.0f * zc * zdy;

	// Temporary variables to calculate mean curvature
	float ex = 0.5f * zdx * ddx - zdx2 * d;
	float ey = 0.5f * zdy * ddy - zdy2 * d;


  
	// Finally, mean curvature
	float h = 0.5f * ((cy * ex + cx * ey) / pow(d, 1.5f));


	const float dt = 0.00013f;
  const float dzt = 1000.0f;
  
  vec3 dxyz = vec3(zdx, zdy, h);

		// Vary contribution with absolute depth differential - trick from pySPH
  float outDepth = (zc) +( dxyz.z * dt * (1.0f + (abs(dxyz.x) + abs(dxyz.y)) * dzt));

	return vec4(vec3(2*outDepth),1);
}

void main()
{
    if(texture(frontdepthMap,uv).x==0){
        discard;
    }
    if(Blur > 0){
      if(Blur == 1){

        blur_depth_front = smooth_Custom(frontdepthMap);
        blur_depth_back = smooth_Custom(backdepthMap);
        blur_depth_back = texture(backdepthMap,uv);

      }
      else if (Blur == 2){ 
        blur_depth_front = bilateral_new(frontdepthMap);
        blur_depth_back = texture(backdepthMap,uv);
        blur_depth_back = bilateral_new(backdepthMap);
      }
      else if (Blur == 3) {
        blur_depth_front = vec4(vec3(meanCurvature(frontdepthMap).z),1);
        blur_depth_back =  vec4(vec3(meanCurvature(backdepthMap).z),1);
      }

        blur_thick = vec4(smooth_Custom(thickness).xyz,1); 
        //blur_thick = texture(thickness,uv);
      
    }
    else{
        blur_depth_front = texture(frontdepthMap,uv);
        blur_depth_back = texture(backdepthMap,uv);
        blur_thick = texture(thickness,uv);
    }

}