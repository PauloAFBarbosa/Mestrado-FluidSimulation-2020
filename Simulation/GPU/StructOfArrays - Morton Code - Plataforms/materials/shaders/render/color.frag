#version 430

uniform sampler2D backgroundTex;
uniform sampler2D fluidColor;

uniform sampler2D front_depth;
uniform sampler2D back_depth;
uniform sampler2D depthMap;

uniform samplerCube texUnit;

uniform vec3 lightDir;
uniform vec4 Particle_Color;
uniform vec4 Thin_Color;
uniform float eta;
uniform int DebugParticles;
uniform vec4 m_campos;
uniform mat4 m_inv_proj;
uniform mat4 m_inv_view;
uniform vec2 viewport_resolution;
uniform vec2 fluid_resolution;
uniform int thicktype;

in vec2 uv;

out vec4 color;
out float gl_FragDepth;

uniform float Particle_Thickness;

vec2 resolution = vec2(viewport_resolution);

vec4 getWorldPosBack(vec2 offset)
{
	vec2 coords = uv + (offset / resolution);
	vec4 tex = (texture(front_depth, coords));
	vec4 clip_space_pos;

	clip_space_pos = vec4(coords * 2.0 - 1.0, (tex.g), 1.0);

	vec4 viewspace = m_inv_proj * clip_space_pos;
	vec4 worldspacepos = m_inv_view * ((viewspace / viewspace.w));

	return worldspacepos;
}

vec4 getWorldPosFront(vec2 offset)
{
	vec2 coords = uv + (offset / resolution);
	vec4 tex = (texture(front_depth, coords));
	vec4 clip_space_pos;

	clip_space_pos = vec4(coords * 2.0 - 1.0, (tex.r), 1.0);

	//vec4 clip_space_pos = vec4(coords * 2.0 - 1.0, texture(front_depth,coords).x, 1.0);
	vec4 viewspace = m_inv_proj * clip_space_pos;
	vec4 worldspacepos = m_inv_view * (viewspace / viewspace.w);

	return worldspacepos;
}

vec4 normalFromFrontDepth()
{

	vec4 t = getWorldPosFront(vec2(1, 0));
	vec4 b = getWorldPosFront(vec2(-1, 0));
	vec4 r = getWorldPosFront(vec2(0, 1));
	vec4 l = getWorldPosFront(vec2(0, -1));
	return vec4(normalize(cross(normalize((t - b).xyz), normalize((r - l).xyz))), 0);

	//return -1*vec4(normalize(cross((t-c).xyz,(l-c).xyz)),0);
}

//vec4 normalFromBackDepth(){
//  vec4 t = getWorldPosBack(vec2( 1,0));
//  vec4 b = getWorldPosBack(vec2(-1,0));
//  vec4 r = getWorldPosBack(vec2(0, 1));
//  vec4 l = getWorldPosBack(vec2(0,-1));
//
//  return  vec4(normalize(cross(normalize((t-b).xyz),normalize((r-l).xyz))),0);
//}
vec4 worldspace(vec3 screenSpace)
{

	vec4 clip_space_pos;

	clip_space_pos = vec4(screenSpace.xy, screenSpace.z, 1);
	;

	vec4 viewspace = m_inv_proj * clip_space_pos;
	vec4 worldspacepos = m_inv_view * ((viewspace / viewspace.w));

	return worldspacepos;
}

vec4 getWorldPosBack(vec2 pos, vec2 offset)
{
	vec2 coords = pos + (offset / resolution);
	vec4 tex = (texture(front_depth, coords));
	vec4 clip_space_pos;

	clip_space_pos = vec4(coords * 2.0 - 1.0, (tex.g), 1.0);

	vec4 viewspace = m_inv_proj * clip_space_pos;
	vec4 worldspacepos = m_inv_view * ((viewspace / viewspace.w));

	return worldspacepos;
}

vec4 normalFromBackDepth(vec2 pos)
{
	vec4 t = getWorldPosBack(pos, vec2(1, 0));
	vec4 b = getWorldPosBack(pos, vec2(-1, 0));
	vec4 r = getWorldPosBack(pos, vec2(0, 1));
	vec4 l = getWorldPosBack(pos, vec2(0, -1));

	return vec4(normalize(cross(normalize((t - b).xyz), normalize((r - l).xyz))), 0);
}

uniform int maxSteps;
uniform float stepSize;
uniform int power;
uniform int specular_power;

int z = 0;
vec3 val;
vec3 raycast()
{

	float j = 0;
	vec4 t = vec4(uv + vec2(1, 0) / fluid_resolution, texture(front_depth, uv + vec2(1, 0) / fluid_resolution).x, 1);
	vec4 b = vec4(uv + vec2(-1, 0) / fluid_resolution, texture(front_depth, uv + vec2(-1, 0) / fluid_resolution).x, 1);
	vec4 r = vec4(uv + vec2(0, 1) / fluid_resolution, texture(front_depth, uv + vec2(0, 1) / fluid_resolution).x, 1);
	vec4 l = vec4(uv + vec2(0, -1) / fluid_resolution, texture(front_depth, uv + vec2(0, -1) / fluid_resolution).x, 1);
	vec3 normal = normalize(cross(normalize((t - b).xyz), normalize((r - l).xyz)));
	vec3 raydir = -1 * normalize(refract(vec3(0, 0, 1), normal, eta));

	vec3 rayPos = vec3(uv, texture(front_depth, uv).x);
	vec3 oldPos = vec3(uv, texture(front_depth, uv).x);

	val = texture(front_depth, rayPos.xy).xyz;
	float stepLength = stepSize ;//* (val.y - val.x); // stepSize;// max(,stepSize);
	if (stepLength == 0)
		return rayPos;
	while (true)
	{

		val = texture(front_depth, rayPos.xy).xyz;

		vec3 rayStep = (stepLength)*raydir;
		if ((rayPos.x > 1.0 || rayPos.x < 0.0 || rayPos.y > 1.0 || rayPos.y < 0.0))
			return oldPos;
		//if(rayPos.z > val.x)
		//	return oldPos;
		if (val.z == 0)
		{
			rayPos = oldPos + (rayPos - oldPos) / 2;
			return rayPos;
		}
		else
		{
			if (rayPos.z > val.y)
			{
				return oldPos + (rayPos - oldPos) / 2;
			}
			oldPos = rayPos;
			rayPos += raydir * (stepLength);
		}
		z++;
		if (z > maxSteps)
		{
			return rayPos;
		}
	}
	return rayPos;
}

vec3 getSkyboxColor(vec3 dir)
{
	return texture(texUnit, (-1 * normalize(dir))).rgb;
	//return (dir * 0.5 + vec3(1));
}

float etaInverse = 1 / eta;

vec3 refratBackColor(vec3 inciding, vec3 ldir, out vec3 r2)
{
	vec3 raypos = raycast();
	r2 = raypos;
	vec3 n = normalize(normalFromBackDepth(raypos.xy).xyz);
	vec3 e = normalize(inciding);
	vec3 t = normalize(reflect(e, n));
	vec3 ref = getSkyboxColor(t);
	vec3 r = normalize(refract(e, n, etaInverse));
	//float k = 1.0 - eta * eta * (1.0 - dot(n, e) * dot(n, e));
	vec3 refrat;
	refrat = getSkyboxColor(r);

	float f = pow(((1 - etaInverse) * (1 - etaInverse)) / ((1 + etaInverse) * (1 + etaInverse)),2);
	float shlick = f + (1 - f) * pow(1 - dot(-e, n), power);
	vec3 refRefac = mix(refrat, ref, shlick);
	vec3 h = normalize(ldir - e);
	// compute the specular intensity
	float intSpec = max(dot(h, n), 0.0);
	// compute the specular term into spec
	vec4 spec = vec4(1) * pow(intSpec, specular_power);
	vec3 outvar = (refRefac);
	// if(outvar.x<0 || outvar.y<0 || outvar.z <0)

	return max(vec3(0), outvar + spec.xyz);
}

uniform int Surface_Num;

void main()
{

	// color = texture(spherePosTex,uv);  return;
	float depth = texture(front_depth, uv).z;
	if (DebugParticles == 0)
		color = texture(backgroundTex, uv);
	//vec2 resolution = vec2(1920,1080);
	vec4 thickness;
	//thickness.xyz = iterate();

	vec3 normal = normalFromFrontDepth().xyz; //normalize((smooth_Custom(normalMapSphere,uv,resolution)));

	//if(depth == 0)
	//  normal = normalFromBackDepth().xyz;

	vec4 worldpos = getWorldPosFront(vec2(0));

	vec4 lookAtVec = normalize(worldpos - m_campos);

	vec4 particleColor = texture(fluidColor, uv);
	vec3 n = normalize(normal.xyz);
	vec3 ldir = normalize(-lightDir);

	float diff = (0.2 + 0.8 * max(dot(n, ldir), 0.0));
	if (DebugParticles != 0)
	{
		if (depth == 1)
		{
			if (DebugParticles == 1){
				color = vec4(vec3(diff) * Particle_Color.xyz, 1);
				return;
			}
			else if (DebugParticles == 2){
				color = vec4(vec3((worldpos)), 1);
				return;
			}
			else if (DebugParticles == 3){
				color = vec4(vec3((normal + vec3(1)) / 2), 1);
				return;
			}
			else if (DebugParticles == 5)
			{
				color = vec4(vec3(worldspace(raycast().xyz)), 1);
				return;
			}
		}
	}

	if (depth == 1)
	{
		vec3 n = normalize(normal.xyz);
		vec3 e = normalize(lookAtVec.xyz);
		vec3 t = normalize(reflect(e, n));
		vec3 ref = getSkyboxColor(t);

		//vec3 refratnormal = normalize(normalFromBackDepth().xyz);

		vec3 ldir = normalize(-lightDir);

		vec3 refrat;
		vec3 r = normalize(refract(e, n, eta));

		thickness = texture(fluidColor, uv);
		float concentration = Particle_Thickness;
		if (Surface_Num == 1)
		{
			refrat = getSkyboxColor(r);
		}
		else if (Surface_Num == 2)
		{
			vec3 r2;
			refrat = refratBackColor(r.xyz, ldir,r2); // vec3(texture(texUnit, (-1*r)));
			vec3 backworld = worldspace(vec3(r2.xy,texture(front_depth,r2.xy).y)).xyz;
			vec3 frontworld = worldspace(vec3(uv,texture(front_depth,r2.xy).x)).xyz;
			if (thicktype == 1)
			{
				thickness = vec4(distance(frontworld, backworld));
				//concentration /= 4;
			}
			else if (thicktype == 3)
			{
				thickness = vec4(distance(getWorldPosFront(vec2(0)).xyz, getWorldPosBack(vec2(0)).xyz));
			//
				//thickness = vec4(
				//		min(vec4(distance(getWorldPosFront(vec2(0)).xyz, backworld)),
				//			texture(fluidColor, uv).x));
				//concentration /= 4;
			}
		}
		
		thickness *= concentration;
		if (DebugParticles == 4){
			color = vec4((thickness.xyz), 1);
			return;
		}

		vec3 h = normalize(ldir - e);
		// compute the specular intensity
		float intSpec = max(dot(h, n), 0.0);
		// compute the specular term into spec
		vec4 spec = vec4(1) * pow(intSpec, specular_power);

		refrat = max(vec3(0), refrat - (thickness.x * (1 - Particle_Color.xyz))); //mix(refrat,refrat*(Particle_Color.xyz),smoothstep(0,1,(thickness.x)));
		float f = pow(((1 - eta) * (1 - eta)) / ((1 + eta) * (1 + eta)),2);
		float shlick = f + (1 - f) * pow(1 - dot(-e, n), power);
		vec3 refRefac = mix(refrat, ref, shlick);
		//if(refrat.x==0)
		//	refRefac = ref;
		refRefac = mix(refRefac, Thin_Color.xyz * refRefac, smoothstep(0, 1, (Thin_Color.a)));

		color = vec4(
			(refRefac + spec.xyz), 1.0);
	}
}