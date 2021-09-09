#version 330

uniform samplerCube CMHouse;


uniform vec2  view_size;

in data {
    float light_height;
    vec4 position;
} o;

out vec4 color;

void main()
{
    //float dist = length(gl_FragCoord.xy - sun_screenPos * view_size);
    
    /* Simple circle */
    //float sun_value = dist < sun_cutoff ? 1.0 : 0.0;
    
    /* Fade based on distance */
    //float sun_value = (1 / dist) * sun_cutoff;
    
    /* Gaussian Distribution */
    //float sun_value = pow(2.17, -(dist*dist) / (2 * pow(sun_cutoff, 2))) * sun_fade;
    //vec3 sky_color = mix(sky_orange, sky_blue, o.light_height);

    color = texture(CMHouse, (-1)*o.position.xyz);
}