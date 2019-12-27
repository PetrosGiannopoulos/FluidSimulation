#version 330 core
layout (location = 0) out vec3 gPosition;
layout (location = 1) out vec3 gNormal;
layout (location = 2) out vec4 gAlbedoSpec;
layout (location = 3) out vec3 gDepth;
layout (location = 4) out vec3 gColor;

in vec2 TexCoords;
in vec3 FragPos;
in vec3 Normal;
in mat4 invProjection;
in mat4 viewMat;
in vec4 Position;

uniform sampler2D texture_diffuse1;
uniform sampler2D texture_specular1;
uniform sampler2D textureFrame;

float near = 0.1; 
float far  = 1000.0; 

uniform vec3 camPos;

float LinearizeDepth(float depth) 
{
    float z = depth * 2.0 - 1.0; // back to NDC 
    return (2.0 * near * far) / (far + near - z * (far - near));	
}

vec3 LinearizeDepthV(vec3 depth){

	depth.x = LinearizeDepth(depth.x);
	depth.y = LinearizeDepth(depth.y);
	depth.z = LinearizeDepth(depth.z);

	return depth;
}

vec3 PositionFromDepth(float depth) {
    float z = depth*2.0-1.0;

	//z = LinearizeDepth(z)/far;

    vec4 clipSpacePosition = vec4((gl_FragCoord.xy/vec2(1000,800))*2.0-1.0, z, 1.0);
    vec4 viewSpacePosition = invProjection * clipSpacePosition;

    // Perspective division
    viewSpacePosition /= viewSpacePosition.w;

    return viewSpacePosition.xyz;
}


void main()
{    
    // store the fragment position vector in the first gbuffer texture
	
    //gPosition = vec3(FragPos.z);
	//gPosition.z = LinearizeDepth(FragPos.z)/far;
	gPosition = FragPos;
    // also store the per-fragment normals into the gbuffer
    gNormal = normalize(Normal);
    // and the diffuse per-fragment color
    gAlbedoSpec.rgb = vec3(0.95);//texture(texture_diffuse1, TexCoords).rgb;
    // store specular intensity in gAlbedoSpec's alpha component
    gAlbedoSpec.a = texture(texture_specular1, TexCoords).r;

	gColor.rgb = texture(textureFrame, TexCoords).rgb;

	float viewDepth = FragPos.z;//(viewMat*gl_FragCoord).z;

	float depth = gl_FragCoord.z;//viewDepth/far;//Position.z;
	//depth = LinearizeDepth(depth)/far;
	gDepth = vec3(depth);
}