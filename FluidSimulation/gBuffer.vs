#version 330 core
layout (location = 0) in vec3 aPos;
layout (location = 1) in vec3 aNormal;
layout (location = 2) in vec2 aTexCoords;

out vec3 FragPos;
out vec2 TexCoords;
out vec3 Normal;
out mat4 invProjection;
out mat4 viewMat;
out vec4 Position;

uniform mat4 model;
uniform mat4 view;
uniform mat4 projection;

uniform bool isParticle;
uniform bool reverse_normals;
uniform vec3 camPos;

void main()
{
    vec4 worldPos = model * vec4(aPos, 1.0);
	vec4 viewPos = view*model * vec4(aPos, 1.0);
    FragPos = viewPos.xyz; 
    TexCoords = aTexCoords;

	viewMat = view;
	//Normal = vec3(view*model*vec4(aNormal,0.0));

	mat3 normalMatrix = transpose(inverse(mat3(view * model)));
    Normal = normalMatrix * aNormal;

    gl_Position = projection * view * worldPos;

	Position = gl_Position;
	invProjection = inverse(projection);
}