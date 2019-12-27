#version 330 core
layout(location = 0) in vec3 aPos;

out vec3 TexCoords;

uniform mat4 projection;
uniform mat4 view;
uniform vec3 rotation;

const float PI = 3.14159265359;

void main()
{

	TexCoords = aPos;
	vec3 pPos = aPos;

	vec4 pos = projection * view * vec4(pPos, 1.0);

	

    gl_Position = pos.xyww;
}
