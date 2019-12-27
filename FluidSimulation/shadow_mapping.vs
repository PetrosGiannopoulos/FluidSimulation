#version 330 core
layout (location = 0) in vec3 aPos;
layout (location = 3) in vec3 aOffset;

uniform mat4 model;

uniform bool isParticle;

void main()
{
	vec4 pos = vec4(aPos,1.0);

	if(isParticle)pos = vec4(aPos+aOffset,1.0);

	gl_Position = model * pos;
	
}