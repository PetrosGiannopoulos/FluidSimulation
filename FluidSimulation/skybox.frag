#version 330 core
out vec4 FragColor;

in vec3 TexCoords;

uniform samplerCube skybox;

void main()
{    
	float gamma = 2.2;
    FragColor.xyz = pow(texture(skybox, TexCoords).xyz,vec3(gamma));
	FragColor.w = texture(skybox,TexCoords).w;

	 // apply gamma correction
    
    FragColor.rgb = pow(FragColor.rgb, vec3(1.0/gamma));

	//FragColor.rgb = texture(skybox, TexCoords).rgb;
}

