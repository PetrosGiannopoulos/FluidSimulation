#version 330 core

layout (location = 0) in vec3 aPos;
layout (location = 1) in vec3 aNormal;
layout (location = 2) in vec2 aTexCoords;
layout (location = 3) in vec3 aOffset;
layout (location = 4) in vec3 aTangent;
layout (location = 5) in vec3 aBitangent;

out vec2 TexCoords;

out VS_OUT {
    vec3 FragPos;
    vec3 Normal;
    vec2 TexCoords;

	vec3 TangentLightPos;
    vec3 TangentViewPos;
    vec3 TangentFragPos;

	vec3 reflectedVector;
	vec3 refractedVector;

	vec3 refractedVectorR;
	vec3 refractedVectorG;
	vec3 refractedVectorB;

} vs_out;

uniform mat4 projection;
uniform mat4 view;
uniform mat4 model;

uniform bool isParticle;

uniform bool reverse_normals;

uniform vec4 plane;

uniform vec3 lightPos;
uniform vec3 viewPosition;
uniform vec3 camPos;
uniform mat4 camView;

uniform float water;

uniform samplerCube depthMap;

uniform vec3 cubemapCenter;
uniform vec3 bboxmin;
uniform vec3 bboxmax;

out mat4 invView;
out mat4 invProj;
out mat4 invViewProj;
out mat4 ViewProj;

out vec3 viewRay;
out vec3 ViewNormal;

void main()
{

	vec4 pos = vec4(aPos,1.0);
	vec4 pPos = pos;
	

	if(isParticle)pos = vec4(aPos+aOffset,1.0);

	vec4 worldPosition = model*pos; 
	
    vs_out.FragPos = vec3(model * pos);

	//gl_ClipDistance[0] = dot(worldPosition,plane);

    if(reverse_normals) // a slight hack to make sure the outer large cube displays lighting from the 'inside' instead of the default 'outside'.
        vs_out.Normal = transpose(inverse(mat3(model))) * (-1.0 * aNormal);
    else
        vs_out.Normal = transpose(inverse(mat3(model))) * aNormal;
    vs_out.TexCoords = aTexCoords;

	mat3 normalMatrix = transpose(inverse(mat3(model)));
    vec3 T = normalize(normalMatrix * aTangent);
    vec3 N = normalize(normalMatrix * aNormal);
    T = normalize(T - dot(T, N) * N);
    vec3 B = cross(N, T);
    
    mat3 TBN = transpose(mat3(T, B, N));    
    vs_out.TangentLightPos = TBN * lightPos;
    vs_out.TangentViewPos  = TBN * viewPosition;
    vs_out.TangentFragPos  = TBN * vs_out.FragPos;

    gl_Position = projection * view * model * pos;

	//vs_out.FragPos = vec3(gl_Position);

	float waterRatio = 1./1.333;

	mat3 cameraRotation = mat3(view);

	mat4 invM = transpose(inverse(model));
	mat4 worldToLocal = transpose(inverse(model));
	mat3 invV = transpose(inverse(mat3(view)));
	invView = inverse(view);
	invViewProj = inverse(projection*view);
	ViewProj = view;
	mat3 M = mat3(model);
	mat3 V = mat3(view);

	invProj = (inverse(projection));

	vec3 vNormal = normalize((invM*vec4(aNormal,0.0)).xyz);

	//vNormal = normalize(model*vec4(aNormal,0.0)).xyz;
	mat4 viewNoTranslation = view; 
	viewNoTranslation[3] = vec4(0.0, 0.0, 0.0, 1.0);
	vNormal = normalize(viewNoTranslation*model*vec4(aNormal,0.0)).xyz;

	vec3 viewVector = normalize(worldPosition.xyz-camPos);

	viewVector = normalize((view*model*vec4(aPos-camPos,1.0)).xyz);

	viewRay = vec3(view*model*vec4(aPos,1.0));
	viewRay = vec3(viewRay.xy*(1000/viewRay.z),1000);

	vs_out.reflectedVector = normalize(reflect(viewVector,vNormal)); 

	//vs_out.reflectedVector = normalize(2.0f * vNormal * dot(viewVector, vNormal) - viewVector);

	vs_out.refractedVector = normalize(refract(viewVector,vNormal, waterRatio));

	//chromatic dispersion
	float waterRatioR = 1./1.33;
	float waterRatioG = 1./1.34;
	float waterRatioB = 1./1.35;

	vs_out.refractedVectorR = normalize(refract(viewVector,vNormal, waterRatioR));
	vs_out.refractedVectorG = normalize(refract(viewVector,vNormal, waterRatioG));
	vs_out.refractedVectorB = normalize(refract(viewVector,vNormal, waterRatioB));

	//vec3 RayLS = vec3(worldToLocal*vec4(vs_out.reflectedVector,0.0));
	//vec3 PositionLS = vec3(worldToLocal*vec4(worldPosition.xyz,1.0));

	//float boxSize = max(bboxmin.x,bboxmax.x);

	//vec3 Unitary = vec3(1.0);
	//vec3 FirstPlaneIntersect = (Unitary-PositionLS)/RayLS;
	//vec3 SecondPlaneIntersect = (-Unitary-PositionLS)/RayLS;
	//vec3 FurthestPlane = max(FirstPlaneIntersect, SecondPlaneIntersect);
	//float Distance = min(FurthestPlane.x, min(FurthestPlane.y, FurthestPlane.z));

	//vec3 IntersectPosition = worldPosition.xyz+vs_out.reflectedVector*Distance;
	//vs_out.reflectedVector = IntersectPosition - cubemapCenter;

	//float radius = 20.0f;
	//float EnvMapOffset = 1./radius;
	//vs_out.reflectedVector = EnvMapOffset *(worldPosition.xyz-cubemapCenter)+vs_out.reflectedVector;
}