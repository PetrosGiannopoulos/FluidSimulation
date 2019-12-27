#version 330 core

out vec4 FragColor;

#extension GL_NV_shadow_samplers_cube : enable
in VS_OUT {
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
} fs_in;

uniform sampler2D diffuseTexture;
uniform samplerCube depthMap;
uniform samplerCube skybox;
uniform sampler2D normalMap;

uniform sampler2D gPosition;
uniform sampler2D gNormal;
uniform sampler2D gAlbedoSpec;
uniform sampler2D gPrevFrame;
uniform sampler2D gDepth;
uniform sampler2D ssao;

uniform vec3 lightPos;
uniform vec3 viewPosition;
uniform vec3 camPos;
uniform vec4 colorVar;

uniform float water;

uniform float far_plane;
uniform bool shadows;

uniform vec3 sphereCenter;
uniform bool isSSR;
uniform mat4 cameraToClipMatrix;

uniform float angleInput;
uniform mat4 projection, view;

in mat4 invView;
in mat4 invProj;
in vec3 ViewNormal;
in vec3 viewRay;
in mat4 invViewProj;
in mat4 ViewProj;


#define Scale vec3(.8, .8, .8)
#define K 19.19

const float step = 0.1;
const float minRayStep = 0.1;
const int maxSteps = 30;
const int numBinarySearchSteps = 5;
const float LLimiter = 0.2;
const float reflectionSpecularFalloffExponent = 3.0;

const float PI = 3.14159265359;
// array of offset direction for sampling
vec3 gridSamplingDisk[20] = vec3[]
(
   vec3(1, 1,  1), vec3( 1, -1,  1), vec3(-1, -1,  1), vec3(-1, 1,  1), 
   vec3(1, 1, -1), vec3( 1, -1, -1), vec3(-1, -1, -1), vec3(-1, 1, -1),
   vec3(1, 1,  0), vec3( 1, -1,  0), vec3(-1, -1,  0), vec3(-1, 1,  0),
   vec3(1, 0,  1), vec3(-1,  0,  1), vec3( 1,  0, -1), vec3(-1, 0, -1),
   vec3(0, 1,  1), vec3( 0, -1,  1), vec3( 0, -1, -1), vec3( 0, 1, -1)
);

float rand(vec2 n) { 
	return fract(sin(dot(n, vec2(12.9898, 4.1414))) * 43758.5453);
}

float noise(vec2 n) {
	const vec2 d = vec2(0.0, 1.0);
	vec2 b = floor(n), f = smoothstep(vec2(0.0), vec2(1.0), fract(n));
	return mix(mix(rand(b), rand(b + d.yx), f.x), mix(rand(b + d.xy), rand(b + d.yy), f.x), f.y);
}

float fbm(vec2 n) {
	float total = 0.0, amplitude = 1.0;
	for (int i = 0; i < 7; i++) {
		total += noise(n) * amplitude;
		n += n;
		amplitude *= 0.5;
	}
	return total;
}

vec3 hash(vec3 a)
{
    a = fract(a * Scale);
    a += dot(a, a.yxz + K);
    return fract((a.xxy + a.yxx)*a.zyx);
}

vec3 hash33(vec3 p){     
    float n = sin(dot(p, vec3(7, 157, 113)));    
    return fract(vec3(2097152, 262144, 32768)*n); 
}

float ShadowCalculation(vec3 fragPos)
{
    // get vector between fragment position and light position
    vec3 fragToLight = fragPos - lightPos;
    // use the fragment to light vector to sample from the depth map    
    // float closestDepth = texture(depthMap, fragToLight).r;
    // it is currently in linear range between [0,1], let's re-transform it back to original depth value
    // closestDepth *= far_plane;
    // now get current linear depth as the length between the fragment and light position
    float currentDepth = length(fragToLight);
    // test for shadows
    // float bias = 0.05; // we use a much larger bias since depth is now in [near_plane, far_plane] range
    // float shadow = currentDepth -  bias > closestDepth ? 1.0 : 0.0;
    // PCF
    // float shadow = 0.0;
    // float bias = 0.05; 
    // float samples = 4.0;
    // float offset = 0.1;
    // for(float x = -offset; x < offset; x += offset / (samples * 0.5))
    // {
        // for(float y = -offset; y < offset; y += offset / (samples * 0.5))
        // {
            // for(float z = -offset; z < offset; z += offset / (samples * 0.5))
            // {
                // float closestDepth = texture(depthMap, fragToLight + vec3(x, y, z)).r; // use lightdir to lookup cubemap
                // closestDepth *= far_plane;   // Undo mapping [0;1]
                // if(currentDepth - bias > closestDepth)
                    // shadow += 1.0;
            // }
        // }
    // }
    // shadow /= (samples * samples * samples);
    float shadow = 0.0;
    float bias = 0.15;
    int samples = 20;
    float viewDistance = length(viewPosition - fragPos);
    float diskRadius = (1.0 + (viewDistance / far_plane)) / 25.0;
    for(int i = 0; i < samples; ++i)
    {
        float closestDepth = texture(depthMap, fragToLight + gridSamplingDisk[i] * diskRadius).r;
        closestDepth *= far_plane;   // undo mapping [0;1]
        if(currentDepth - bias > closestDepth)
            shadow += 1.0;
    }
    shadow /= float(samples);
        
    // display closestDepth as debug (to visualize depth cubemap)
    //FragColor = vec4(vec3(closestDepth / far_plane), 1.0);    
        
    return shadow;
}

vec3 BlinnPhong(vec3 normal, vec3 fragPos, vec3 lightPos, vec3 lightColor)
{
    // diffuse
    vec3 lightDir = normalize(lightPos - fragPos);
    float diff = max(dot(lightDir, normal), 0.0);
    vec3 diffuse = diff * lightColor;
    // specular
    vec3 viewDir = normalize(camPos - fragPos);
    vec3 reflectDir = reflect(-lightDir, normal);
    float spec = 0.0;
    vec3 halfwayDir = normalize(lightDir + viewDir);  
    //spec = pow(max(dot(normal, reflectDir), 0.0), 64.0);
	spec = pow(max(dot(normal, halfwayDir), 0.0), 64.0);
    vec3 specular = spec * lightColor;    
    // simple attenuation
    float max_distance = 1.5;
    float distance = length(lightPos - fragPos);
    float attenuation = 1.0 / (distance*distance);
    
	attenuation *= 40.2;
    diffuse *= attenuation;
    specular *= attenuation;
    
    return diffuse + specular;
}

vec3 fresnelSchlick(float cosTheta, vec3 F0)
{
    //return F0 + (1.0 - F0) * pow(1.0 - cosTheta, 5.0);
	return max(F0 + (1.0 - F0) * pow(1.0 - cosTheta, 5.0), 0.0);
}  

vec3 fresnelSchlick2(float cosTheta, vec3 F0)
{
    return F0 + (1.0 - F0) * pow(1.0 - cosTheta, 5.0);
}

float DistributionGGX(vec3 N, vec3 H, float roughness)
{
    float a      = roughness*roughness;
    float a2     = a*a;
    float NdotH  = max(dot(N, H), 0.0);
    float NdotH2 = NdotH*NdotH;
	
    float num   = a2;
    float denom = (NdotH2 * (a2 - 1.0) + 1.0);
    denom = PI * denom * denom;
	
    return num / denom;
}

float GeometrySchlickGGX(float NdotV, float roughness)
{
    float r = (roughness + 1.0);
    float k = (r*r) / 8.0;

    float num   = NdotV;
    float denom = NdotV * (1.0 - k) + k;
	
    return num / denom;
}
float GeometrySmith(vec3 N, vec3 V, vec3 L, float roughness)
{
    float NdotV = max(dot(N, V), 0.0);
    float NdotL = max(dot(N, L), 0.0);
    float ggx2  = GeometrySchlickGGX(NdotV, roughness);
    float ggx1  = GeometrySchlickGGX(NdotL, roughness);
	
    return ggx1 * ggx2;
}

vec2 uvSphere(vec3 p){

	float u = 0.5 + atan(p.z,p.x)/(2*PI);
	float v = 0.5 - asin(p.y)/PI;

	return vec2(u,v);
}

float near = 0.1; 
float far  = 1000.0; 
  
float LinearizeDepth(float depth) 
{
    float z = depth * 2.0 - 1.0; // back to NDC 
    //return (2.0 * near * far) / (far + near - z * (far - near));	

	float linDepth = (2.0 * near) / (far + near - depth * (far - near));
	return linDepth;
	
}

vec2 getScreenSpacePosition(){
	return gl_FragCoord.xy/vec2(1000,800);
}

float ComputeDepth(vec4 clippos){
	return (clippos.z / clippos.w) * 0.5 + 0.5;
}

vec3 CalcViewPosition(in vec2 uv)
{
    // Combine UV & depth into XY & Z (NDC)
    vec3 rawPosition                = vec3(uv, texture(gDepth, uv).z);

    // Convert from (0, 1) range to (-1, 1)
    vec4 ScreenSpacePosition        = vec4( rawPosition*2-1, 1);

    // Undo Perspective transformation to bring into view space
    vec4 ViewPosition               = invProj * ScreenSpacePosition ;

    // Perform perspective divide and return
    return                          ViewPosition.xyz / ViewPosition.w;
}

//Convert something in camera space to screen space
vec3 convertCameraSpaceToScreenSpace(in vec3 cameraSpace)
{
	vec4 clipSpace = cameraToClipMatrix * vec4(cameraSpace, 1);
	vec3 NDCSpace = clipSpace.xyz / clipSpace.w;
	vec3 screenSpace = 0.5 * NDCSpace + 0.5;
	return screenSpace;
}

vec3 PositionFromDepth(float depth) {
    float z = depth * 2.0 - 1.0;

    vec4 clipSpacePosition = vec4(getScreenSpacePosition() * 2.0 - 1.0, z, 1.0);
    vec4 viewSpacePosition = inverse(projection) * clipSpacePosition;

    // Perspective division
    viewSpacePosition /= viewSpacePosition.w;

    return viewSpacePosition.xyz;
}

vec2 BinarySearch(inout vec3 dir, inout vec3 hitCoord, inout float dDepth)
{
    float depth;

    vec4 projectedCoord;
 
    for(int i = 0; i < numBinarySearchSteps; i++){
        projectedCoord = projection * vec4(hitCoord, 1.0);
        projectedCoord.xy /= projectedCoord.w;
        projectedCoord.xy = projectedCoord.xy * 0.5 + 0.5;

        depth =  textureLod(gPosition,projectedCoord.xy,0).z;
		dDepth = (hitCoord.z - depth);

        dir *= 0.5;
        if(dDepth > 0.0)
            hitCoord += dir;
        else
            hitCoord -= dir;    
    }

    projectedCoord = projection * vec4(hitCoord, 1.0);
    projectedCoord.xy /= projectedCoord.w;
    projectedCoord.xy = projectedCoord.xy * 0.5 + 0.5;
 
    return vec2(projectedCoord.xy);
}

float GetAngle(vec3 v1, vec3 v2) {
		
	//v1.v2 = |v1|*|v2|*cos(theta);
	return acos(dot(v1, v2)/(length(v1)*length(v2)))*(180/PI);
}

float random(vec2 co) {
    return fract(sin(dot(co.xy, vec2(12.9898, 78.233))) * 43758.5453123);
}

vec2 RayMarch(vec3 dir, inout vec3 hitCoord,out float dDepth, out bool success)
{

    dir *= step;
  
	vec4 projectedCoord;

	const float df = 1.2;
	const float angleDf = 168;

	vec3 P = hitCoord;

    for(int i = 0; i < maxSteps; i++){
        hitCoord += dir;
 
        projectedCoord = projection * vec4(hitCoord, 1.0);
        projectedCoord.xy /= projectedCoord.w;
        projectedCoord.xy = projectedCoord.xy * 0.5 + 0.5;

        float depth = textureLod(gPosition,projectedCoord.xy,0).z;
		dDepth = (hitCoord.z - depth);
		
		if(depth>1000)continue;

		vec3 B = textureLod(gPosition,projectedCoord.xy,0).xyz;
		vec3 A = hitCoord;

		vec3 PA = P-A;
		vec3 BP = B-P;

		float angleV = GetAngle(PA,BP);

		//if(angleV>angleDf){
        if((dir.z-dDepth) < df)
		if(dDepth <= 0.0) {
			success = true;
			return BinarySearch(dir, hitCoord, dDepth);
		}
        //if(dDepth < 0.0) return projectedCoord.xy;
    }
	
    return projectedCoord.xy;
}

vec3 GetSpecular(vec3 N, vec3 V, vec3 L, vec3 H, vec3 F){

	float roughness = 0.7;
	float NDF = DistributionGGX(N, H, roughness);       
	float G   = GeometrySmith(N, V, L, roughness); 

	vec3 numerator    = NDF * G * F;
	float denominator = 4.0 * max(dot(N, V), 0.0) * max(dot(N, L), 0.0);
	vec3 specular     = numerator / max(denominator, 0.001); 

	return specular;
}

void main()
{           
	vec2 sphereCoords = uvSphere(normalize(fs_in.FragPos-sphereCenter));

    vec3 color = texture(diffuseTexture, fs_in.TexCoords).rgb;
						
	vec3 sphereColor = texture(diffuseTexture,sphereCoords).rgb;

    vec3 normal = normalize(fs_in.Normal);
    vec3 lightColor = vec3(0.7);
	//vec3 lightColor  = vec3(23.47, 21.31, 20.79);
	
    // ambient
	float AmbientOcclusion;
	vec3 waterColor = colorVar.xyz;
	//waterColor.xyz = vec3(0,0,0);
    vec3 ambient = 0.7 * (color)*AmbientOcclusion;//+waterColor);
	if(isSSR){
		AmbientOcclusion = texture(ssao, getScreenSpacePosition()).r;
		ambient *= AmbientOcclusion;
	}
   
	if(water==1){
		if(isSSR)ambient = 0.7*lightColor*AmbientOcclusion;
		else ambient = 0.7*lightColor;
		color = lightColor;
	}
    // diffuse
    vec3 lightDir = normalize(lightPos - fs_in.FragPos);
	//float distance = length(lightPos-fs_in.FragPos);
	//float attenuation = 1.0 / (distance*distance); 
	//color *= attenuation;
    float diff = max(dot(lightDir, normal), 0.0);
    //vec3 diffuse = diff * lightColor;
    // specular
    vec3 viewDir = normalize(viewPosition - fs_in.FragPos);
    vec3 reflectDir = reflect(-lightDir, normal);
    float spec = 0.0;
    vec3 halfwayDir = normalize(lightDir + viewDir);  
    spec = pow(max(dot(normal, halfwayDir), 0.0), 64.0);
    //vec3 specular = spec * lightColor;    

	//normal mapping
	vec3 normalPlane = texture(normalMap, fs_in.TexCoords).rgb;
    // transform normal vector to range [-1,1]
    normalPlane = normalize(normalPlane * 2.0 - 1.0);

	vec3 specularP;
	vec3 diffuseP;
	if(water != 1.0 && water != 2.0){
		normal = normalPlane;
		lightDir = normalize(fs_in.TangentLightPos - fs_in.TangentFragPos);
		diff = max(dot(lightDir, normal), 0.0);
		diffuseP = diff * color;
		viewDir = normalize(fs_in.TangentViewPos - fs_in.TangentFragPos);
		reflectDir = reflect(-lightDir, normal);
        halfwayDir = normalize(lightDir + viewDir);  
        spec = pow(max(dot(normal, halfwayDir), 0.0), 32.0);

		specularP = vec3(0.2) * spec;
		
	}

	vec3 blinnPhong = BlinnPhong(normal,fs_in.FragPos,lightPos,lightColor);
	vec3 blinnPhongSphere = BlinnPhong(normal,fs_in.FragPos,lightPos,lightColor);
    // calculate shadow
    float shadow = shadows ? ShadowCalculation(fs_in.FragPos) : 0.0;                      
    
	vec3 lighting = (ambient+ (1.0 - shadow) * (blinnPhong))*color;
	vec3 lightingSphere = (ambient + (1.0 - shadow) * (blinnPhongSphere))*sphereColor;   
	//lighting = pow(lighting,vec3(1/2.2));
	//vec3 lighting = (ambient + (1.0 - shadow) * (diffuse + specular))*color;    
	vec3 N = normalize(fs_in.Normal);
	if(water != 1.0 && water != 2.0)N = normalPlane;
	vec3 V = normalize(camPos - fs_in.FragPos);
	vec3 L = normalize(lightPos - fs_in.FragPos);
	vec3 H = normalize(V + L);

	float distance    = length(lightPos - fs_in.FragPos);
	float attenuation = 1.0 / (distance*distance );
	//attenuation *= 40;
	vec3 radiance     = lightColor * attenuation;

	vec4 reflectedColour, refractedColour;
	if(isSSR){
		//SSR Algorithm
		vec2 screenSpacePosition2D = getScreenSpacePosition();
		vec2 uv = screenSpacePosition2D;

		vec3 viewNormal = vec3(texture(gNormal, uv));
		//vec3 viewNormal = vec3(texture2D(gNormal, fs_in.TexCoords));
		vec3 viewPos_ = texture(gPosition,uv).xyz;
		vec3 albedo_ = texture(gPrevFrame, uv).rgb;
		float spec_ = 0.7;//texture(gAlbedoSpec, uv).a;

		vec3 F0_ = vec3(0.04); 
		F0_      = mix(F0_, albedo_, 1.0);
		vec3 Fresnel = fresnelSchlick2(max(dot(normalize(viewNormal), normalize(viewPos_)), 0.0), F0_);
	
		vec3 reflected = normalize(reflect((viewPos_), normalize(viewNormal)));
		//reflected = cameraSpaceVector;

		vec3 hitPos = viewPos_;
		float dDepth;
		bool success = false;

		vec3 floorSpecular = GetSpecular(N,V,L,H, Fresnel);
	
 
		vec3 wp = vec3(vec4(viewPos_, 1.0)*invView);
		vec3 jitt = vec3(hash(viewPos_*floorSpecular));//mix(vec3(0.0), vec3(hash(wp)), spec_);
		vec2 coords = RayMarch((jitt+reflected * max(-viewPos_.z,minRayStep)), hitPos, dDepth,success);
		vec2 dCoords = smoothstep(0.2, 0.6, abs(vec2(0.5, 0.5) - coords));
		float screenEdgefactor = clamp(1.0 - (dCoords.x + dCoords.y), 0.0, 1.0);
		float ReflectionMultiplier = pow(1.0, reflectionSpecularFalloffExponent) * screenEdgefactor * -reflected.z;
	
		// Get color
		//vec3 SSR = textureLod(gPrevFrame, coords.xy, 0).rgb * clamp(ReflectionMultiplier, 0.0, 0.9) * 1.0;//Fresnel;
		vec3 SSR = texture(gPrevFrame, coords).rgb * clamp(ReflectionMultiplier, 0.0, 0.9) * 1.0;//Fresnel;
	
		float Le = length(texture(gPosition,coords).xyz - viewPos_);
		Le = clamp(Le * LLimiter, 0, 1);
		float error = 1 - Le;
		//SSR *= error;


		//SSR = mix(vec4(SSR,1.0), texture(gPrevFrame, coords)* clamp(ReflectionMultiplier, 0.0, 0.9) * 1.0,0.5).rgb;
		if(success==false){
			//SSR = texture(gPrevFrame, coords).rgb;// * clamp(ReflectionMultiplier, 0.0, 0.9) * 1.0;
		}

		reflectedColour = vec4(SSR,1.0);
	}
	else{
		reflectedColour = textureCube(skybox, fs_in.reflectedVector);
		refractedColour = vec4(0.0);
	
		refractedColour.r = textureCube(skybox, fs_in.refractedVectorR).r;
		refractedColour.g = textureCube(skybox, fs_in.refractedVectorG).g;
		refractedColour.b = textureCube(skybox, fs_in.refractedVectorB).b;
		refractedColour.a = 1.0;
	}
	//vec4 refractedColour = textureCube(skybox,fs_in.refractedVector);

	vec4 reflEnvWater2;
	
	if(isSSR)reflEnvWater2 = reflectedColour;
	else reflEnvWater2 = mix(reflectedColour, refractedColour, 0.9);

	//Fresnel mix equation (glass like material)
	float cos_theta =  abs(dot(-halfwayDir, normal));//abs(dot(normalize(-halfwayDir), normal));
    float cos_theta2 = cos_theta*cos_theta;// normal is normalized
    float r0 = 0.0204;
    float r = clamp(r0 +(1-r0)*cos_theta*cos_theta2*cos_theta2*cos_theta, 0.0, 1.0);


	vec4 reflEnvWater3 = mix(refractedColour, reflectedColour, r);
	

	if(water==1.0){
		
		vec3 albedo = waterColor.xyz/(1+shadow);//reflEnv.xyz;
		//vec3 albedo = lightColor/(1+shadow);
		float metallic =  0.250;//0.142;
		//0.2
		float roughness = 0.5;//0.255;//0.01;//0.285;
		float ao =  1;
		vec3 F0 = vec3(0.82); 
		F0 = mix(F0, albedo, metallic);
		
		//vec3 F = fresnelSchlick(cos_theta, F0);
		vec3 F  = fresnelSchlick(max(dot(H, V), 0.0), F0);

		
		float NDF = DistributionGGX(N, H, roughness);       
		float G   = GeometrySmith(N, V, L, roughness); 

		vec3 numerator    = NDF * G * F;
		float denominator = 4.0 * max(dot(N, V), 0.0) * max(dot(N, L), 0.0);
		vec3 specular     = numerator / max(denominator, 0.001); 
		

		vec3 kS = F;
		vec3 kD = vec3(1.0) - kS;
  
		kD *= 1.0 - metallic;	
  
		float NdotL = max(dot(N, L), 0.0);   
		vec3 Lo = vec3(0.0);     
		Lo += (kD * albedo / PI + specular) * radiance * NdotL;

		vec3 ambient_ = vec3(0.03) *albedo * ao;
		vec3 color_   = ( ambient_ + Lo);  

		
		//simple tone mapping
		color_ = color_ / (color_ + vec3(1.0));
		//color_ += reflEnv.xyz*0.2;
		
		color_ = pow(color_, vec3(1.0/2.2))*waterColor.xyz;  

		
		//lighting = pow(lighting, vec3(1.0/2.2));
		FragColor = vec4(lighting,1.0);
		FragColor = mix(FragColor, vec4(color_,1.0), 0.5);
		FragColor = mix(FragColor,reflEnvWater2,0.3);
		

		//FragColor = reflEnvWater2;
		
		
		//FragColor = vec4(lighting,1.0);
	    //FragColor = mix(FragColor,vec4(color_,colorVar.w),1);
		

		//FragColor = mix(vec4(lighting, 1.0), reflEnvWater2, 0.5);
		//only reflections
		//FragColor = reflEnvWater2;
		
		
		//FragColor = mix(FragColor,reflEnvWater,0.5);
		
		//FragColor.w = colorVar.w;
		//FragColor = pow(FragColor, vec4(1.0/2.2));
		//reflEnvWater.xyz *= color_;
		
		//reflEnvWater.xyz = pow(reflEnvWater.xyz, vec3(1.0/2.2));  

		//FragColor = reflEnvWater;
		//FragColor = mix(FragColor, reflEnvWater2, 0.5);
		 
		//FragColor = mix(FragColor, vec4(ambient_ + (1-shadow)*(diffuseP + specularP), 1.0),0.5);

		//FragColor.w = colorVar.w;
		
		
	}
	else if (water==2.0){
		
		vec3 albedo = lightingSphere/(1+shadow);//vec3(0.82, 0.11,0.11)/(1+shadow);//reflEnv.xyz;
		//albedo *= albedo;
		float metallic =  1;//0.142;
		//0.2
		float roughness = 0.5;//0.9;//0.01;//0.285;
		float ao = 1;
		vec3 F0 = vec3(0.47); 
		
		F0 = mix(F0, albedo, metallic);
		
		//vec3 F  = fresnelSchlick(max(dot(H, V), 0.0), F0);
		vec3 F  = fresnelSchlick(max(dot(H, V), 0.0), F0);
		
		//F = mix(F, reflEnvSphere.xyz, metallic);

		float NDF = DistributionGGX(N, H, roughness);       
		float G   = GeometrySmith(N, V, L, roughness); 

		vec3 numerator    = NDF * G * F;
		float denominator = 4.0 * max(dot(N, V), 0.0) * max(dot(N, L), 0.0);
		vec3 specular     = numerator / max(denominator, 0.001); 

		vec3 kS = F;
		vec3 kD = vec3(1.0) - kS;
  
		kD *= 1.0 - metallic;	
  
		float NdotL = max(dot(N, L), 0.0);   
		vec3 Lo = vec3(0.0);     
		Lo += (kD * albedo / PI + specular) * radiance * NdotL;

		//Lo = mix(Lo, reflEnvSphere.xyz, metallic);
		vec3 ambient_ = vec3(0.03)*albedo * ao;
		vec3 color_   =  (ambient_ + Lo);  

		

		//color_ *= ambient_+(1-shadow)*Lo;
		
		color_ = color_ / (color_ + vec3(1.0));
		//color_ += reflEnv.xyz*0.2;
		
		color_ = pow(color_, vec3(1.0/2.2));  
		lightingSphere = pow(lightingSphere, vec3(1.0/2.2));
		FragColor = vec4(lightingSphere,1.0);
		FragColor = mix(FragColor, vec4(color_,1.0), 0.5);
		//FragColor = mix(FragColor,reflEnvWater2,0.5);
		

		if(isSSR == false)FragColor = reflEnvWater2;
		
	}
	else {
		
		vec3 albedo = lighting;//reflEnv.xyz;
		float metallic =  1;//0.250;//0.142;
		//0.2
		float roughness = 0.75;//0.9;//0.01;//0.285;
		float ao =  1;
		vec3 F0 = vec3(0.03); 
		F0 = mix(F0, albedo, metallic);
		
		//vec3 F = fresnelSchlick(cos_theta, F0);
		vec3 F  = fresnelSchlick(max(dot(H, V), 0.0), F0);

		
		float NDF = DistributionGGX(N, H, roughness);       
		float G   = GeometrySmith(N, V, L, roughness); 

		vec3 numerator    = NDF * G * F;
		float denominator = 4.0 * max(dot(N, V), 0.0) * max(dot(N, L), 0.0);
		vec3 specular     = numerator / max(denominator, 0.001); 

		vec3 kS = F;
		vec3 kD = vec3(1.0) - kS;
  
		kD *= 1.0 - metallic;	
  
		float NdotL = max(dot(N, L), 0.0);   
		vec3 Lo = vec3(0.0);     
		Lo += (kD * albedo / PI + specular) * radiance * NdotL;

		vec3 ambient_ = vec3(0.23)*albedo * ao;
		vec3 color_   =  ambient_ + Lo;  

		

		//color_ *= ambient_+(1-shadow)*Lo;
		
		color_ = color_ / (color_ + vec3(1.0));
		//color_ += reflEnv.xyz*0.2;
		
		color_ = pow(color_, vec3(1.0/2.2));  

		//use of lighting
		FragColor = vec4(color_, 1.0);

		//no transparency to water
		//FragColor = vec4(lighting,1.0);

		//no lighting
		//FragColor = vec4(color, colorVar.w);

		//no lighting, no transparency
		//FragColor = vec4(color, 1.0);

		//if (coords.xy == vec2(-1.0))FragColor = mix(FragColor, vec4(ambient_ + (1-shadow)*(diffuseP + specularP), 1.0),0.5);
		//else FragColor = reflEnvWater2;

		if(isSSR) {
			FragColor = reflEnvWater2;
			//if (coords.xy == vec2(-1.0))FragColor = mix(FragColor, vec4(ambient_ + (1-shadow)*(diffuseP + specularP), 1.0),0.5);
			//else FragColor = reflEnvWater2;
		}
		else FragColor = mix(FragColor, vec4(ambient_ + (1-shadow)*(diffuseP + specularP), 1.0),0.5);

		//if(shadow != 0)FragColor = vec4(vec3(0),1);
	}

	// apply gamma correction
   
    //FragColor.rgb = pow(FragColor.rgb, vec3(1.0/2.2));
    //float depth = LinearizeDepth(gl_FragCoord.z) / far;
    //FragColor = vec4(vec3(depth), 1.0);
}