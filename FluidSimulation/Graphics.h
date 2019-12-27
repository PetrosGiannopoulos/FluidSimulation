#pragma once

#include <glad\glad.h>
#include <GLFW\glfw3.h>

#include "Shader.h"
#include "CShader.h"

#include "stb_image.h"

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

#include <iostream>

#include "Camera.h"
#include "Models.h"
#include "MarchingCubes.h"
#include "CollisionObject.h"

#include <omp.h>
#include <CL/cl.hpp>

#include <thread>

class Graphics {

public:
	GLFWwindow *window;
	unsigned int VBO, planeVAO, lightVAO, VBO2, EBO;
	unsigned int skyboxVAO, skyboxVBO;
	unsigned int screenQuadVAO, screenQuadVBO;
	unsigned int cubemapTexture;


	unsigned int instanceVBO;

	Shader shader, lampShader,simpleDepthShader;
	Shader skyboxShader,reflSkyboxShader;
	Shader gBufferShader;
	Shader screenShader;
	Shader ssaoShader, ssaoBlurShader;
	CShader computeShader;
	unsigned int planeTexture;
	unsigned int planeNormalMap;
	unsigned int waterTexture;
	unsigned int depthMapFBO;
	unsigned int depthMap;
	unsigned int depthCubemap;

	unsigned int gBuffer;
	unsigned int gPosition;
	unsigned int gNormal;
	unsigned int gAlbedoSpec;
	unsigned int gColor;

	//SSAO
	unsigned int ssaoFBO, ssaoBlurFBO;
	std::vector<glm::vec3> ssaoKernel;
	std::vector<glm::vec3> ssaoNoise;
	unsigned int noiseTexture;
	//SSAO textures
	unsigned int ssaoColorBuffer, ssaoColorBufferBlur;

	Vector4 clippingPlane;

	// lighting
	glm::vec3 lightPos;
	glm::mat4 lightProjection;

	// camera
	Camera camera;
	float lastX = 1000 / 2.0f;
	float lastY = 800 / 2.0f;
	bool firstMouse = true;

	// timing
	float deltaTime = 0.0f;	// time between current frame and last frame
	float lastFrame = 0.0f;

	SphGrid *grid;

	bool shadows = true;
	bool shadowsKeyPressed = false;

	bool isIrridescent = false;

	unsigned int diffuseMap, specularMap;
	const unsigned int SHADOW_WIDTH = 1024, SHADOW_HEIGHT = 1024;

	struct Vertex {
		// position
		Vector3 Position;
		// normal
		Vector3 Normal;
		// texCoords
		Vector2 TexCoords;
	};

	glm::vec3 lightPositions[4] = {
		glm::vec3(-3.0f, 0.0f, 0.0f),
		glm::vec3(-1.0f, 0.0f, 0.0f),
		glm::vec3(1.0f, 0.0f, 0.0f),
		glm::vec3(3.0f, 0.0f, 0.0f)
	};
	glm::vec3 lightColors[4] = {
		glm::vec3(0.25),
		glm::vec3(0.50),
		glm::vec3(0.75),
		glm::vec3(1.00)
	};

	Models model;

	MarchingCubes mcCubes;

	bool startSim = false;
	bool startCollSim = false;
	
	CollisionObject collisionSphere;
	CollisionObject collisionWood;
	CollisionObject collisionCube;
	CollisionObject collisionMirror;

	int cubemapWidth, cubemapHeight;
	unsigned int dcmap, framebuffer, prevFrameBuffer, texColorBuffer, copyBuffer;
	unsigned int noSphereMap, noSphereFrameBuffer;
	unsigned int dcRenderBuffer, noSphereRenderBuffer;
	unsigned int depthTexture;

	unsigned int computeShaderTexture;

	unsigned int ssbo_sphere;
	
	bool isdcmapDrawn = false;
	glm::mat4 model_;

	//1024
	const unsigned int cubemapSize = 2048;

	vector<CollisionObject> collObjects;

	float rotX = 0;
	float rotY = 0;
	float rotZ = 0;

	float tX = 10, tY=-17,tZ=8;

	glm::vec3 pCM;

	float scaleAway = 0;

	float camDist1;
	float camDist2;

	glm::vec3 cubemapCenter1;
	glm::vec3 cubemapCenter2;

	glm::vec3 bboxmin1, bboxmin2;
	glm::vec3 bboxmax1, bboxmax2;

	glm::vec3 objectDir;

	bool executeMC;
	bool drawSphere;
	bool useSSR;
	bool useDirectRaytracing;
	bool renderWater;

	glm::vec3 customModelPos;

	float fBias = 0.74, fScale = 2.65, fPower = -1.34;

	float gScattering = 0;
	float scatteringAmount = 0;
	
public:

	~Graphics() {
		delete grid;
		delete this;
	}

	float lerp(float a, float b, float f)
	{
		return a + f * (b - a);
	}

	Graphics(SphGrid *grid) {

		mcCubes = MarchingCubes(grid);
		mcCubes.setupCL();
		//setupCL();

		//this->grid = grid;
		string runPlatformCollisionTest;
		cout << "Do you want to run platform collision test (Y/N)?";
		cin >> runPlatformCollisionTest;
		cout << endl;

		string ssrString;
		cout << "Do you want to perform Screen Space Reflection (Y/N)?";
		cin >> ssrString;
		cout << endl;

		if (ssrString == "Y" || ssrString == "y") {
			useSSR = true;
		}
		else useSSR = false;

		if (useSSR == false) {
			string raytracingString;
			cout << "Do you want to perform direct raytracing (Y/N)?";
			cin >> raytracingString;
			cout << endl;

			if (raytracingString == "y" || raytracingString == "Y")useDirectRaytracing = true;
			else useDirectRaytracing = false;
		}
		else {
			useDirectRaytracing = false;
		}

		if (useDirectRaytracing)renderWater = false;
		else renderWater = true;
		//renderWater = true;

		/*if (useDirectRaytracing == true) {

			grid->vars.kXMax = 5;
			grid->vars.kYMax = 5;
			grid->vars.kZMax = 5;
			grid->vars.kNumParticles = 100;
			grid->vars.recalculate();
			grid->recalculate();
		}
		*/

		glfwInit();
		glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
		glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 6);
		glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
		//glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
		glfwWindowHint(GLFW_SAMPLES, 4);
		
		

		window = glfwCreateWindow(1000, 800, "Fluid Simulation OpenGL", NULL, NULL);
		if (window == NULL)
		{
			std::cout << "Failed to create GLFW window" << std::endl;
			glfwTerminate();
		}
		glfwMakeContextCurrent(window);
		
		glfwSetWindowPos(window, 800, 200);

		glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_DISABLED);
		if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress))
		{
			std::cout << "Failed to initialize GLAD" << std::endl;

		}
		//glEnable(GL_MULTISAMPLE);
		glEnable(GL_DEPTH_TEST);
		//glEnable(GL_CULL_FACE);
		glViewport(0, 0, 1000, 800);

		//glDisable(GL_CULL_FACE);

		shader = Shader("vertexShader.vs", "fragmentShader.fs");
		lampShader = Shader("lampShader.vs", "lampShader.fs");
		simpleDepthShader = Shader("shadow_mapping.vs","shadow_mapping.fs", "shadow_mapping.gs");
		skyboxShader = Shader("skybox.vert","skybox.frag");
		reflSkyboxShader = Shader("skyboxrefl.vert", "skybox.frag");
		gBufferShader = Shader("gBuffer.vs", "gBuffer.fs");
		screenShader = Shader("screenShader.vs", "screenShader.fs");
		ssaoShader = Shader("ssao.vs", "ssao.fs");
		ssaoBlurShader = Shader("ssaoBlur.vs", "ssaoBlur.fs");


		computeShader = CShader("computeShader.glslcs");
		computeShaderTexture = computeShader.createFrameBufferTexture();
		
		computeShader.use();
		computeShader.setInt("skybox",0);
		computeShader.setInt("diffuseTexture",1);
		computeShader.setInt("ssao", 2);
		computeShader.setInt("shadowMap", 3);

		//glLineWidth(5);
		glEnable(GL_BLEND);
		
		//glEnable(GL_TEXTURE_CUBE_MAP_SEAMLESS);
		
		//glBlendFunc(GL_DST_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
		//glBlendFunc(GL_DST_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
		glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

		//glEnable(GL_CLIP_DISTANCE0);

		//cout << glGetString(GL_VERSION) << endl;

		float lvertices[] = {
			-0.5f, -0.5f, -0.5f,  0.0f,  0.0f, -1.0f,
			0.5f, -0.5f, -0.5f,  0.0f,  0.0f, -1.0f,
			0.5f,  0.5f, -0.5f,  0.0f,  0.0f, -1.0f,
			0.5f,  0.5f, -0.5f,  0.0f,  0.0f, -1.0f,
			-0.5f,  0.5f, -0.5f,  0.0f,  0.0f, -1.0f,
			-0.5f, -0.5f, -0.5f,  0.0f,  0.0f, -1.0f,

			-0.5f, -0.5f,  0.5f,  0.0f,  0.0f,  1.0f,
			0.5f, -0.5f,  0.5f,  0.0f,  0.0f,  1.0f,
			0.5f,  0.5f,  0.5f,  0.0f,  0.0f,  1.0f,
			0.5f,  0.5f,  0.5f,  0.0f,  0.0f,  1.0f,
			-0.5f,  0.5f,  0.5f,  0.0f,  0.0f,  1.0f,
			-0.5f, -0.5f,  0.5f,  0.0f,  0.0f,  1.0f,

			-0.5f,  0.5f,  0.5f, -1.0f,  0.0f,  0.0f,
			-0.5f,  0.5f, -0.5f, -1.0f,  0.0f,  0.0f,
			-0.5f, -0.5f, -0.5f, -1.0f,  0.0f,  0.0f,
			-0.5f, -0.5f, -0.5f, -1.0f,  0.0f,  0.0f,
			-0.5f, -0.5f,  0.5f, -1.0f,  0.0f,  0.0f,
			-0.5f,  0.5f,  0.5f, -1.0f,  0.0f,  0.0f,

			0.5f,  0.5f,  0.5f,  1.0f,  0.0f,  0.0f,
			0.5f,  0.5f, -0.5f,  1.0f,  0.0f,  0.0f,
			0.5f, -0.5f, -0.5f,  1.0f,  0.0f,  0.0f,
			0.5f, -0.5f, -0.5f,  1.0f,  0.0f,  0.0f,
			0.5f, -0.5f,  0.5f,  1.0f,  0.0f,  0.0f,
			0.5f,  0.5f,  0.5f,  1.0f,  0.0f,  0.0f,

			-0.5f, -0.5f, -0.5f,  0.0f, -1.0f,  0.0f,
			0.5f, -0.5f, -0.5f,  0.0f, -1.0f,  0.0f,
			0.5f, -0.5f,  0.5f,  0.0f, -1.0f,  0.0f,
			0.5f, -0.5f,  0.5f,  0.0f, -1.0f,  0.0f,
			-0.5f, -0.5f,  0.5f,  0.0f, -1.0f,  0.0f,
			-0.5f, -0.5f, -0.5f,  0.0f, -1.0f,  0.0f,

			-0.5f,  0.5f, -0.5f,  0.0f,  1.0f,  0.0f,
			0.5f,  0.5f, -0.5f,  0.0f,  1.0f,  0.0f,
			0.5f,  0.5f,  0.5f,  0.0f,  1.0f,  0.0f,
			0.5f,  0.5f,  0.5f,  0.0f,  1.0f,  0.0f,
			-0.5f,  0.5f,  0.5f,  0.0f,  1.0f,  0.0f,
			-0.5f,  0.5f, -0.5f,  0.0f,  1.0f,  0.0f
		};

		float vertices[] = {
			// positions            // normals         // texcoords
			10.0f, -0.5f,  10.0f,  0.0f, 1.0f, 0.0f,  10.0f,  0.0f,
			-10.0f, -0.5f,  10.0f,  0.0f, 1.0f, 0.0f,   0.0f,  0.0f,
			-10.0f, -0.5f, -10.0f,  0.0f, 1.0f, 0.0f,   0.0f, 10.0f,

			10.0f, -0.5f,  10.0f,  0.0f, 1.0f, 0.0f,  10.0f,  0.0f,
			-10.0f, -0.5f, -10.0f,  0.0f, 1.0f, 0.0f,   0.0f, 10.0f,
			10.0f, -0.5f, -10.0f,  0.0f, 1.0f, 0.0f,  10.0f, 10.0f
		};

		float skyboxVertices[] = {
			// positions          
			-1.0f,  1.0f, -1.0f,
			-1.0f, -1.0f, -1.0f,
			1.0f, -1.0f, -1.0f,
			1.0f, -1.0f, -1.0f,
			1.0f,  1.0f, -1.0f,
			-1.0f,  1.0f, -1.0f,

			-1.0f, -1.0f,  1.0f,
			-1.0f, -1.0f, -1.0f,
			-1.0f,  1.0f, -1.0f,
			-1.0f,  1.0f, -1.0f,
			-1.0f,  1.0f,  1.0f,
			-1.0f, -1.0f,  1.0f,

			1.0f, -1.0f, -1.0f,
			1.0f, -1.0f,  1.0f,
			1.0f,  1.0f,  1.0f,
			1.0f,  1.0f,  1.0f,
			1.0f,  1.0f, -1.0f,
			1.0f, -1.0f, -1.0f,

			-1.0f, -1.0f,  1.0f,
			-1.0f,  1.0f,  1.0f,
			1.0f,  1.0f,  1.0f,
			1.0f,  1.0f,  1.0f,
			1.0f, -1.0f,  1.0f,
			-1.0f, -1.0f,  1.0f,

			-1.0f,  1.0f, -1.0f,
			1.0f,  1.0f, -1.0f,
			1.0f,  1.0f,  1.0f,
			1.0f,  1.0f,  1.0f,
			-1.0f,  1.0f,  1.0f,
			-1.0f,  1.0f, -1.0f,

			-1.0f, -1.0f, -1.0f,
			-1.0f, -1.0f,  1.0f,
			1.0f, -1.0f, -1.0f,
			1.0f, -1.0f, -1.0f,
			-1.0f, -1.0f,  1.0f,
			1.0f, -1.0f,  1.0f
		};

		// skybox VAO
		
		glGenVertexArrays(1, &skyboxVAO);
		glGenBuffers(1, &skyboxVBO);
		glBindVertexArray(skyboxVAO);
		glBindBuffer(GL_ARRAY_BUFFER, skyboxVBO);
		glBufferData(GL_ARRAY_BUFFER, sizeof(skyboxVertices), &skyboxVertices, GL_STATIC_DRAW);
		glEnableVertexAttribArray(0);
		glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);

		// load textures
		// -------------
		vector<std::string> faces
		{
			"Textures/s_right.jpg",
			"Textures/s_left.jpg",
			"Textures/s_top.jpg",
			"Textures/s_bottom.jpg",
			"Textures/s_front.jpg",
			"Textures/s_back.jpg",

		};
		cubemapTexture = loadCubemap(faces);

		glGenVertexArrays(1, &planeVAO);
		glGenBuffers(1, &VBO);
		glGenBuffers(1, &VBO2);
		//glGenBuffers(1, &EBO);


		glBindVertexArray(planeVAO);
		glBindBuffer(GL_ARRAY_BUFFER, VBO);
		glBufferData(GL_ARRAY_BUFFER, sizeof(vertices), vertices, GL_STATIC_DRAW);

		/*glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
		glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(indices), indices, GL_STATIC_DRAW);*/

		
		glGenBuffers(1, &instanceVBO);


		// position attribute
		glEnableVertexAttribArray(0);
		glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 8 * sizeof(float), (void*)0);

		//normal attribute
		glEnableVertexAttribArray(1);
		glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 8 * sizeof(float), (void*)(3 * sizeof(float)));

		// texture coord attribute
		glEnableVertexAttribArray(2);
		glVertexAttribPointer(2, 2, GL_FLOAT, GL_FALSE, 8 * sizeof(float), (void*)(6 * sizeof(float)));

		glBindVertexArray(0);
		//light
		glGenVertexArrays(1, &lightVAO);
		glBindVertexArray(lightVAO);

		// we only need to bind to the VBO (to link it with glVertexAttribPointer), no need to fill it; the VBO's data already contains all we need (it's already bound, but we do it again for educational purposes)
		glBindBuffer(GL_ARRAY_BUFFER, VBO2);
		glBufferData(GL_ARRAY_BUFFER, sizeof(lvertices), &lvertices, GL_STATIC_DRAW);

		glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void*)0);
		glEnableVertexAttribArray(0);


		planeTexture = loadTexture("Textures/floor.jpg");
		planeNormalMap = loadTexture("Textures/floorNormal.jpg");
		waterTexture = loadTexture("Textures/TexturesCom_WaterFoam0001_2_S.jpg");
		// load and create a texture 
		// -------------------------
		//diffuse Map: texture with light
		//diffuseMap = loadTexture("Textures/WoodFlooring044_COL_4K.jpg");
		//specularMap = loadTexture("Textures/WoodFlooring044_REFL_4K.jpg");

		// configure depth map FBO
		// -----------------------
		glGenFramebuffers(1, &depthMapFBO);
		// create depth cubemap texture
		
		glGenTextures(1, &depthCubemap);
		glBindTexture(GL_TEXTURE_CUBE_MAP, depthCubemap);
		for (unsigned int i = 0; i < 6; ++i)
			glTexImage2D(GL_TEXTURE_CUBE_MAP_POSITIVE_X + i, 0, GL_DEPTH_COMPONENT, SHADOW_WIDTH, SHADOW_HEIGHT, 0, GL_DEPTH_COMPONENT, GL_FLOAT, NULL);
		glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
		glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
		glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
		glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
		glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_WRAP_R, GL_CLAMP_TO_EDGE);
		// attach depth texture as FBO's depth buffer
		glBindFramebuffer(GL_FRAMEBUFFER, depthMapFBO);
		glFramebufferTexture(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, depthCubemap, 0);
		glDrawBuffer(GL_NONE);
		glReadBuffer(GL_NONE);
		glBindFramebuffer(GL_FRAMEBUFFER, 0);

		//configure g-buffer
		glGenFramebuffers(1, &gBuffer);
		glBindFramebuffer(GL_FRAMEBUFFER, gBuffer);
		// position color buffer
		glGenTextures(1, &gPosition);
		glBindTexture(GL_TEXTURE_2D, gPosition);
		glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB16F, 1000, 800, 0, GL_RGB, GL_FLOAT, NULL);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
		glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, gPosition, 0);
		// normal color buffer
		glGenTextures(1, &gNormal);
		glBindTexture(GL_TEXTURE_2D, gNormal);
		
		glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB16F, 1000, 800, 0, GL_RGB, GL_FLOAT, NULL);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
		glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT1, GL_TEXTURE_2D, gNormal, 0);
		// color + specular color buffer
		glGenTextures(1, &gAlbedoSpec);
		glBindTexture(GL_TEXTURE_2D, gAlbedoSpec);
		glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, 1000, 800, 0, GL_RGBA, GL_UNSIGNED_BYTE, NULL);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
		glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT2, GL_TEXTURE_2D, gAlbedoSpec, 0);

		glGenTextures(1, &depthTexture);
		glBindTexture(GL_TEXTURE_2D, depthTexture);
		glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB32F, 1000, 800, 0, GL_RGB, GL_UNSIGNED_BYTE, NULL);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
		glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT3, GL_TEXTURE_2D, depthTexture, 0);

		// tell OpenGL which color attachments we'll use (of this framebuffer) for rendering 
		unsigned int attachments[4] = { GL_COLOR_ATTACHMENT0, GL_COLOR_ATTACHMENT1, GL_COLOR_ATTACHMENT2,GL_COLOR_ATTACHMENT3};
		glDrawBuffers(4, attachments);
		// create and attach depth buffer (renderbuffer)
		unsigned int rboDepth;
		glGenRenderbuffers(1, &rboDepth);
		glBindRenderbuffer(GL_RENDERBUFFER, rboDepth);
		glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH_COMPONENT, 1000, 800);
		glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_RENDERBUFFER, rboDepth);

		
		/*glGenTextures(1, &depthTexture);
		glBindTexture(GL_TEXTURE_2D, depthTexture);
		glTexImage2D(GL_TEXTURE_2D, 0, GL_DEPTH_COMPONENT32, 1000, 800, 0, GL_DEPTH_COMPONENT, GL_FLOAT, NULL);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
		glFramebufferTexture2D(GL_DRAW_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_TEXTURE_2D, depthTexture, 0);
		*/


		// finally check if framebuffer is complete
		if (glCheckFramebufferStatus(GL_FRAMEBUFFER) != GL_FRAMEBUFFER_COMPLETE)
			std::cout << "Framebuffer not complete!" << std::endl;
		glBindFramebuffer(GL_FRAMEBUFFER, 0);

		
		glGenFramebuffers(1, &prevFrameBuffer);
		glBindFramebuffer(GL_FRAMEBUFFER, prevFrameBuffer);

		glGenTextures(1, &texColorBuffer);
		glBindTexture(GL_TEXTURE_2D, texColorBuffer);
		glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, 1000, 800, 0, GL_RGB, GL_UNSIGNED_BYTE, NULL);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
		glBindTexture(GL_TEXTURE_2D, 0);
		// attach it to currently bound framebuffer object
		glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, texColorBuffer, 0);
		
		unsigned int prbo;
		glGenRenderbuffers(1, &prbo);
		glBindRenderbuffer(GL_RENDERBUFFER, prbo);
		glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH24_STENCIL8, 1000, 800);
		glBindRenderbuffer(GL_RENDERBUFFER, 0);
		glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_STENCIL_ATTACHMENT, GL_RENDERBUFFER, prbo);

		

		if (glCheckFramebufferStatus(GL_FRAMEBUFFER) != GL_FRAMEBUFFER_COMPLETE)
			std::cout << "ERROR::FRAMEBUFFER:: Framebuffer is not complete!" << std::endl;
		glBindFramebuffer(GL_FRAMEBUFFER, 0);


		glGenFramebuffers(1, &copyBuffer);
		glBindFramebuffer(GL_FRAMEBUFFER, copyBuffer);

		glGenTextures(1, &gColor);
		glBindTexture(GL_TEXTURE_2D, gColor);
		glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, 1000, 800, 0, GL_RGB, GL_UNSIGNED_BYTE, NULL);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
		glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, gColor, 0);

		unsigned int cprbo;
		glGenRenderbuffers(1, &cprbo);
		glBindRenderbuffer(GL_RENDERBUFFER, cprbo);
		glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH24_STENCIL8, 1000, 800);
		glBindRenderbuffer(GL_RENDERBUFFER, 0);
		glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_STENCIL_ATTACHMENT, GL_RENDERBUFFER, cprbo);



		if (glCheckFramebufferStatus(GL_FRAMEBUFFER) != GL_FRAMEBUFFER_COMPLETE)
			std::cout << "ERROR::FRAMEBUFFER:: Framebuffer is not complete!" << std::endl;
		glBindFramebuffer(GL_FRAMEBUFFER, 0);

		//SSAO
		glGenFramebuffers(1, &ssaoFBO);  glGenFramebuffers(1, &ssaoBlurFBO);
		glBindFramebuffer(GL_FRAMEBUFFER, ssaoFBO);
		
		// SSAO color buffer
		glGenTextures(1, &ssaoColorBuffer);
		glBindTexture(GL_TEXTURE_2D, ssaoColorBuffer);
		glTexImage2D(GL_TEXTURE_2D, 0, GL_RED, 1000, 800, 0, GL_RGB, GL_FLOAT, NULL);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
		glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, ssaoColorBuffer, 0);
		if (glCheckFramebufferStatus(GL_FRAMEBUFFER) != GL_FRAMEBUFFER_COMPLETE)
			std::cout << "SSAO Framebuffer not complete!" << std::endl;
		// and blur stage
		glBindFramebuffer(GL_FRAMEBUFFER, ssaoBlurFBO);
		glGenTextures(1, &ssaoColorBufferBlur);
		glBindTexture(GL_TEXTURE_2D, ssaoColorBufferBlur);
		glTexImage2D(GL_TEXTURE_2D, 0, GL_RED, 1000, 800, 0, GL_RGB, GL_FLOAT, NULL);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
		glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, ssaoColorBufferBlur, 0);
		if (glCheckFramebufferStatus(GL_FRAMEBUFFER) != GL_FRAMEBUFFER_COMPLETE)
			std::cout << "SSAO Blur Framebuffer not complete!" << std::endl;
		glBindFramebuffer(GL_FRAMEBUFFER, 0);

		// generate sample kernel
		// ----------------------
		std::uniform_real_distribution<GLfloat> randomFloats(0.0, 1.0); // generates random floats between 0.0 and 1.0
		std::default_random_engine generator;
		for (unsigned int i = 0; i < 64; ++i)
		{
			glm::vec3 sample(randomFloats(generator) * 2.0 - 1.0, randomFloats(generator) * 2.0 - 1.0, randomFloats(generator));
			sample = glm::normalize(sample);
			sample *= randomFloats(generator);
			float scale = float(i) / 64.0;

			// scale samples s.t. they're more aligned to center of kernel
			scale = lerp(0.1f, 1.0f, scale * scale);
			sample *= scale;
			ssaoKernel.push_back(sample);
		}

		for (unsigned int i = 0; i < 16; i++)
		{
			glm::vec3 noise(randomFloats(generator) * 2.0 - 1.0, randomFloats(generator) * 2.0 - 1.0, 0.0f); // rotate around z-axis (in tangent space)
			ssaoNoise.push_back(noise);
		}
		glGenTextures(1, &noiseTexture);
		glBindTexture(GL_TEXTURE_2D, noiseTexture);
		glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB32F, 4, 4, 0, GL_RGB, GL_FLOAT, &ssaoNoise[0]);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);


		float screenQuadVertices[] = { 
			-1.0f,  1.0f,  0.0f, 1.0f,
			-1.0f, -1.0f,  0.0f, 0.0f,
			1.0f, -1.0f,  1.0f, 0.0f,

			-1.0f,  1.0f,  0.0f, 1.0f,
			1.0f, -1.0f,  1.0f, 0.0f,
			1.0f,  1.0f,  1.0f, 1.0f
		};

		
		glGenVertexArrays(1, &screenQuadVAO);
		glGenBuffers(1, &screenQuadVBO);
		glBindVertexArray(screenQuadVAO);
		glBindBuffer(GL_ARRAY_BUFFER, screenQuadVBO);
		glBufferData(GL_ARRAY_BUFFER, sizeof(screenQuadVertices), &screenQuadVertices, GL_STATIC_DRAW);
		glEnableVertexAttribArray(0);
		glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 4 * sizeof(float), (void*)0);
		glEnableVertexAttribArray(1);
		glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, 4 * sizeof(float), (void*)(2 * sizeof(float)));

		// shader configuration
		// --------------------
		shader.use();
		shader.setInt("diffuseTexture", 0);
		shader.setInt("depthMap", 1);
		shader.setInt("skybox", 2);
		shader.setInt("normalMap", 3);

		shader.setInt("gPosition",4);
		shader.setInt("gNormal", 5);
		shader.setInt("gAlbedoSpec",6);
		shader.setInt("gPrevFrame", 7);
		shader.setInt("gDepth", 8);
		shader.setInt("ssao", 9);

		ssaoShader.use();
		ssaoShader.setInt("gPosition", 0);
		ssaoShader.setInt("gNormal", 1);
		ssaoShader.setInt("texNoise", 2);
		ssaoBlurShader.use();
		ssaoBlurShader.setInt("ssaoInput", 0);
		
		skyboxShader.use();
		skyboxShader.setInt("skybox",0);
		
		screenShader.use();
		screenShader.setInt("screenTexture", 0);

		

		camera = Camera(glm::vec3(10.0f, 8.0f, 20.0f));
		camera.Pitch -= 15;
		camera.updateCameraVectors();
		
		if(useDirectRaytracing==false)lightPos = glm::vec3(10.0f, 16.0f, -5.0f);
		else lightPos = glm::vec3(10.0f, 30.0f, -5.0f);

		model = Models();
		model.loadModel("icosahedron.obj");
		model.refineModel(2);
		//model.redefineStructure();
		model.rebuildStructure();
		model.rescaleModel(grid->particleSize);

		grid->setTriangleData(model.m_triangles);
		//grid->setData(model.m_vertices, model.m_indices);
		grid->setData(model.m_vertices, model.m_indices, model.m_vertexes);
		//grid->setupParticles();
		//grid->setupParticlesInstanced(instanceVBO);

		this->grid = grid;
		
		grid->particles[0].setupBuffersInstanced();
		//mcCubes.setBuffersCL();
		//mcCubes.setupRenderBuffers();

		mcCubes.setupGridCells();

		//Sphere
		//cout << glGetString(GL_VERSION) << endl;
		collisionSphere = CollisionObject(grid);
		Models collModels = Models();
		collModels.loadModel("icosahedron.obj");
		
		if(useDirectRaytracing==false)collModels.refineModel(7);
		//model.redefineStructure();
		collModels.rebuildStructure();
		collModels.rescaleModel(collisionSphere.radius);
		
		collModels.fillSSBO();

		//collModels.smoothNormals();
		//collModels.vertexNormals();

		collisionSphere.type = "Sphere";
		collisionSphere.setData(collModels.m_vertices, collModels.m_indices, collModels.m_vertexes);
		collisionSphere.setTriangleData(collModels.m_triangles, collModels.ssbo_triangles);
		collisionSphere.setupBuffers();

		if (useDirectRaytracing == false && useSSR == false) {
			dcmap = initDynamicCubemap();
			noSphereMap = initDynamicNoSphereMap();
		}

		collisionSphere.id = 2;
		collisionSphere.bounded = false;
		
		customModelPos = collisionSphere.getPos();
		customModelPos.x -= 10;

		//wood
		collisionWood = CollisionObject(grid);
		
		collModels.clearAll();
		collModels.loadModel("wood.obj");
		//collModels.refineModel(2);
		//model.redefineStructure();
		collModels.rebuildStructure();
		collModels.rescaleModel(10);

		collModels.fillSSBO();

		collisionWood.setData(collModels.m_vertices, collModels.m_indices, collModels.m_vertexes);
		collisionWood.setTriangleData(collModels.m_triangles,collModels.ssbo_triangles);
		collisionWood.setupBuffers();

		collisionWood.id = 2;
		//collisionWood.texture = planeTexture;
		collisionWood.bounded = true;
		collisionWood.type = "Random";
		collisionWood.pos += Vector3(10,-10,-10);

		//cube
		collisionCube = CollisionObject(grid);

		collModels.clearAll();
		collModels.loadModel("cube.obj");
		//collModels.refineModel(2);
		//model.redefineStructure();
		collModels.rebuildStructure();
		//collModels.rescaleModel(10);

		collModels.ssbo_triangles;

		collisionCube.type = "Cube";
		collisionCube.setData(collModels.m_vertices, collModels.m_indices, collModels.m_vertexes);
		collisionCube.setTriangleData(collModels.m_triangles, collModels.ssbo_triangles);
		collisionCube.setupBuffers();

		collisionCube.id = 2;
		collisionCube.DRT = useDirectRaytracing;
		//collisionCube.texture = planeTexture;
		
		collisionCube.cubeSize = collisionCube.getCubeSize();
		//collisionCube.pos += Vector3(10, -10, -10);

		//mirror
		collisionMirror = CollisionObject(grid);

		collModels.clearAll();
		collModels.loadModel("mirror.obj");
		//collModels.refineModel(2);
		//model.redefineStructure();
		collModels.rebuildStructure();
		//collModels.rescaleModel(10);

		collModels.fillSSBO();

		collisionMirror.type = "Plane";
		collisionMirror.setData(collModels.m_vertices, collModels.m_indices, collModels.m_vertexes);
		collisionMirror.setTriangleData(collModels.m_triangles,collModels.ssbo_triangles);
		collisionMirror.setupBuffers();

		collisionMirror.id = 2;
		collisionMirror.texture = planeTexture;
		


		collisionSphere.collId = 0;
		collisionWood.collId = 1;
		collisionCube.collId = 2;
		collisionMirror.collId = 3;
		
		collisionSphere.pos += Vector3(5,0,0);
		collisionMirror.pos += Vector3(12,0,0);

		collObjects.push_back(collisionSphere);
		collObjects.push_back(collisionCube);
		collObjects.push_back(collisionMirror);
		
		if(runPlatformCollisionTest == "y" || runPlatformCollisionTest == "Y")collObjects.push_back(collisionWood);
		
		for (int i = 0; i < collObjects.size();i++) {
			collObjects[i].setCollObjects(&collObjects);
		}
		
		//remote distance
		/*camera.Position = glm::vec3(69.8066, 58.8452, 70.7498);
		camera.Yaw = -122.4;
		camera.Pitch = -27.6;
		camera.Zoom = 45.0f;
		*/

		mcCubes.setupNormalPointLines();
		
		pCM = glm::vec3(0,0,0);

		objectDir = glm::normalize(collObjects[0].getPos() - collObjects[1].getPos());

		if(useDirectRaytracing)initSSBOData();

		//cout << collObjects[0].m_vertices.size() << endl;
		
	}

	void initSSBOData() {

		//GL_SHADER_STORAGE_BUFFER
		//compute Shader Buffer Data
		glCreateBuffers(1, &ssbo_sphere);
		glBindBuffer(GL_SHADER_STORAGE_BUFFER, ssbo_sphere);
		
		glNamedBufferData(ssbo_sphere, collObjects[0].ssbo_triangles.size() * sizeof(Models::TriangleSSBO), collObjects[0].ssbo_triangles.data(), GL_STATIC_DRAW);

		//glNamedBufferData(ssbo_sphere, mcCubes.triangulator.mesh.ssbo_triangles.size() * sizeof(Models::TriangleSSBO), mcCubes.triangulator.mesh.ssbo_triangles.data(), GL_STATIC_DRAW);
		//glNamedBufferData(ssbo_sphere, 1000 * sizeof(Models::TriangleSSBO), NULL, GL_DYNAMIC_DRAW);
		
		//glBufferData(GL_SHADER_STORAGE_BUFFER, collObjects[0].ssbo_triangles.size() * sizeof(Models::TriangleSSBO), collObjects[0].ssbo_triangles.data(), GL_DYNAMIC_COPY);
		glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 1, ssbo_sphere);//binding number = 1
		glBindBuffer(GL_SHADER_STORAGE_BUFFER, 0); // unbind

		

		
	}

	void renderSSBOData(unsigned int ssbo, const GLchar* name, vector<Models::TriangleSSBO> ssbo_tris) {
		//Only needed if another buffer has been bound to the binding point 0
		//glBindBuffer(GL_SHADER_STORAGE_BUFFER, ssbo);
		
		glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 1, ssbo);
		//Only needed if the buffer data needs to be changed
		//glNamedBufferData(ssbo_sphere, ssbo_tris.size() * sizeof(Models::TriangleSSBO), ssbo_tris.data(), GL_DYNAMIC_DRAW);

		//glBufferSubData(GL_SHADER_STORAGE_BUFFER, 0, sizeof(unsigned int), (unsigned int)ssbo_tris.size());
		glBufferSubData(GL_SHADER_STORAGE_BUFFER, 0, ssbo_tris.size() * sizeof(Models::TriangleSSBO), &(ssbo_tris[0]));

		

		//If you use the uniform option
		glUniform1i(glGetUniformLocation(computeShader.ID, name), ssbo_tris.size());
	}

	void unBindSSBO() {
		glBindBuffer(GL_SHADER_STORAGE_BUFFER, 0);
	}

	void updateSSBOData(unsigned int ssbo_update, vector<Models::TriangleSSBO> &ssbo_tris, int bindingIndex, const GLchar* name) {
		glBindBuffer(GL_SHADER_STORAGE_BUFFER, ssbo_update);
		struct Models::TriangleSSBO* p = (struct Models::TriangleSSBO*)glMapBufferRange(GL_SHADER_STORAGE_BUFFER, 0, ssbo_tris.size() * sizeof(Models::TriangleSSBO), GL_MAP_WRITE_BIT | GL_MAP_INVALIDATE_BUFFER_BIT);
		//memcpy(p, &shader_data, sizeof(shader_data));
		//i.e. 
		memcpy(p, ssbo_tris.data(), ssbo_tris.size()*sizeof(Models::TriangleSSBO));
		glUnmapBuffer(GL_SHADER_STORAGE_BUFFER);
		glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 1, ssbo_update);

		unsigned int block_index = 0;
		block_index = glGetProgramResourceIndex(computeShader.ID, GL_SHADER_STORAGE_BLOCK, name);

		unsigned int ssbo_binding_point_index = 1;
		glShaderStorageBlockBinding(computeShader.ID, block_index, ssbo_binding_point_index);
	}
	
	unsigned int createFrameBuffer() {
		unsigned int frameBuffer;
		glGenFramebuffers(1, &frameBuffer);
		glBindFramebuffer(GL_FRAMEBUFFER,frameBuffer);
		glDrawBuffer(GL_COLOR_ATTACHMENT0);
		return frameBuffer;
		
	}

	unsigned int createTextureAttachment(int width, int height) {
		unsigned int textureAttach;
		glGenTextures(1, &textureAttach);
		glBindTexture(GL_TEXTURE_2D, textureAttach);

		glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, width, height, 0, GL_RGB, GL_UNSIGNED_BYTE, NULL);

		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

		glBindTexture(GL_TEXTURE_2D, 0);

		
		glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0,GL_TEXTURE_2D, textureAttach, 0);
		
		return textureAttach;
	}

	unsigned int createDepthTextureAttachment(int width, int height) {
		unsigned int texture;
		glGenTextures(1, &texture);
		glBindTexture(GL_TEXTURE_2D, texture);
		glTexImage2D(GL_TEXTURE_2D, 0, GL_DEPTH24_STENCIL8, width, height, 0,GL_DEPTH_STENCIL, GL_UNSIGNED_INT_24_8, NULL);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
		
		glFramebufferTexture2D(GL_FRAMEBUFFER, GL_DEPTH_STENCIL_ATTACHMENT, GL_TEXTURE_2D, texture, 0);
		return texture;
	}

	unsigned int createDepthBufferAttachment(int width, int height) {
		unsigned int depthBuffer;
		glGenRenderbuffers(1, &depthBuffer);
		glBindRenderbuffer(GL_RENDERBUFFER, depthBuffer);
		glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH_COMPONENT, width,height);
		glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT,GL_RENDERBUFFER, depthBuffer);
		return depthBuffer;
	}
	

	void framebuffer_resize_callback(GLFWwindow* window, int width, int height)
	{
		glViewport(0, 0, width, height);
	}

	void mouse_callback(GLFWwindow* window, double xpos, double ypos)
	{
		if (firstMouse)
		{
			lastX = xpos;
			lastY = ypos;
			firstMouse = false;
		}

		float xoffset = xpos - lastX;
		float yoffset = lastY - ypos; // reversed since y-coordinates go from bottom to top

		lastX = xpos;
		lastY = ypos;

		camera.ProcessMouseMovement(xoffset, yoffset);
	}

	// glfw: whenever the mouse scroll wheel scrolls, this callback is called
	// ----------------------------------------------------------------------
	void scroll_callback(GLFWwindow* window, double xoffset, double yoffset)
	{
		camera.ProcessMouseScroll(yoffset);
	}

	
	void mainLoop() {

		while (!glfwWindowShouldClose(window))
		{

			float currentFrame = glfwGetTime();
			deltaTime = currentFrame - lastFrame;
			grid->clk.deltaTime = deltaTime;
			lastFrame = currentFrame;
			//cout <<"deltaTime: " << deltaTime << endl;
			//Input
			processInput();

			//cout << grid->isSimulating << endl;
			//simulate
			if (startSim) {
				grid->simulate();
				//grid->simTime += deltaTime;
			}

			if (startCollSim) {

				for (int i = 0; i < collObjects.size(); i++) {
					//collisionSphere.simulate();
					collObjects[i].simulate();
				}


			}
			
			//Render
			render();

			glfwSwapBuffers(window);
			glfwPollEvents();
		}
	}

	void render() {

		//cout << camera.Position.x << ", "<< camera.Position.y << ", "<< camera.Position.z << ", " << camera.Yaw << ", " << camera.Pitch<< ", " << camera.Zoom<<endl;
		// move light position over time
		
		//lightPos.y -= sin(glfwGetTime() * 0.5)*0.1;


		float near_plane;
		float far_plane;

		glm::mat4 projection = glm::mat4(1.0f);
		glm::mat4 view = glm::mat4(1.0f);

		if (useSSR == false && useDirectRaytracing == false) {
			glEnable(GL_BLEND);
			//render scene to cubemap
			near_plane = 0.1f;
			far_plane = 100.f;

			glViewport(0, 0, cubemapSize, cubemapSize);


			//projection = glm::perspective(glm::radians(camera.Zoom), 1.0f, near_plane, far_plane);
			glm::mat4 currentCubeMapView = glm::mat4(0.0);

			GLenum cubemapSides[6] = {
				GL_TEXTURE_CUBE_MAP_POSITIVE_X,
				GL_TEXTURE_CUBE_MAP_NEGATIVE_X,
				GL_TEXTURE_CUBE_MAP_POSITIVE_Y,
				GL_TEXTURE_CUBE_MAP_NEGATIVE_Y,
				GL_TEXTURE_CUBE_MAP_POSITIVE_Z,
				GL_TEXTURE_CUBE_MAP_NEGATIVE_Z,
			};

			glm::vec3 targetVectors[6] = {
				glm::vec3(1.0f, 0.0f, 0.0f),
				glm::vec3(-1.0f, 0.0f, 0.0f),
				glm::vec3(0.0f, 1.0f, 0.0f),
				glm::vec3(0.0f, -1.0f, 0.0f),
				glm::vec3(0.0f, 0.0f, 1.0f),
				glm::vec3(0.0f, 0.0f, -1.0f)
			};

			glm::vec3 upVectors[6] = {
				glm::vec3(0.0f, -1.0f, 0.0f),
				glm::vec3(0.0f, -1.0f, 0.0f),
				glm::vec3(0.0f, 0.0f, 1.0f),
				glm::vec3(0.0f, 0.0f, -1.0f),
				glm::vec3(0.0f, -1.0f, 0.0f),
				glm::vec3(0.0f, -1.0f, 0.0f)
			};

			projection = glm::mat4(1.0);
			projection = glm::perspective(glm::radians(90.0f), 1.0f, near_plane, far_plane);

			glm::vec3 cm = collObjects[2].getPos();// pCM;
			//projection = glm::perspective(glm::radians(camera.Zoom), (float)1000 / (float)800, 0.1f, 1000.0f);
			glBindFramebuffer(GL_FRAMEBUFFER, framebuffer);
			glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
			glm::vec3 move = glm::vec3(0, -2.5, -grid->vars.kZMax);
			for (GLuint face = 0; face < 6; ++face)
			{


				glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, cubemapSides[face], dcmap, 0);
				glClear(GL_DEPTH_BUFFER_BIT);


				glm::vec3 modelOrigin = glm::vec3(0, 0, 0) + cm + move;

				float camDist = glm::length(modelOrigin - camera.Position);
				camDist1 = camDist;

				glm::vec3 object = collObjects[0].getPos() + move;

				glm::vec3 objectToWaterSurface = modelOrigin - object;

				glm::vec3 cubemapDir = glm::normalize(objectToWaterSurface);



				float cubemap1Length = objectToWaterSurface.length();

				cubemapCenter1 = object + 0.5f*cubemap1Length*cubemapDir;
				cubemapCenter1 = glm::vec3(0, 0, 0);//object;

				bboxmin1 = cubemapCenter1 - cubemap1Length*0.5f;
				bboxmax1 = cubemapCenter1 + cubemap1Length*0.5f;

				currentCubeMapView = glm::mat4(0.0);
				currentCubeMapView = glm::lookAt(modelOrigin, modelOrigin + targetVectors[face], upVectors[face]);

				// Render skybox

				glDepthFunc(GL_LEQUAL);  // change depth function so depth test passes when values are equal to depth buffer's content
				reflSkyboxShader.use();
				view = glm::mat4(glm::mat3(camera.GetViewMatrix())); // remove translation from the view matrix
				//currentCubeMapView = glm::mat4(glm::mat3(currentCubeMapView));
				//currentCubeMapView2 = glm::mat4(glm::mat3(currentCubeMapView2));
				reflSkyboxShader.setMat4("view", glm::mat4(glm::mat3(currentCubeMapView)));
				//reflSkyboxShader.setMat4("view", currentCubeMapView);
				reflSkyboxShader.setMat4("projection", projection);
				// skybox cube
				glBindVertexArray(skyboxVAO);
				glActiveTexture(GL_TEXTURE0);
				glBindTexture(GL_TEXTURE_CUBE_MAP, cubemapTexture);
				glDrawArrays(GL_TRIANGLES, 0, 36);
				glBindVertexArray(0);
				glDepthFunc(GL_LESS); // set depth function back to default



				shader.use();
				shader.setInt("diffuseTexture", 0);
				shader.setInt("depthMap", 1);
				shader.setInt("skybox", 2);
				shader.setInt("normalMap", 3);
				//shader.setInt("depth2D",4);
				shader.setVec3("lightPos", lightPos);
				shader.setMat4("view", currentCubeMapView);
				shader.setMat4("projection", projection);
				shader.setFloat("far_plane", 25.0f);
				shader.setVec3("viewPosition", camera.Position);
				shader.setVec3("camPos", camera.Position);



				renderCollisionObjectAt(shader, 0, camera, true);
				renderCollisionObjectAt(shader, 1, camera, true);




			}

			projection = glm::mat4(1.0);
			projection = glm::perspective(glm::radians(90.0f), 1.0f, near_plane, far_plane);
			glBindFramebuffer(GL_FRAMEBUFFER, 0); // Unbind FBO, set default framebuffer

			glViewport(0, 0, cubemapSize, cubemapSize);
			glBindFramebuffer(GL_FRAMEBUFFER, noSphereFrameBuffer);
			glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
			//projection = camera.perspectiveProjection(90, 1, near_plane, far_plane);
			for (GLuint face = 0; face < 6; ++face)
			{


				glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, cubemapSides[face], noSphereMap, 0);
				glClear(GL_DEPTH_BUFFER_BIT);


				glm::vec3 modelOrigin = glm::vec3(0, 0, 0) + collObjects[2].getPos() + move;

				float camDist = glm::length(modelOrigin - camera.Position);
				camDist2 = camDist;

				glm::vec3 objectToMirror = modelOrigin - collObjects[0].getPos() - move;

				float cubemap2Length = objectToMirror.length();

				cubemapCenter2 = collObjects[0].getPos() + move + 0.5f*cubemap2Length*glm::normalize(objectToMirror);

				bboxmin2 = cubemapCenter2 - cubemap2Length*5.5f;
				bboxmax2 = cubemapCenter2 + cubemap2Length*5.5f;



				currentCubeMapView = glm::mat4(1.0);
				currentCubeMapView = glm::lookAt(modelOrigin, modelOrigin + targetVectors[face], upVectors[face]);

				// Render skybox
				glDepthFunc(GL_LEQUAL);  // change depth function so depth test passes when values are equal to depth buffer's content
				reflSkyboxShader.use();
				view = glm::mat4(glm::mat3(camera.GetViewMatrix())); // remove translation from the view matrix

				reflSkyboxShader.setMat4("view", glm::mat4(glm::mat3(currentCubeMapView)));
				reflSkyboxShader.setMat4("projection", projection);

				// skybox cube
				glBindVertexArray(skyboxVAO);
				glActiveTexture(GL_TEXTURE0);
				glBindTexture(GL_TEXTURE_CUBE_MAP, cubemapTexture);
				glDrawArrays(GL_TRIANGLES, 0, 36);
				glBindVertexArray(0);
				glDepthFunc(GL_LESS); // set depth function back to default


				shader.use();
				shader.setInt("diffuseTexture", 0);
				shader.setInt("depthMap", 1);
				shader.setInt("skybox", 2);
				shader.setInt("normalMap", 3);
				//shader.setInt("depth2D", 4);
				shader.setVec3("lightPos", lightPos);
				shader.setMat4("view", currentCubeMapView);
				shader.setMat4("projection", projection);
				shader.setVec3("viewPosition", camera.Position);
				shader.setVec3("camPos", camera.Position);
				shader.setFloat("far_plane", 25.0f);

				renderCollisionObjectAt(shader, 0, camera, true);
				renderCollisionObjectAt(shader, 1, camera, true);

				//shader.use();




				glActiveTexture(GL_TEXTURE0);
				glBindTexture(GL_TEXTURE_2D, planeTexture);

				glActiveTexture(GL_TEXTURE3);
				glBindTexture(GL_TEXTURE_2D, planeNormalMap);

				shader.setVec4("colorVar", glm::vec4(0, 0, 0, 1.0));
				shader.setFloat("water", 0.0);
				shader.setBool("isParticle", false);
				shader.setInt("shadows", shadows);


				glm::mat4 model = glm::mat4(1.0);
				model = glm::scale(model, glm::vec3(5.0f));
				shader.setMat4("model", model);
				//glDisable(GL_CULL_FACE); // note that we disable culling here since we render 'inside' the cube instead of the usual 'outside' which throws off the normal culling methods.
										 //shader.setInt("reverse_normals", 1); // A small little hack to invert normals when drawing cube from the inside so lighting still works.
										 //renderCube();
										 //renderPlane();
				shader.setInt("reverse_normals", 0);

				renderQuad();


				glActiveTexture(GL_TEXTURE0);
				glBindTexture(GL_TEXTURE_2D, 0);

				glActiveTexture(GL_TEXTURE3);
				glBindTexture(GL_TEXTURE_2D, 0);

			}

			glBindFramebuffer(GL_FRAMEBUFFER, 0);
		}

		glViewport(0, 0, 1000, 800);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
		// 0. create depth cubemap transformation matrices
		// -----------------------------------------------
		near_plane = 1.0f;
		far_plane = 25.0f;
		glm::mat4 shadowProj = glm::perspective(glm::radians(90.0f), (float)SHADOW_WIDTH / (float)SHADOW_HEIGHT, near_plane, far_plane);
		std::vector<glm::mat4> shadowTransforms;
		shadowTransforms.push_back(shadowProj * glm::lookAt(lightPos, lightPos + glm::vec3(1.0f, 0.0f, 0.0f), glm::vec3(0.0f, -1.0f, 0.0f)));
		shadowTransforms.push_back(shadowProj * glm::lookAt(lightPos, lightPos + glm::vec3(-1.0f, 0.0f, 0.0f), glm::vec3(0.0f, -1.0f, 0.0f)));
		shadowTransforms.push_back(shadowProj * glm::lookAt(lightPos, lightPos + glm::vec3(0.0f, 1.0f, 0.0f), glm::vec3(0.0f, 0.0f, 1.0f)));
		shadowTransforms.push_back(shadowProj * glm::lookAt(lightPos, lightPos + glm::vec3(0.0f, -1.0f, 0.0f), glm::vec3(0.0f, 0.0f, -1.0f)));
		shadowTransforms.push_back(shadowProj * glm::lookAt(lightPos, lightPos + glm::vec3(0.0f, 0.0f, 1.0f), glm::vec3(0.0f, -1.0f, 0.0f)));
		shadowTransforms.push_back(shadowProj * glm::lookAt(lightPos, lightPos + glm::vec3(0.0f, 0.0f, -1.0f), glm::vec3(0.0f, -1.0f, 0.0f)));

		// 1. render scene to depth cubemap
		// --------------------------------
		glViewport(0, 0, SHADOW_WIDTH, SHADOW_HEIGHT);
		glBindFramebuffer(GL_FRAMEBUFFER, depthMapFBO);
		glClear(GL_DEPTH_BUFFER_BIT);
		simpleDepthShader.use();
		for (unsigned int i = 0; i < 6; ++i)
			simpleDepthShader.setMat4("shadowMatrices[" + std::to_string(i) + "]", shadowTransforms[i]);
		simpleDepthShader.setFloat("far_plane", far_plane);
		simpleDepthShader.setVec3("lightPos", lightPos);
		renderScene(simpleDepthShader);
		glBindFramebuffer(GL_FRAMEBUFFER, 0);


		
		glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
		glDisable(GL_BLEND);
		//Deferred Rendering - Render offscreen variables to gBuffer
		// ------------------------
		glBindFramebuffer(GL_FRAMEBUFFER, gBuffer);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
		projection = glm::perspective(glm::radians(camera.Zoom), (float)1000 / (float)800, 0.1f, 1000.0f);
		view = camera.GetViewMatrix();
		glm::mat4 model = glm::mat4(1.0f);
		gBufferShader.use();
		gBufferShader.setMat4("projection", projection);
		gBufferShader.setMat4("view", view);
		gBufferShader.setVec3("camPos", camera.Position);
		executeMC = true;
		drawSphere = false;

		glActiveTexture(GL_TEXTURE1);
		glBindTexture(GL_TEXTURE_CUBE_MAP, depthCubemap);
		//glActiveTexture(GL_TEXTURE2);
		
		if(useSSR)renderScene(gBufferShader);
		//executeMC = false;
		glBindFramebuffer(GL_FRAMEBUFFER, 0);
		//glEnable(GL_BLEND);

		// generate ssao texture
		// ------------------------
		glBindFramebuffer(GL_FRAMEBUFFER, ssaoFBO);
		glClear(GL_COLOR_BUFFER_BIT);
		ssaoShader.use();
		// Send kernel + rotation 
		for (unsigned int i = 0; i < 64; ++i)
			ssaoShader.setVec3("samples[" + std::to_string(i) + "]", ssaoKernel[i]);
		ssaoShader.setMat4("projection", projection);
		glActiveTexture(GL_TEXTURE0);
		glBindTexture(GL_TEXTURE_2D, gPosition);
		glActiveTexture(GL_TEXTURE1);
		glBindTexture(GL_TEXTURE_2D, gNormal);
		glActiveTexture(GL_TEXTURE2);
		glBindTexture(GL_TEXTURE_2D, noiseTexture);
		if (useSSR)renderScreenQuad();
		glBindFramebuffer(GL_FRAMEBUFFER, 0);

		// 3. blur SSAO texture to remove noise
		// ------------------------------------
		glBindFramebuffer(GL_FRAMEBUFFER, ssaoBlurFBO);
		glClear(GL_COLOR_BUFFER_BIT);
		ssaoBlurShader.use();
		glActiveTexture(GL_TEXTURE0);
		glBindTexture(GL_TEXTURE_2D, ssaoColorBuffer);
		if (useSSR)renderScreenQuad();
		glBindFramebuffer(GL_FRAMEBUFFER, 0);


		glBindFramebuffer(GL_FRAMEBUFFER, prevFrameBuffer);
		glEnable(GL_DEPTH_TEST);
		glEnable(GL_BLEND);
		// 2. render scene as normal 
		// -------------------------
		glViewport(0, 0, 1000, 800);
		glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
		shader.use();
		projection = glm::perspective(glm::radians(camera.Zoom), (float)1000 / (float)800, 0.1f, 1000.0f);
		view = camera.GetViewMatrix();
		shader.setMat4("projection", projection);
		shader.setMat4("view", view);
		// set lighting uniforms
		shader.setVec3("lightPos", lightPos);
		if (useSSR)shader.setBool("isSSR", true);
		else shader.setBool("isSSR", false);
		shader.setVec3("viewPosition", camera.Position);
		shader.setVec3("camPos", camera.Position);
		shader.setInt("shadows", shadows); // enable/disable shadows by pressing 'SPACE'
		shader.setFloat("far_plane", far_plane);
		glActiveTexture(GL_TEXTURE0);
		glBindTexture(GL_TEXTURE_2D, planeTexture);
		glActiveTexture(GL_TEXTURE1);
		glBindTexture(GL_TEXTURE_CUBE_MAP, depthCubemap);
		glActiveTexture(GL_TEXTURE2);
		if(useDirectRaytracing == false && useSSR == false)glBindTexture(GL_TEXTURE_CUBE_MAP, dcmap);
		else glBindTexture(GL_TEXTURE_CUBE_MAP, cubemapTexture);
		
		glActiveTexture(GL_TEXTURE3);
		glBindTexture(GL_TEXTURE_2D, planeNormalMap);

		glActiveTexture(GL_TEXTURE4);
		glBindTexture(GL_TEXTURE_2D, gPosition);

		glActiveTexture(GL_TEXTURE5);
		glBindTexture(GL_TEXTURE_2D, gNormal);

		glActiveTexture(GL_TEXTURE6);
		glBindTexture(GL_TEXTURE_2D, gAlbedoSpec);

		glActiveTexture(GL_TEXTURE7);
		glBindTexture(GL_TEXTURE_2D, gColor);

		glActiveTexture(GL_TEXTURE8);
		glBindTexture(GL_TEXTURE_2D, depthTexture);

		glActiveTexture(GL_TEXTURE9);
		glBindTexture(GL_TEXTURE_2D, ssaoColorBufferBlur);

		shader.setVec4("colorVar", glm::vec4(0, 0, 0, 1.0));
		shader.setFloat("water", 0.0);
		shader.setBool("isParticle", false);

		glm::mat4 cameraClip = camera.getCameraToClipMatrix(90.0f, near_plane, far_plane);
		shader.setMat4("cameraToClipMatrix", cameraClip);
		shader.setFloat("angleInput", scaleAway);
		clippingPlane = Vector4(0, -1, 0, 15);
		//glDisable(GL_CLIP_DISTANCE0);
		shader.setVec4("plane", glm::vec4(clippingPlane.x, clippingPlane.y, clippingPlane.z, clippingPlane.w));
		if(useDirectRaytracing==false)renderScene(shader);



		glEnable(GL_DEPTH_TEST);

		//render skybox
		//glDepthMask(GL_FALSE);
		glDepthFunc(GL_LEQUAL);  // change depth function so depth test passes when values are equal to depth buffer's content
		skyboxShader.use();
		view = glm::mat4(glm::mat3(camera.GetViewMatrix())); // remove translation from the view matrix
		skyboxShader.setMat4("view", view);
		skyboxShader.setMat4("projection", projection);
		// skybox cube
		glBindVertexArray(skyboxVAO);
		glActiveTexture(GL_TEXTURE0);
		glBindTexture(GL_TEXTURE_CUBE_MAP, cubemapTexture);
		glDrawArrays(GL_TRIANGLES, 0, 36);
		glBindVertexArray(0);
		glDepthFunc(GL_LESS); // set depth function back to default
		//glDepthMask(GL_TRUE);
		//cout << camera.Position.x << ", " << camera.Position.y << ", " <<camera.Position.z << ", " << endl;
		glBindFramebuffer(GL_FRAMEBUFFER, 0);

		if (useSSR) {
			glBindFramebuffer(GL_READ_FRAMEBUFFER, prevFrameBuffer);
			glBindFramebuffer(GL_DRAW_FRAMEBUFFER, copyBuffer); // write to default framebuffer
													   // blit to default framebuffer. Note that this may or may not work as the internal formats of both the FBO and default framebuffer have to match.
													   // the internal formats are implementation defined. This works on all of my systems, but if it doesn't on yours you'll likely have to write to the
													   // depth buffer in another shader stage (or somehow see to match the default framebuffer's internal format with the FBO's internal format).
			glBlitFramebuffer(0, 0, 1000, 800, 0, 0, 1000, 800, GL_COLOR_BUFFER_BIT, GL_NEAREST);
			glBindFramebuffer(GL_FRAMEBUFFER, 0);
		}
		
		if (useDirectRaytracing == true) {
			
			computeShaderTrace(projection*camera.GetViewMatrix());
		}

		
		//debug g-buffer
		screenShader.use();
		glBindVertexArray(screenQuadVAO);
		glDisable(GL_DEPTH_TEST);
		glActiveTexture(GL_TEXTURE0);
		if (useDirectRaytracing==true){
			
			glBindTexture(GL_TEXTURE_2D, computeShaderTexture);
		}
		else glBindTexture(GL_TEXTURE_2D, texColorBuffer);
		//glDisable(GL_DEPTH_TEST);
		glDrawArrays(GL_TRIANGLES, 0, 6);
		glEnable(GL_DEPTH_TEST);
	}

	int nextPowerOfTwo(int x) {
		x--;
		x |= x >> 1; // handle 2 bit numbers
		x |= x >> 2; // handle 4 bit numbers
		x |= x >> 4; // handle 8 bit numbers
		x |= x >> 8; // handle 16 bit numbers
		x |= x >> 16; // handle 32 bit numbers
		x++;
		return x;
	}

	void computeShaderTrace(glm::mat4 projectionview) {
		glm::mat4 invProjectionView = glm::inverse(projectionview);

		float lightLength = 25;
		glm::mat4 lightViewMatrix = glm::lookAt(lightPos, lightPos + glm::vec3(0,-1,0),glm::vec3(0,1,0));
		lightProjection = glm::ortho(lightPos.x- lightLength, lightPos.x+ lightLength, lightPos.y- lightLength, lightPos.y+ lightLength, lightPos.z+ lightLength,lightPos.z- lightLength)*lightViewMatrix;

		computeShader.use();

		computeShader.setVec3("eye", camera.Position);
		computeShader.setVec3("ray00", camera.GetEyeRay(-1,-1, invProjectionView));
		computeShader.setVec3("ray01", camera.GetEyeRay(-1, 1, invProjectionView));
		computeShader.setVec3("ray10", camera.GetEyeRay( 1, -1, invProjectionView));
		computeShader.setVec3("ray11", camera.GetEyeRay( 1,  1, invProjectionView));

		glm::vec3 move = glm::vec3(0,-2.5,-30);
		
		computeShader.setVec4("sphere", glm::vec4(collObjects[0].getPos(),collObjects[0].radius));
		computeShader.setVec4("cube", glm::vec4(collObjects[1].getPos(), collObjects[1].cubeSize));
		computeShader.setVec4("mirror", glm::vec4(collObjects[2].getPos(), 10));
		computeShader.setVec4("customModel", glm::vec4(customModelPos, 1));
		computeShader.setVec4("lightPos", glm::vec4(lightPos, 1));
		computeShader.setMat4("lightProjection", lightProjection);

		computeShader.setFloat("fBias", fBias);
		computeShader.setFloat("fScale", fScale);
		computeShader.setFloat("fPower", fPower);

		computeShader.setBool("isIrridescent", isIrridescent);
		computeShader.setFloat("gScattering", gScattering);
		computeShader.setFloat("scatteringAmount", scatteringAmount);

		
		//add triangles of sphere
		computeShader.setInt("sphereTriangleSize", collObjects[0].ssbo_triangles.size());
		
		//updateSSBOData(ssbo_sphere, mcCubes.triangulator.mesh.ssbo_triangles, 1, "SphereTriangles");
		//renderSSBOData(ssbo_sphere, "sphereTriangleSize", mcCubes.triangulator.mesh.ssbo_triangles);
		
		


		
		glBindImageTexture(0, computeShaderTexture,0, GL_FALSE, 0, GL_WRITE_ONLY, GL_RGBA32F);

		/* Compute appropriate invocation dimension. */
		unsigned int worksizeX = nextPowerOfTwo(1000);
		unsigned int worksizeY = nextPowerOfTwo(800);
		
		glActiveTexture(GL_TEXTURE0);
		glBindTexture(GL_TEXTURE_CUBE_MAP, cubemapTexture);
		glActiveTexture(GL_TEXTURE1);
		glBindTexture(GL_TEXTURE_2D, planeTexture);
		glActiveTexture(GL_TEXTURE2);
		glBindTexture(GL_TEXTURE_2D, ssaoColorBufferBlur);
		glActiveTexture(GL_TEXTURE3);
		glBindTexture(GL_TEXTURE_2D, depthTexture);
		//glBindTexture(GL_TEXTURE_CUBE_MAP, depthCubemap);
		

		/* Invoke the compute shader. */
		//glDispatchCompute(worksizeX / computeShader.workGroupSizeX, worksizeY / computeShader.workGroupSizeY, 1);
		//glDispatchCompute(worksizeX, worksizeY, 1);
		glDispatchCompute(worksizeX/32, worksizeY/32, 1);
		//glDispatchCompute(worksizeX / 32, worksizeY / 32, 1);
		
		//glDispatchCompute(32, 32, 1);
		//glDispatchCompute(computeShader.work_grp_count[0], computeShader.work_grp_count[1], 1);

		/* Reset image binding. */
		glBindImageTexture(0, 0, 0, GL_FALSE, 0, GL_READ_WRITE, GL_RGBA32F);
		glMemoryBarrier(GL_SHADER_IMAGE_ACCESS_BARRIER_BIT);

		

		//glMemoryBarrier(GL_ALL_BARRIER_BITS);
		//glUseProgram(0);
		//unBindSSBO();
	}

	// renders the 3D scene
	// --------------------
	void renderScene(const Shader &shader, bool refl=false)
	{

		

		glm::mat4 model = glm::mat4(1.0f);
		if (refl == false) {
			// room cube

			model = glm::scale(model, glm::vec3(5.0f));
			shader.setMat4("model", model);
			glDisable(GL_CULL_FACE); // note that we disable culling here since we render 'inside' the cube instead of the usual 'outside' which throws off the normal culling methods.
			//shader.setInt("reverse_normals", 1); // A small little hack to invert normals when drawing cube from the inside so lighting still works.
			//renderCube();
			//renderPlane();
			shader.setInt("reverse_normals",0);

			renderQuad();
			//shader.setInt("reverse_normals", 0); // and of course disable it
			//glEnable(GL_CULL_FACE);
			// cubes
			model = glm::mat4(1.0f);
			model = glm::translate(model, glm::vec3(4.0f, -3.5f, 0.0));
			model = glm::scale(model, glm::vec3(0.5f));
			shader.setMat4("model", model);
			//renderCube();
			model = glm::mat4(1.0f);
			model = glm::translate(model, glm::vec3(2.0f, 3.0f, 1.0));
			model = glm::scale(model, glm::vec3(0.75f));
			shader.setMat4("model", model);
			//renderCube();
			model = glm::mat4(1.0f);
			model = glm::translate(model, glm::vec3(-3.0f, -1.0f, 0.0));
			model = glm::scale(model, glm::vec3(0.5f));
			shader.setMat4("model", model);
			//renderCube();
			model = glm::mat4(1.0f);
			model = glm::translate(model, glm::vec3(-1.5f, 1.0f, 1.5));
			model = glm::scale(model, glm::vec3(0.5f));
			shader.setMat4("model", model);
			//renderCube();
			model = glm::mat4(1.0f);
			model = glm::translate(model, glm::vec3(-1.5f, 2.0f, -3.0));
			model = glm::rotate(model, glm::radians(60.0f), glm::normalize(glm::vec3(1.0, 0.0, 1.0)));
			model = glm::scale(model, glm::vec3(0.75f));
			shader.setMat4("model", model);
			//renderCube();

			model = glm::mat4(1.0f);
			/*model = glm::translate(model, glm::vec3(-1.5f, 2.0f, -3.0));
			model = glm::rotate(model, glm::radians(60.0f), glm::normalize(glm::vec3(1.0, 0.0, 1.0)));
			model = glm::scale(model, glm::vec3(0.75f));*/
			shader.setMat4("model", model);
		}

		renderCollisionObjectAt(shader, 0, camera);
		renderCollisionObjectAt(shader, 1, camera);
		
		//renderCollisionObjectAt(shader, 2, camera);
		//renderCollisionObjects(shader, refl);
		
		if (refl == false) {
			//glActiveTexture(GL_TEXTURE0);
			//glBindTexture(GL_TEXTURE_2D, waterTexture);
			//renderParticles(shader);
			//renderParticlesInstanced(shader);
			
			glEnable(GL_CULL_FACE);
			if(renderWater)if (shader.ID != 10)renderSurface(shader);
			glDisable(GL_CULL_FACE);
		}
		shader.setBool("isParticle", false);

		
	}

	void renderCollisionObjectAt(const Shader &shader, int index, Camera camera_, bool refl=false) {
		

		//if (index == 0 && drawSphere==false)return;

		if (refl==false) {
			//glActiveTexture(GL_TEXTURE2);
			//if(index!=2)glBindTexture(GL_TEXTURE_CUBE_MAP, noSphereMap);
			//else glBindTexture(GL_TEXTURE_CUBE_MAP, noSphereMap);

			//shader.setMat4("camView", camera.GetViewMatrix());


			

		}
		else {
			//glActiveTexture(GL_TEXTURE2);
			//glBindTexture(GL_TEXTURE_CUBE_MAP, cubemapTexture);
		}
		if (useDirectRaytracing == false && useSSR == false) {
			glActiveTexture(GL_TEXTURE2);
			glBindTexture(GL_TEXTURE_CUBE_MAP, noSphereMap);
		}
		glEnable(GL_CULL_FACE);
		
		Vector3 move = Vector3(0, -2.5, -30);// grid->vars.kZMax);

		glm::mat4 model;
		model = glm::mat4(1.0);

		

		//else model_ = glm::mat4(1.0f);

		float scaleX = grid->vars.defaultSize / grid->vars.kXMax;
		float scaleY = grid->vars.defaultSize / grid->vars.kYMax;
		float scaleZ = grid->vars.defaultSize / grid->vars.kZMax;
		//model = glm::scale(model, glm::vec3(scaleX, scaleY, scaleZ));

		model = glm::translate(model, glm::vec3(collObjects[index].pos.x + move.x, collObjects[index].pos.y + move.y, collObjects[index].pos.z + move.z));
		shader.setVec4("colorVar", glm::vec4(collObjects[index].color.x, collObjects[index].color.y, collObjects[index].color.z, 0.3));
		shader.setFloat("water", collObjects[index].id);
		shader.setVec3("lightPos", lightPos);

		shader.setVec3("viewPosition", camera.Position);
		shader.setVec3("camPos", camera.Position);
		
		
		//shader.setBool("isParticle", false);
		shader.setInt("diffuseTexture", 0);
		shader.setVec3("sphereCenter", glm::vec3(collObjects[index].pos.x + move.x, collObjects[index].pos.y + move.y, collObjects[index].pos.z + move.z));
		//shader.setInt("reverse_normals", 1);
		//if (refl)model = glm::scale(model, glm::vec3(0.5));
		

		shader.setMat4("model", model);
		//else shader.setMat4("model", model_);
		collObjects[index].texture = planeTexture;

		//collisionSphere.render();
		collObjects[index].render();
		glDisable(GL_CULL_FACE);

	}

	void renderCollisionObjects(const Shader &shader, bool refl=false) {
		
		
		
		

		glEnable(GL_CULL_FACE);
		for (int i = 0; i < collObjects.size(); i++) {

			if (i == 0) {
				glActiveTexture(GL_TEXTURE2);
				glBindTexture(GL_TEXTURE_CUBE_MAP, noSphereMap);
			}
			else {
				glActiveTexture(GL_TEXTURE2);
				glBindTexture(GL_TEXTURE_CUBE_MAP, 0);
			}
			Vector3 move = Vector3(0, -2.5, -30);// grid->vars.kZMax);

			glm::mat4 model;
			model = glm::mat4(1.0f);
		
			//else model_ = glm::mat4(1.0f);

			float scaleX = grid->vars.defaultSize / grid->vars.kXMax;
			float scaleY = grid->vars.defaultSize / grid->vars.kYMax;
			float scaleZ = grid->vars.defaultSize / grid->vars.kZMax;
			//model = glm::scale(model, glm::vec3(scaleX, scaleY, scaleZ));

			model = glm::translate(model, glm::vec3(collObjects[i].pos.x + move.x, collObjects[i].pos.y + move.y, collObjects[i].pos.z + move.z));
			shader.setVec4("colorVar", glm::vec4(collObjects[i].color.x, collObjects[i].color.y, collObjects[i].color.z, 0.3));
			shader.setFloat("water", collObjects[i].id);
			shader.setVec3("lightPos", lightPos);
			
			shader.setVec3("viewPosition", camera.Position);
			
			shader.setVec3("camPos", camera.Position);
			shader.setInt("diffuseTexture", 0);
			shader.setVec3("sphereCenter",glm::vec3(collObjects[i].pos.x+move.x, collObjects[i].pos.y+move.y, collObjects[i].pos.z+move.z));
			//shader.setInt("reverse_normals", 1);
			//if (refl)model = glm::scale(model, glm::vec3(0.5));
			shader.setMat4("model", model);
			//else shader.setMat4("model", model_);
			//collObjects[i].texture = planeTexture;

			//collisionSphere.render();
			collObjects[i].render();
		}
		glDisable(GL_CULL_FACE);
	}

	void renderSurface(const Shader &shader) {
		
		//glActiveTexture(GL_TEXTURE2);
		//glBindTexture(GL_TEXTURE_CUBE_MAP, dcmap);

		//mcCubes.generateGeometry();
		//float time0 = glfwGetTime();
		//executeMC = true;
		if (executeMC || (useDirectRaytracing==false && useSSR == false)) {

			if (mcCubes.runOpenCL == 0) {
				//CPU implementation
				mcCubes.generateSurface();
			}
			else {
				//GPU implementation
				mcCubes.generateSurfaceCL();
			}
		}
		//float time1 = glfwGetTime();

		//cout << "Time to complete: " << (time1-time0) << endl;

		shader.setVec3("cubemapCenter", cubemapCenter1);
		shader.setVec3("bboxmin", bboxmin1);
		shader.setVec3("bboxmax", bboxmax1);
		shader.setVec3("camPos", camera.Position);


		//if(mcCubes.triangulator.mesh.m_vertices.size()!=0)cout << mcCubes.triangulator.mesh.m_vertices.size() << endl;
		//cout << mcCubes.triangulator.mesh.m_vertices.size() << endl;
		//mcCubes.generateGeometry();
		//glm::vec3 cm = mcCubes.getCenterMass();
		Vector3 move = Vector3(0, -2.5, -grid->vars.kZMax);
		glm::mat4 model = glm::mat4(1.0f);

		float scaleX = grid->vars.defaultSize / grid->vars.kXMax;
		float scaleY = grid->vars.defaultSize / grid->vars.kYMax;
		float scaleZ = grid->vars.defaultSize / grid->vars.kZMax;
		//model = glm::scale(model, glm::vec3(scaleX,scaleY,scaleZ));
		
		model = glm::translate(model, glm::vec3(move.x, move.y, move.z));
		shader.setVec3("lightPos", lightPos);
		shader.setVec3("camPos", camera.Position);
		shader.setVec3("viewPosition", camera.Position);
		shader.setMat4("model", model);
		
		shader.setVec4("colorVar", glm::vec4(64/255., 164/255., 223/255., 0.8));
		//shader.setVec4("colorVar", glm::vec4(1, 1, 1, 0.9));
		shader.setFloat("water",1.0);
		
		//mcCubes.renderNormalPointLines();
		//mcCubes.texture = waterTexture;
		//mcCubes.renderGridCells();
		//mcCubes.generateGeometry();
		//mcCubes.sortTriangles(Vector3(camera.Position.x, camera.Position.y, camera.Position.z));
		
		//mcCubes.diffuseTexture = planeTexture;
		mcCubes.texture = cubemapTexture;
		mcCubes.render();
		
		//pCM = mcCubes.getCenterMass();
	}


	// renderCube() renders a 1x1 3D cube in NDC.
	// -------------------------------------------------
	unsigned int cubeVAO = 0;
	unsigned int cubeVBO = 0;
	void renderCube()
	{
		// initialize (if necessary)
		if (cubeVAO == 0)
		{
			float vertices[] = {
				// back face
				-1.0f, -1.0f, -1.0f,  0.0f,  0.0f, -1.0f, 0.0f, 0.0f, // bottom-left
				1.0f,  1.0f, -1.0f,  0.0f,  0.0f, -1.0f, 1.0f, 1.0f, // top-right
				1.0f, -1.0f, -1.0f,  0.0f,  0.0f, -1.0f, 1.0f, 0.0f, // bottom-right         
				1.0f,  1.0f, -1.0f,  0.0f,  0.0f, -1.0f, 1.0f, 1.0f, // top-right
				-1.0f, -1.0f, -1.0f,  0.0f,  0.0f, -1.0f, 0.0f, 0.0f, // bottom-left
				-1.0f,  1.0f, -1.0f,  0.0f,  0.0f, -1.0f, 0.0f, 1.0f, // top-left
				// front face
				-1.0f, -1.0f,  1.0f,  0.0f,  0.0f,  1.0f, 0.0f, 0.0f, // bottom-left
				1.0f, -1.0f,  1.0f,  0.0f,  0.0f,  1.0f, 1.0f, 0.0f, // bottom-right
				1.0f,  1.0f,  1.0f,  0.0f,  0.0f,  1.0f, 1.0f, 1.0f, // top-right
				1.0f,  1.0f,  1.0f,  0.0f,  0.0f,  1.0f, 1.0f, 1.0f, // top-right
				-1.0f,  1.0f,  1.0f,  0.0f,  0.0f,  1.0f, 0.0f, 1.0f, // top-left
				-1.0f, -1.0f,  1.0f,  0.0f,  0.0f,  1.0f, 0.0f, 0.0f, // bottom-left
				// left face
				-1.0f,  1.0f,  1.0f, -1.0f,  0.0f,  0.0f, 1.0f, 0.0f, // top-right
				-1.0f,  1.0f, -1.0f, -1.0f,  0.0f,  0.0f, 1.0f, 1.0f, // top-left
				-1.0f, -1.0f, -1.0f, -1.0f,  0.0f,  0.0f, 0.0f, 1.0f, // bottom-left
				-1.0f, -1.0f, -1.0f, -1.0f,  0.0f,  0.0f, 0.0f, 1.0f, // bottom-left
				-1.0f, -1.0f,  1.0f, -1.0f,  0.0f,  0.0f, 0.0f, 0.0f, // bottom-right
				-1.0f,  1.0f,  1.0f, -1.0f,  0.0f,  0.0f, 1.0f, 0.0f, // top-right
				// right face
				1.0f,  1.0f,  1.0f,  1.0f,  0.0f,  0.0f, 1.0f, 0.0f, // top-left
				1.0f, -1.0f, -1.0f,  1.0f,  0.0f,  0.0f, 0.0f, 1.0f, // bottom-right
				1.0f,  1.0f, -1.0f,  1.0f,  0.0f,  0.0f, 1.0f, 1.0f, // top-right         
				1.0f, -1.0f, -1.0f,  1.0f,  0.0f,  0.0f, 0.0f, 1.0f, // bottom-right
				1.0f,  1.0f,  1.0f,  1.0f,  0.0f,  0.0f, 1.0f, 0.0f, // top-left
				1.0f, -1.0f,  1.0f,  1.0f,  0.0f,  0.0f, 0.0f, 0.0f, // bottom-left     
				// bottom face
				-1.0f, -1.0f, -1.0f,  0.0f, -1.0f,  0.0f, 0.0f, 1.0f, // top-right
				1.0f, -1.0f, -1.0f,  0.0f, -1.0f,  0.0f, 1.0f, 1.0f, // top-left
				1.0f, -1.0f,  1.0f,  0.0f, -1.0f,  0.0f, 1.0f, 0.0f, // bottom-left
				1.0f, -1.0f,  1.0f,  0.0f, -1.0f,  0.0f, 1.0f, 0.0f, // bottom-left
				-1.0f, -1.0f,  1.0f,  0.0f, -1.0f,  0.0f, 0.0f, 0.0f, // bottom-right
				-1.0f, -1.0f, -1.0f,  0.0f, -1.0f,  0.0f, 0.0f, 1.0f, // top-right
				// top face
				-1.0f,  1.0f, -1.0f,  0.0f,  1.0f,  0.0f, 0.0f, 1.0f, // top-left
				1.0f,  1.0f , 1.0f,  0.0f,  1.0f,  0.0f, 1.0f, 0.0f, // bottom-right
				1.0f,  1.0f, -1.0f,  0.0f,  1.0f,  0.0f, 1.0f, 1.0f, // top-right     
				1.0f,  1.0f,  1.0f,  0.0f,  1.0f,  0.0f, 1.0f, 0.0f, // bottom-right
				-1.0f,  1.0f, -1.0f,  0.0f,  1.0f,  0.0f, 0.0f, 1.0f, // top-left
				-1.0f,  1.0f,  1.0f,  0.0f,  1.0f,  0.0f, 0.0f, 0.0f  // bottom-left        
			};
			glGenVertexArrays(1, &cubeVAO);
			glGenBuffers(1, &cubeVBO);
			// fill buffer
			glBindBuffer(GL_ARRAY_BUFFER, cubeVBO);
			glBufferData(GL_ARRAY_BUFFER, sizeof(vertices), vertices, GL_STATIC_DRAW);
			// link vertex attributes
			glBindVertexArray(cubeVAO);
			glEnableVertexAttribArray(0);
			glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 8 * sizeof(float), (void*)0);
			glEnableVertexAttribArray(1);
			glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 8 * sizeof(float), (void*)(3 * sizeof(float)));
			glEnableVertexAttribArray(2);
			glVertexAttribPointer(2, 2, GL_FLOAT, GL_FALSE, 8 * sizeof(float), (void*)(6 * sizeof(float)));
			glBindBuffer(GL_ARRAY_BUFFER, 0);
			glBindVertexArray(0);
		}
		// render Cube
		glBindVertexArray(cubeVAO);
		glDrawArrays(GL_TRIANGLES, 0, 36);
		glBindVertexArray(0);
	}


	unsigned int squadVAO = 0;
	unsigned int squadVBO;
	void renderScreenQuad()
	{
		if (squadVAO == 0)
		{
			float squadVertices[] = {
				// positions        // texture Coords
				-1.0f,  1.0f, 0.0f, 0.0f, 1.0f,
				-1.0f, -1.0f, 0.0f, 0.0f, 0.0f,
				1.0f,  1.0f, 0.0f, 1.0f, 1.0f,
				1.0f, -1.0f, 0.0f, 1.0f, 0.0f,
			};
			// setup plane VAO
			glGenVertexArrays(1, &squadVAO);
			glGenBuffers(1, &squadVBO);
			glBindVertexArray(squadVAO);
			glBindBuffer(GL_ARRAY_BUFFER, squadVBO);
			glBufferData(GL_ARRAY_BUFFER, sizeof(squadVertices), &squadVertices, GL_STATIC_DRAW);
			glEnableVertexAttribArray(0);
			glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 5 * sizeof(float), (void*)0);
			glEnableVertexAttribArray(1);
			glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, 5 * sizeof(float), (void*)(3 * sizeof(float)));
		}
		glBindVertexArray(squadVAO);
		glDrawArrays(GL_TRIANGLE_STRIP, 0, 4);
		glBindVertexArray(0);
	}


	// renderQuad() renders a 1x1 XY quad in NDC
	// -----------------------------------------
	unsigned int quadVAO = 0;
	unsigned int quadVBO;
	void renderQuad()
	{
		if (quadVAO == 0)
		{

			// positions
			glm::vec3 pos1(-10.0f, -0.5f, -10.0f);
			glm::vec3 pos2(-10.0f, -0.5f, 10.0f);
			glm::vec3 pos3(10.0f, -0.5f, 10.0f);
			glm::vec3 pos4(10.0f, -0.5f, -10.0f);
			// texture coordinates
			glm::vec2 uv1(0.0f, 10.0f);
			glm::vec2 uv2(0.0f, 0.0f);
			glm::vec2 uv3(10.0f, 0.0f);
			glm::vec2 uv4(10.0f, 10.0f);
			// normal vector
			glm::vec3 nm(0.0f, 1.0f, 0.0f);

			// calculate tangent/bitangent vectors of both triangles
			glm::vec3 tangent1, bitangent1;
			glm::vec3 tangent2, bitangent2;
			// triangle 1
			// ----------
			glm::vec3 edge1 = pos2 - pos1;
			glm::vec3 edge2 = pos3 - pos1;
			glm::vec2 deltaUV1 = uv2 - uv1;
			glm::vec2 deltaUV2 = uv3 - uv1;

			GLfloat f = 1.0f / (deltaUV1.x * deltaUV2.y - deltaUV2.x * deltaUV1.y);

			tangent1.x = f * (deltaUV2.y * edge1.x - deltaUV1.y * edge2.x);
			tangent1.y = f * (deltaUV2.y * edge1.y - deltaUV1.y * edge2.y);
			tangent1.z = f * (deltaUV2.y * edge1.z - deltaUV1.y * edge2.z);
			tangent1 = glm::normalize(tangent1);

			bitangent1.x = f * (-deltaUV2.x * edge1.x + deltaUV1.x * edge2.x);
			bitangent1.y = f * (-deltaUV2.x * edge1.y + deltaUV1.x * edge2.y);
			bitangent1.z = f * (-deltaUV2.x * edge1.z + deltaUV1.x * edge2.z);
			bitangent1 = glm::normalize(bitangent1);

			// triangle 2
			// ----------
			edge1 = pos3 - pos1;
			edge2 = pos4 - pos1;
			deltaUV1 = uv3 - uv1;
			deltaUV2 = uv4 - uv1;

			f = 1.0f / (deltaUV1.x * deltaUV2.y - deltaUV2.x * deltaUV1.y);

			tangent2.x = f * (deltaUV2.y * edge1.x - deltaUV1.y * edge2.x);
			tangent2.y = f * (deltaUV2.y * edge1.y - deltaUV1.y * edge2.y);
			tangent2.z = f * (deltaUV2.y * edge1.z - deltaUV1.y * edge2.z);
			tangent2 = glm::normalize(tangent2);


			bitangent2.x = f * (-deltaUV2.x * edge1.x + deltaUV1.x * edge2.x);
			bitangent2.y = f * (-deltaUV2.x * edge1.y + deltaUV1.x * edge2.y);
			bitangent2.z = f * (-deltaUV2.x * edge1.z + deltaUV1.x * edge2.z);
			bitangent2 = glm::normalize(bitangent2);

			/*float quadVertices[] = {
				// positions            // normals         // texcoords
				10.0f, -0.5f,  10.0f,  0.0f, 1.0f, 0.0f,  10.0f,  0.0f,
				-10.0f, -0.5f,  10.0f,  0.0f, 1.0f, 0.0f,   0.0f,  0.0f,
				-10.0f, -0.5f, -10.0f,  0.0f, 1.0f, 0.0f,   0.0f, 10.0f,

				10.0f, -0.5f,  10.0f,  0.0f, 1.0f, 0.0f,  10.0f,  0.0f,
				-10.0f, -0.5f, -10.0f,  0.0f, 1.0f, 0.0f,   0.0f, 10.0f,
				10.0f, -0.5f, -10.0f,  0.0f, 1.0f, 0.0f,  10.0f, 10.0f
			};*/

			float quadVertices[] = {
				// positions            // normal         // texcoords  // tangent                          // bitangent
				pos1.x, pos1.y, pos1.z, nm.x, nm.y, nm.z, uv1.x, uv1.y, tangent1.x, tangent1.y, tangent1.z, bitangent1.x, bitangent1.y, bitangent1.z,
				pos2.x, pos2.y, pos2.z, nm.x, nm.y, nm.z, uv2.x, uv2.y, tangent1.x, tangent1.y, tangent1.z, bitangent1.x, bitangent1.y, bitangent1.z,
				pos3.x, pos3.y, pos3.z, nm.x, nm.y, nm.z, uv3.x, uv3.y, tangent1.x, tangent1.y, tangent1.z, bitangent1.x, bitangent1.y, bitangent1.z,

				pos1.x, pos1.y, pos1.z, nm.x, nm.y, nm.z, uv1.x, uv1.y, tangent2.x, tangent2.y, tangent2.z, bitangent2.x, bitangent2.y, bitangent2.z,
				pos3.x, pos3.y, pos3.z, nm.x, nm.y, nm.z, uv3.x, uv3.y, tangent2.x, tangent2.y, tangent2.z, bitangent2.x, bitangent2.y, bitangent2.z,
				pos4.x, pos4.y, pos4.z, nm.x, nm.y, nm.z, uv4.x, uv4.y, tangent2.x, tangent2.y, tangent2.z, bitangent2.x, bitangent2.y, bitangent2.z
			};

			// setup plane VAO
			glGenVertexArrays(1, &quadVAO);
			glGenBuffers(1, &quadVBO);
			glBindVertexArray(quadVAO);
			glBindBuffer(GL_ARRAY_BUFFER, quadVBO);
			glBufferData(GL_ARRAY_BUFFER, sizeof(quadVertices), &quadVertices, GL_STATIC_DRAW);
			
			
			// position attribute
			glEnableVertexAttribArray(0);
			glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 14 * sizeof(float), (void*)0);

			//normal attribute
			glEnableVertexAttribArray(1);
			glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 14 * sizeof(float), (void*)(3 * sizeof(float)));

			// texture coord attribute
			glEnableVertexAttribArray(2);
			glVertexAttribPointer(2, 2, GL_FLOAT, GL_FALSE, 14 * sizeof(float), (void*)(6 * sizeof(float)));

			glEnableVertexAttribArray(4);
			glVertexAttribPointer(4, 3, GL_FLOAT, GL_FALSE, 14 * sizeof(float), (void*)(8 * sizeof(float)));
			glEnableVertexAttribArray(5);
			glVertexAttribPointer(5, 3, GL_FLOAT, GL_FALSE, 14 * sizeof(float), (void*)(11 * sizeof(float)));
			
		
		}
		glBindVertexArray(quadVAO);
		//glDrawArrays(GL_TRIANGLE_STRIP, 0, 4);
		glDrawArrays(GL_TRIANGLES, 0, 6);
		glBindVertexArray(0);
	}

	void renderParticles(const Shader &shader) {
		int N = grid->nParticles;

		Vector3 move = Vector3(0,-2.5,-grid->vars.kZMax);
		for (int i = 0; i < N;i++) {

			//if (grid->particles[i].isRendering == false)continue;

			glm::mat4 model = glm::mat4(1.0f);
			
			model = glm::translate(model, glm::vec3(grid->particles[i].pos.x + move.x, grid->particles[i].pos.y + move.y, grid->particles[i].pos.z + move.z));
			//model = glm::translate(model, glm::vec3(move.x, move.y, move.z));
			
			//grid->particles[i].updateMesh();
			shader.setMat4("model", model);
			grid->particles[i].texture = waterTexture;
			grid->particles[i].render();

		}
	}

	void renderParticlesInstanced(const Shader &shader) {
		shader.setBool("isParticle", true);

		int N = grid->nParticles;


		Vector3 move = Vector3(0, -2.5, -grid->vars.kZMax);
		
		glm::mat4 model = glm::mat4(1.0f);
		//if(shader.ID==10)model = glm::scale(model, glm::vec3(grid->particleSize));
		//else model = glm::scale(model, glm::vec3(grid->particleSize));
		//model = glm::translate(model, glm::vec3(grid->particles[i].pos.x + move.x, grid->particles[i].pos.y + move.y, grid->particles[i].pos.z + move.z));
		shader.setMat4("model", model);
		
		for (int i = 0; i < N;i++) {


			grid->offsets[i] = move;

			grid->offsets[i] += grid->particles[i].pos;
			grid->particles[i].texture = waterTexture;
		}

		glBindVertexArray(grid->particles[0].VAO);
		glBindBuffer(GL_ARRAY_BUFFER, grid->particles[0].OVBO);
		glBufferData(GL_ARRAY_BUFFER, sizeof(Vector3) * grid->nParticles, &grid->offsets[0], GL_STREAM_DRAW);
		


		glDrawElementsInstanced(GL_TRIANGLES, grid->particles[0].indices.size(), GL_UNSIGNED_INT, 0, grid->nParticles);
		
		glBindBuffer(GL_ARRAY_BUFFER, 0);
	
	}

	void terminate() {

	

		

		glfwTerminate();

		/*glDeleteVertexArrays(1, &planeVAO);
		glDeleteVertexArrays(1, &lightVAO);
		glDeleteBuffers(1, &VBO);
		glDeleteBuffers(1, &VBO2);*/
		
		grid->terminate();
	}

	void processInput()
	{
		if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS)
			glfwSetWindowShouldClose(window, true);

		if (glfwGetKey(window, GLFW_KEY_W) == GLFW_PRESS)
			camera.ProcessKeyboard(FORWARD, deltaTime);
		if (glfwGetKey(window, GLFW_KEY_S) == GLFW_PRESS)
			camera.ProcessKeyboard(BACKWARD, deltaTime);
		if (glfwGetKey(window, GLFW_KEY_A) == GLFW_PRESS)
			camera.ProcessKeyboard(LEFT, deltaTime);
		if (glfwGetKey(window, GLFW_KEY_D) == GLFW_PRESS)
			camera.ProcessKeyboard(RIGHT, deltaTime);

		if (glfwGetKey(window, GLFW_KEY_SPACE) == GLFW_PRESS && !shadowsKeyPressed)
		{
			shadows = !shadows;
			shadowsKeyPressed = true;
		}
		if (glfwGetKey(window, GLFW_KEY_SPACE) == GLFW_RELEASE)
		{
			shadowsKeyPressed = false;
		}

		if (glfwGetKey(window, GLFW_KEY_ENTER) == GLFW_PRESS && !startSim) {
			startSim = true;
		}
		
		if (glfwGetKey(window, GLFW_KEY_Y) == GLFW_PRESS && !startCollSim) {
			startCollSim = true;
		}

		if (glfwGetKey(window, GLFW_KEY_O) == GLFW_PRESS) {
			//scaleAway++;
			gScattering += 10;
			cout << "gScattering: " << gScattering << endl;
		}

		if (glfwGetKey(window, GLFW_KEY_P) == GLFW_PRESS) {
			//scaleAway--;
			gScattering -= 10;
			cout << "gScattering: " << gScattering << endl;
		}

		if (glfwGetKey(window, GLFW_KEY_K) == GLFW_PRESS) {
			//scaleAway++;
			scatteringAmount += 10;
			cout << "scatteringAmount: " << scatteringAmount << endl;
		}

		if (glfwGetKey(window, GLFW_KEY_L) == GLFW_PRESS) {
			//scaleAway--;
			scatteringAmount -= 10;
			cout << "scatteringAmount: " << scatteringAmount << endl;
		}

		if (glfwGetKey(window, GLFW_KEY_N) == GLFW_PRESS) {
			//scaleAway++;
			fPower += 0.01;
			cout << "fPower: " << fPower << endl;
		}

		if (glfwGetKey(window, GLFW_KEY_M) == GLFW_PRESS) {
			//scaleAway--;
			fPower -= 0.01;
			cout << "fPower: " << fPower << endl;
		}

		if (glfwGetKey(window, GLFW_KEY_I) == GLFW_PRESS) {
			isIrridescent = !isIrridescent;
			cout << "isIrridescent: " << isIrridescent << endl;
		}
		
	}

	// utility function for loading a 2D texture from file
	// ---------------------------------------------------
	unsigned int loadTexture(char const * path)
	{
		
		//stbi_set_flip_vertically_on_load(1);
		unsigned int textureID;
		glGenTextures(1, &textureID);

		int width, height, nrComponents;
		unsigned char *data = stbi_load(path, &width, &height, &nrComponents, 0);

		if (data)
		{
			GLenum format;
			GLenum dataFormat;
			if (nrComponents == 1) {
				format = GL_RED;
				dataFormat = GL_RED;
			}
			else if (nrComponents == 3) {
				format = GL_SRGB;
				dataFormat = GL_RGB;
			}
			else if (nrComponents == 4) {
				format = GL_SRGB_ALPHA;
				dataFormat = GL_RGBA;
			}

			glBindTexture(GL_TEXTURE_2D, textureID);
			glTexImage2D(GL_TEXTURE_2D, 0, format, width, height, 0, dataFormat, GL_UNSIGNED_BYTE, data);
			glGenerateMipmap(GL_TEXTURE_2D);

			glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, format == GL_RGBA ? GL_CLAMP_TO_EDGE : GL_REPEAT); 
			glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, format == GL_RGBA ? GL_CLAMP_TO_EDGE : GL_REPEAT);
			glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
			glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

			stbi_image_free(data);
		}
		else
		{
			std::cout << "Texture failed to load at path: " << path << std::endl;
			stbi_image_free(data);
		}

		return textureID;
	}

	unsigned int loadCubemap(vector<std::string> faces)
	{
		unsigned int textureID;
		glGenTextures(1, &textureID);
		glBindTexture(GL_TEXTURE_CUBE_MAP, textureID);

		int width, height, nrComponents;
		for (unsigned int i = 0; i < faces.size(); i++)
		{
			unsigned char *data = stbi_load(faces[i].c_str(), &width, &height, &nrComponents, 0);
			if (data)
			{
				glTexImage2D(GL_TEXTURE_CUBE_MAP_POSITIVE_X + i, 0, GL_RGB, width, height, 0, GL_RGB, GL_UNSIGNED_BYTE, data);
				cubemapWidth = width;
				cubemapHeight = height;
				stbi_image_free(data);
			}
			else
			{
				std::cout << "Cubemap texture failed to load at path: " << faces[i] << std::endl;
				stbi_image_free(data);
			}
		}
		glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
		glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
		glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
		glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
		glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_WRAP_R, GL_CLAMP_TO_EDGE);

		return textureID;
	}

	unsigned int initDynamicCubemap() {

		//glEnable(GL_DEPTH_TEST);
		glEnable(GL_CULL_FACE);

		unsigned int depthRenderbuffer, dynamicCubeTex;
		int width, height, nrComponents;

		// Create empty cubemap
		glGenTextures(1, &dynamicCubeTex);
		glBindTexture(GL_TEXTURE_CUBE_MAP, dynamicCubeTex);

		glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
		glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
		glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
		glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
		glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_WRAP_R, GL_CLAMP_TO_EDGE);

		// Allocate space for each side of the cube map
		for (GLuint i = 0; i < 6; ++i)
		{
			glTexImage2D(GL_TEXTURE_CUBE_MAP_POSITIVE_X + i, 0, GL_RGB, cubemapSize,
				cubemapSize, 0, GL_RGB, GL_UNSIGNED_BYTE, NULL);
		}

		// Create framebuffer
		glGenFramebuffers(1, &framebuffer);
		glBindFramebuffer(GL_FRAMEBUFFER, framebuffer);
		glGenRenderbuffers(1, &dcRenderBuffer);
		glBindRenderbuffer(GL_RENDERBUFFER, dcRenderBuffer);
		glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH_COMPONENT16, cubemapSize, cubemapSize);
		// Attach one of the faces of the cubemap texture to current framebuffer
		glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0,GL_TEXTURE_CUBE_MAP_POSITIVE_X, dynamicCubeTex, 0);
		// Attach depth buffer to framebuffer
		glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_RENDERBUFFER, dcRenderBuffer);
		// Attach only the +X cubemap texture (for completeness)
		glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_CUBE_MAP_POSITIVE_X, dynamicCubeTex, 0);
		
		// Check if current configuration of framebuffer is correct
		if (glCheckFramebufferStatus(GL_FRAMEBUFFER) != GL_FRAMEBUFFER_COMPLETE)
			std::cout << "Framebuffer not complete!" << std::endl;

		// Set default framebuffer
		glBindFramebuffer(GL_FRAMEBUFFER, 0);

		

		return dynamicCubeTex;
	}

	unsigned int initDynamicNoSphereMap() {

		//glEnable(GL_DEPTH_TEST);
		glEnable(GL_CULL_FACE);

		
		unsigned int depthRenderbuffer, dynamicCubeTex;
		int width, height, nrComponents;

		// Create empty cubemap
		glGenTextures(1, &dynamicCubeTex);
		glBindTexture(GL_TEXTURE_CUBE_MAP, dynamicCubeTex);
		glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
		glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
		glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
		glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
		glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_WRAP_R, GL_CLAMP_TO_EDGE);

		// Allocate space for each side of the cube map
		for (GLuint i = 0; i < 6; ++i)
		{
			
			glTexImage2D(GL_TEXTURE_CUBE_MAP_POSITIVE_X + i, 0, GL_RGB, cubemapSize,
				cubemapSize, 0, GL_RGB, GL_UNSIGNED_BYTE, NULL);
		}

		// Create framebuffer
		glGenFramebuffers(1, &noSphereFrameBuffer);
		glBindFramebuffer(GL_FRAMEBUFFER, noSphereFrameBuffer);
		glGenRenderbuffers(1, &noSphereRenderBuffer);
		glBindRenderbuffer(GL_RENDERBUFFER, noSphereRenderBuffer);
		glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH_COMPONENT16, cubemapSize, cubemapSize);
		// Attach one of the faces of the cubemap texture to current framebuffer
		glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_CUBE_MAP_POSITIVE_X, dynamicCubeTex, 0);
		// Attach depth buffer to framebuffer
		glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_RENDERBUFFER, noSphereRenderBuffer);
		// Attach only the +X cubemap texture (for completeness)
		glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_CUBE_MAP_POSITIVE_X, dynamicCubeTex, 0);

		// Check if current configuration of framebuffer is correct
		if (glCheckFramebufferStatus(GL_FRAMEBUFFER) != GL_FRAMEBUFFER_COMPLETE)
			std::cout << "Framebuffer not complete!" << std::endl;

		// Set default framebuffer
		glBindFramebuffer(GL_FRAMEBUFFER, 0);

		return dynamicCubeTex;
	}
};