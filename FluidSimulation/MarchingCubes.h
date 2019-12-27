#pragma once

#include <glad\glad.h>
#include <GLFW\glfw3.h>

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

#include "Triangle.h"
#include "Models.h"
#include "Particle.h"
#include "Triangulator.h"
#include "Noise.h"
#include <vector>

#include <omp.h>
#include <memory>
#include <array>
#include <algorithm>
#include <functional>
#include <map>

#include <CL\cl.hpp>

class MarchingCubes {

public:


	SphGrid *grid;

	//float isolevel;
	float gridSize;

	//OpenCL
	//Device Info
	cl::Device device;
	//Program Info
	cl::Program program;
	//Context Info
	cl::Context* context;

	//Buffer Info
	cl::CommandQueue queue;
	cl::CommandQueue queueInitializeArrays;
	cl::Buffer isoBuffer_CL;
	cl::Buffer bufferPos_CL;
	cl::Buffer xBuff_CL;
	cl::Buffer yBuff_CL;
	cl::Buffer zBuff_CL;
	cl::Buffer xyzBuff_CL;

	cl::Buffer parBufferPos_CL;
	cl::Buffer particles_CL;
	cl::Kernel kernel;
	cl::Kernel zeroKernel;

	typedef struct Vec3f {
		float x;
		float y;
		float z;

	} Vec3f;

	Vec3f *xyzBuff;
	Vec3f *parPosBuff;
	Vec3f *particlePosBuff;


	float *zeroBuff;
	float *isoBuff;
	float *xBuff;
	float *yBuff;
	float *zBuff;

	unsigned int VAO, VBO, EBO;
	unsigned int texture, diffuseTexture;

	unsigned int LVAO, LVBO, LEBO;

	bool rendered = false;
	bool linerendered = false;
	struct Rect3
	{
		Vector3 min;
		Vector3 size;
	};

	Rect3 domain;

	Triangulator triangulator;

	vector<Vector3> voxelGrads;

	vector<float> voxelGradX;
	vector<float> voxelGradY;
	vector<float> voxelGradZ;
	vector<float> voxels;
	vector<float> prevoxels;

	vector<float> clVertX;
	vector<float> clVertY;
	vector<float> clVertZ;

	vector<float> clGradX;
	vector<float> clGradY;
	vector<float> clGradZ;

	cl::Buffer posXBuf;
	cl::Buffer posYBuf;
	cl::Buffer posZBuf;

	cl::Buffer startsBuf;
	cl::Buffer endsBuf;

	vector<int> gridIndices;
	vector<Vector3> gridVertices;
	vector<Models::Vertex> gridVertexes;

	vector<int> normalPointIndices;
	vector<Vector3> normalPointVertices;
	vector<Models::Vertex> normalPointVertexes;

	int edge_indices[12];
	float isoLevel;

	int counter = 0;

	vector<float> A;
	vector<float> B;
	vector<float> C;
	vector<float> D;
	vector<float> E;
	vector<float> F;
	vector<float> G;
	vector<float> H;

	cl::Buffer voxelBuf;
	cl::Buffer voxelGradXBuf;
	cl::Buffer voxelGradYBuf;
	cl::Buffer voxelGradZBuf;


	cl::Buffer clVertXBuf;
	cl::Buffer clVertYBuf;
	cl::Buffer clVertZBuf;

	cl::Buffer clGradXBuf;
	cl::Buffer clGradYBuf;
	cl::Buffer clGradZBuf;

	int gridNum;

	cl::Kernel metaballKernel;
	cl::Kernel mcKernel;

	Vector3 camPos;

	int runOpenCL = 0;

public:

	MarchingCubes() {

	}

	MarchingCubes(SphGrid *grid) {
		this->grid = grid;

		domain.min = { -10.f, -10.f, -10.f };
		domain.size = { 20.f, 20.f, 20.f };

		triangulator = Triangulator();
		triangulator.radius = 1;

		gridNum = ((grid->vars.kMcNx + 1)*(grid->vars.kMcNy + 1)*(grid->vars.kMcNz + 1));
	}

	void setupCL() {

		string rOpenCL;
		cout << "Do you want to run OpenCL (Y/N)?" ;
		cin >> rOpenCL;
		cout << endl;

		if (rOpenCL == "y" || rOpenCL == "Y") {
			runOpenCL = 1;
		}
		else {
			runOpenCL = 0;
			return;
		}

		//------------platforms---------------
		vector< cl::Platform > platforms;
		cl::Platform::get(&platforms);
		cout << "Platform number is: " << platforms.size() << std::endl;
		string s;
		for (int i = 0; i < platforms.size(); i++) {
			platforms[i].getInfo(CL_PLATFORM_NAME, &s);
			cout << i << " : " << s << endl;
		}

		cout << "Select a platform" << endl;
		int platform_selection = -1;
		//cin >> platform_selection;
		platform_selection = 0;

		cl::Platform selected_platform = platforms[platform_selection];

		//----------contexts-----------------
		cl_int err;
		cl_context_properties temp_cp[3] = { CL_CONTEXT_PLATFORM, (cl_context_properties)(selected_platform)(), 0 };
		context = new cl::Context(
			CL_DEVICE_TYPE_ALL,
			temp_cp,
			NULL,
			NULL,
			&err);


		//-----------devices----------------


		string type_s;

		string platformVendor;
		cl_uint cores;
		cl_uint clockf;
		cl_device_type  type;
		cl_ulong memsize;
		size_t workersize[3];

		type_s = "";

		vector<cl::Device> devices;
		devices = context->getInfo<CL_CONTEXT_DEVICES>();

		//-------device info-----------
		for (int k = 0; k < devices.size(); k++) {
			cl::Device temp_device = devices[k];
			temp_device.getInfo(CL_DEVICE_VENDOR, &platformVendor);
			temp_device.getInfo(CL_DEVICE_NAME, &s);
			temp_device.getInfo(CL_DEVICE_TYPE, &type);
			temp_device.getInfo(CL_DEVICE_MAX_COMPUTE_UNITS, &cores);
			temp_device.getInfo(CL_DEVICE_GLOBAL_MEM_SIZE, &memsize);
			temp_device.getInfo(CL_DEVICE_MAX_CLOCK_FREQUENCY, &clockf);
			temp_device.getInfo(CL_DEVICE_MAX_WORK_ITEM_SIZES, &workersize);
			type_s = "CPU ";
			if (type == CL_DEVICE_TYPE_GPU)
				type_s = "GPU ";
			if (type == CL_DEVICE_TYPE_ACCELERATOR)
				type_s = "accelerator ";
			if (type == CL_DEVICE_TYPE_DEFAULT)
				type_s = "default ";
			cout << k << " : " << "Device is by: " << platformVendor << "\n";

			cout << "device " << k << " is " << s << endl;
			cout << "type " << type_s << endl;
			cout << "cores " << cores << endl;
			cout << "clock frequency " << clockf / 1000.0 << " GHZ" << endl;
			cout << "memory size " << memsize / 1000000.0 << " MB" << endl;
			cout << "worker size x:" << workersize[0] << " y:" << workersize[1] << " z:" << workersize[2] << endl;
			cout << endl;

		}

		cout << "Select a device" << endl;
		int device_selection = -1;
		//cin >> device_selection;

		device_selection = 0;

		device = devices[device_selection];

		//----------programm----------------

		std::ifstream file("MetaballField.cl");
		std::string prog(
			std::istreambuf_iterator<char>(file), (std::istreambuf_iterator<char>())
		);
		cl::Program::Sources source(1, std::make_pair(prog.c_str(), prog.length() + 1));
		program = cl::Program(*context, source);


		string buildlog;
		program.build(devices, "");
		program.getBuildInfo(device, (cl_program_build_info)CL_PROGRAM_BUILD_LOG, &buildlog);
		cout << buildlog << endl;

		
		voxels.clear();
		voxels.resize(gridNum, -1);

		voxelGradX.clear();
		voxelGradX.resize(gridNum, 0);

		voxelGradY.clear();
		voxelGradY.resize(gridNum, 0);

		voxelGradZ.clear();
		voxelGradZ.resize(gridNum, 0);

		voxelGrads.clear();
		voxelGrads.resize(gridNum, Vector3());

		clVertX.clear();
		clVertX.resize(gridNum * 15+16, -10000);

		clVertY.clear();
		clVertY.resize(gridNum * 15+16, -10000);

		clVertZ.clear();
		clVertZ.resize(gridNum * 15+16, -10000);

		clGradX.clear();
		clGradX.resize(gridNum * 15+16, -10000);

		clGradY.clear();
		clGradY.resize(gridNum * 15+16, -10000);

		clGradZ.clear();
		clGradZ.resize(gridNum * 15+16, -10000);


		allocArraysCL();

		
		//initArraysCL();

		//mapArraysCL();

		prepareMetaballFieldKernelCL();
		prepareMCKernelCL();
	}

	glm::vec3 getClosestPointToCamera(Vector3 camera) {

		Vector3 camerPosition = Vector3(camera.x, camera.y, camera.z);

		float minDist = FLT_MAX;

		glm::vec3 minVec;
		int N = triangulator.mesh.m_vertices.size();
		for (int i = 0; i < N;i++) {

			Vector3 v = triangulator.mesh.m_vertices[i];

			float dist = (v - camerPosition).lengthSq();

			if (dist < minDist) {
				minDist = dist;

				minVec = glm::vec3(v.x,v.y,v.z);
			}

		}

		return minVec;

	}

	glm::vec3 getCenterMass() {

		Vector3 cm_ = triangulator.mesh.cm();
		return glm::vec3(cm_.x, cm_.y, cm_.z);
	}

	Vector3 getCM() {

		return triangulator.mesh.cm();


		Vector3 cm = Vector3(0, 0, 0);
		/*for (int i = 0; i < grid->nParticles;i++) {

			cm += grid->particles[i].pos;

		}

		cm /= (float)grid->nParticles;*/

		/*for (int i = 0; i < triangulator.mesh.m_vertices.size();i++) {

			cm += triangulator.mesh.m_vertices[i];

		}

		cm /= (float)triangulator.mesh.m_vertices.size();*/

		for (int z = 0; z < (grid->vars.kMcNz + 1); z++) {
			for (int y = 0; y < (grid->vars.kMcNy + 1); y++) {
				for (int x = 0; x < (grid->vars.kMcNx + 1); x++) {

					Vector3 cell = Vector3(x,y,z)*grid->vars.kMcStep;

					cm += cell;

				}
			}
		}

		cm /= (float)voxels.size();

		return cm;

		
	}

	

	void setBuffersCL() {

		//createGrid();
		/*
		//---------buffers-------------------

		int isoLength = 8 * grids.size();
		isoBuff = new float[isoLength];
		xyzBuff = new Vec3f[isoLength];
		zeroBuff = new float[isoLength];
		particlePosBuff = new Vec3f[grid->nParticles];
		int pN = grid->nParticles;
		for (int i = 0; i < pN; i++) {
			particlePosBuff[i] = getVec3f(grid->particles[i].pos);
		}

		isoBuffer_CL = cl::Buffer(*context, CL_MEM_READ_WRITE, sizeof(float)*isoLength);
		xyzBuff_CL = cl::Buffer(*context, CL_MEM_READ_WRITE, sizeof(Vec3f)*isoLength);
		//bufferPos_CL = cl::Buffer(*context, CL_MEM_READ_WRITE, sizeof(Vec3f) * particles.size());
		parBufferPos_CL = cl::Buffer(*context, CL_MEM_READ_WRITE, sizeof(Vec3f)*grid->nParticles);
		particles_CL = cl::Buffer(*context, CL_MEM_READ_WRITE, sizeof(Vec3f)*grid->nParticles);

		//reset isovalues
		for (int i = 0; i < grids.size(); i++) {
			for (int j = 0; j < 8; j++) {
				isoBuff[i * 8 + j] = 0;
				zeroBuff[i * 8 + j] = 0;
				xyzBuff[i * 8 + j].x = grids[i].p[j].x;
				xyzBuff[i * 8 + j].y = grids[i].p[j].y;
				xyzBuff[i * 8 + j].z = grids[i].p[j].z;
			}
		}

		int N = grids.size();


		float radius = grid->h;// grid->particleSize / 2;


		int maxSize = (isoLength);




		//---------CommandQueue--------------


		cl::Event* event = new cl::Event();

		queue = cl::CommandQueue(*context, device);

		queue.enqueueWriteBuffer(xyzBuff_CL, CL_TRUE, 0, sizeof(Vec3f)*isoLength, xyzBuff);
		queue.enqueueWriteBuffer(isoBuffer_CL, CL_TRUE, 0, sizeof(float)*isoLength, isoBuff);
		queue.enqueueWriteBuffer(particles_CL, CL_TRUE, 0, sizeof(Vec3f)*pN, particlePosBuff);

		//---------Kernel--------------

		kernel = cl::Kernel(program, "MarchingCubeCL");
		kernel.setArg(0, isoBuffer_CL);
		kernel.setArg(1, xyzBuff_CL);
		kernel.setArg(2, particles_CL);
		kernel.setArg(3, radius);
		kernel.setArg(4, maxSize);
		kernel.setArg(5, pN);

		zeroKernel = cl::Kernel(program, "zeroBuffer");
		zeroKernel.setArg(0, isoBuffer_CL);
		zeroKernel.setArg(1, maxSize);

		setupRenderBuffers();

		*/
	}

	

	Vec3f getVec3f(Vector3 p) {
		Vec3f v = Vec3f();
		v.x = p.x;
		v.y = p.y;
		v.z = p.z;

		return v;
	}

	void calcTangentSpace() {

		int N = triangulator.mesh.m_triangles.size();
		for (int i = 0; i < N;i++) {

			Triangle t = triangulator.mesh.m_triangles[i];
			int i1 = t.i1;
			int i2 = t.i2;
			int i3 = t.i3;

			Vector3 v1 = t.v1;
			Vector3 v2 = t.v2;
			Vector3 v3 = t.v3;

			Vector2 w1 = triangulator.mesh.m_vertexes[i1].TexCoords;
			Vector2 w2 = triangulator.mesh.m_vertexes[i2].TexCoords;
			Vector2 w3 = triangulator.mesh.m_vertexes[i3].TexCoords;

			float x1 = v2.x - v1.x;
			float x2 = v3.x - v1.x;
			float y1 = v2.y - v1.y;
			float y2 = v3.y - v1.y;
			float z1 = v2.z - v1.z;
			float z2 = v3.z - v1.z;

			float s1 = w2.x - w1.x;
			float s2 = w3.x - w1.x;
			float t1 = w2.y - w1.y;
			float t2 = w3.y - w1.y;

			float r = 1.0f / (s1 * t2 - s2 * t1);
			Vector3 sdir((t2 * x1 - t1 * x2) * r, (t2 * y1 - t1 * y2) * r,
				(t2 * z1 - t1 * z2) * r);
			Vector3 tdir((s1 * x2 - s2 * x1) * r, (s1 * y2 - s2 * y1) * r,
				(s1 * z2 - s2 * z1) * r);

			triangulator.mesh.m_vertexes[i1].Tangent += sdir;
			triangulator.mesh.m_vertexes[i2].Tangent += sdir;
			triangulator.mesh.m_vertexes[i3].Tangent += sdir;

			triangulator.mesh.m_vertexes[i1].BiTangent += tdir;
			triangulator.mesh.m_vertexes[i2].BiTangent += tdir;
			triangulator.mesh.m_vertexes[i3].BiTangent += tdir;
		}

		N = triangulator.mesh.m_vertices.size();
		for (int i = 0; i < N;i++) {
			Vector3 n = triangulator.mesh.m_vertexes[i].Normal;
			Vector3 t = triangulator.mesh.m_vertexes[i].Tangent;

			triangulator.mesh.m_vertexes[i].Tangent = (t-n*n.dot(t)).normalize();
			bool cond = n.cross(t).dot(triangulator.mesh.m_vertexes[i].BiTangent);
			triangulator.mesh.m_vertexes[i].wTangent = (cond < 0.0F) ? -1.0F : 1.0F;

			triangulator.mesh.m_vertexes[i].BiTangent = (n.cross(triangulator.mesh.m_vertexes[i].Tangent))*triangulator.mesh.m_vertexes[i].wTangent;
		}
	}

	void setupRenderBuffers() {
		glGenVertexArrays(1, &VAO);
		glGenBuffers(1, &VBO);
		glGenBuffers(1, &EBO);

		glBindVertexArray(VAO);
		// load data into vertex buffers
		glBindBuffer(GL_ARRAY_BUFFER, VBO);

		//vector<float> aVertices = getVertices();
		//cout << vertexes.size() << endl;

		glBufferData(GL_ARRAY_BUFFER, triangulator.mesh.m_vertexes.size() * sizeof(Models::Vertex), triangulator.mesh.m_vertexes.data(), GL_STREAM_DRAW);
		//glBufferData(GL_ARRAY_BUFFER, sizeof(aVertices.data()), aVertices.data(), GL_STATIC_DRAW);

		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
		glBufferData(GL_ELEMENT_ARRAY_BUFFER, triangulator.mesh.m_indices.size() * sizeof(unsigned int), triangulator.mesh.m_indices.data(), GL_STREAM_DRAW);
		//glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(indices.data())*sizeof(unsigned int), &indices.data()[0], GL_STATIC_DRAW);

		// set the vertex attribute pointers
		// Vector3 Positions
		glEnableVertexAttribArray(0);
		glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(Models::Vertex), (void*)0);
		//glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 9*sizeof(float), (void*)0);

		// Vector3 Normals
		glEnableVertexAttribArray(1);
		glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, sizeof(Models::Vertex), (void*)offsetof(Models::Vertex, Models::Vertex::Normal));
		//glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 9*sizeof(float), (void*)(3 * sizeof(float)));

		// Vector3 (x,y,0) Texture Coords
		glEnableVertexAttribArray(2);
		glVertexAttribPointer(2, 2, GL_FLOAT, GL_FALSE, sizeof(Models::Vertex), (void*)offsetof(Models::Vertex, Models::Vertex::TexCoords));
		//glVertexAttribPointer(2, 3, GL_FLOAT, GL_FALSE, 9*sizeof(float), (void*)(6 * sizeof(float)));

		glEnableVertexAttribArray(4);
		glVertexAttribPointer(4, 3, GL_FLOAT, GL_FALSE, sizeof(Models::Vertex), (void*)offsetof(Models::Vertex, Models::Vertex::Tangent));

		glEnableVertexAttribArray(5);
		glVertexAttribPointer(5, 3, GL_FLOAT, GL_FALSE, sizeof(Models::Vertex), (void*)offsetof(Models::Vertex, Models::Vertex::BiTangent));

		//glBindBuffer(GL_ARRAY_BUFFER, 0);
		//glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
		//(void*)(3 * sizeof(float))

		glBindVertexArray(0);
	}

	void showme(Vector3 v, string s = "") {

		cout << s << "|"<< v.x << "|" << v.y << "|" << v.z << endl;
	}

	inline int offset_3d(Vector3 p)
	{
		//p += Vector3(-kMcMargin, -kMcMargin, -kMcMargin);

		return (p.z * (grid->vars.kMcNy + 1) + p.y) * (grid->vars.kMcNx + 1) + p.x;
	}

	int offset_3d(float i, float j, float k){

		int y = grid->vars.kMcNy;
		int x = grid->vars.kMcNx;
		int z = grid->vars.kMcNz;

		//int y = kMcNx;
		//int x = kMcNy;

		if (i < 0)i = 0;
		if (j < 0)j = 0;
		if (k < 0)k = 0;

		if (i > grid->vars.kMcNx)i = grid->vars.kMcNx;
		if (j > grid->vars.kMcNy)j = grid->vars.kMcNy;
		if (k > grid->vars.kMcNz)k = grid->vars.kMcNz;

		return (int)((k*(y+1) + j)*(x+1) + i);
		//return i*y*z + j*z + k;

		//return (j*x+k)*z + i;
	}

	

	void generateVoxels() {
		Noise n2d(0);
		for (int z = 0; z < (grid->vars.kMcNz+1); z++) {
			for (int y = 0; y < (grid->vars.kMcNy + 1); y++) {
				for (int x = 0; x < (grid->vars.kMcNx + 1); x++) {
					float fy = (float)y / (grid->vars.kMcNy + 1);
					int offset = offset_3d(Vector3( x, y, z ));
					float v = n2d.get(x / 16.0f, z / 16.0f) * 0.25f;
					//cout << v << endl;
					//float v_ = calcDensityAt(x, y, z);
					voxels[offset] = fy - 0.25f - v;
					//voxels[offset] = v_;// fy - 0.25f - v_;// -1;
				}
			}
		}
	}

	void normalizeVoxels(float newmin = -1, float newmax = 1) {

		float minV = FLT_MAX;
		float maxV = -FLT_MAX;

		for (int i = 0; i < voxels.size();i++) {

			float v = voxels[i];

			if (v < minV)minV = v;
			if (v > maxV)maxV = v;
		}

		for (int i = 0; i < voxels.size(); i++) {

			voxels[i] = (voxels[i] - minV)*(newmax - newmin) / (maxV - minV) + newmin;
		}

	}

	void setupGridCells() {
		int count = 0;
		gridVertexes.clear();
		gridVertices.clear();
		gridIndices.clear();
		for (int z = 0; z < (grid->vars.kMcNz); z++) {
			for (int y = 0; y < (grid->vars.kMcNy); y++) {
				for (int x = 0; x < (grid->vars.kMcNx); x++) {

					Vector3 pos = Vector3(x, y, z)*grid->vars.kMcStep;
					
					Vector3 posX = Vector3(x + 1, y, z)*grid->vars.kMcStep;
					Vector3 posY = Vector3(x, y+1, z)*grid->vars.kMcStep;
					Vector3 posZ = Vector3(x, y, z+1)*grid->vars.kMcStep;

					addGridCellPointLine(pos, count);
					addGridCellPointLine(posX, count);

					addGridCellPointLine(pos, count);
					addGridCellPointLine(posY, count);

					addGridCellPointLine(pos, count);
					addGridCellPointLine(posZ, count);
				}
			}
		}


		glGenVertexArrays(1, &LVAO);
		glGenBuffers(1, &LVBO);
		glGenBuffers(1, &LEBO);

		glBindVertexArray(LVAO);
		// load data into vertex buffers
		glBindBuffer(GL_ARRAY_BUFFER, LVBO);

		//vector<float> aVertices = getVertices();
		//cout << vertexes.size() << endl;

		glBufferData(GL_ARRAY_BUFFER, gridVertexes.size() * sizeof(Models::Vertex), &gridVertexes[0], GL_STATIC_DRAW);
		//glBufferData(GL_ARRAY_BUFFER, sizeof(aVertices.data()), aVertices.data(), GL_STATIC_DRAW);

		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, LEBO);
		glBufferData(GL_ELEMENT_ARRAY_BUFFER, gridIndices.size() * sizeof(unsigned int), &gridIndices[0], GL_STATIC_DRAW);
		

		// set the vertex attribute pointers
		// Vector3 Positions
		glEnableVertexAttribArray(0);
		glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(Models::Vertex), (void*)0);
		//glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 9*sizeof(float), (void*)0);

		// Vector3 Normals
		glEnableVertexAttribArray(1);
		glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, sizeof(Models::Vertex), (void*)offsetof(Models::Vertex, Models::Vertex::Normal));
		//glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 9*sizeof(float), (void*)(3 * sizeof(float)));

		// Vector3 (x,y,0) Texture Coords
		glEnableVertexAttribArray(2);
		glVertexAttribPointer(2, 2, GL_FLOAT, GL_FALSE, sizeof(Models::Vertex), (void*)offsetof(Models::Vertex, Models::Vertex::TexCoords));
		//glVertexAttribPointer(2, 3, GL_FLOAT, GL_FALSE, 9*sizeof(float), (void*)(6 * sizeof(float)));
		

		//(void*)(3 * sizeof(float))

		glBindVertexArray(0);

	}
	
	void addGridCellPointLine(Vector3 p, int &count) {

		gridVertices.push_back(p);

		Models::Vertex v = Models::Vertex();
		v.Position = p;
		v.Normal = -p;
		gridVertexes.push_back(v);

		gridIndices.push_back(count);

		count++;
	}

	void setupNormalPointLines() {
		int count = 0;
	    normalPointVertexes.clear();
		normalPointVertices.clear();
		normalPointIndices.clear();

		int N = triangulator.mesh.m_gradients.size();
		
		const float lineLength = 5;
		for (int i = 0; i < N;i++) {

			Vector3 v = triangulator.mesh.m_vertices[i];
			Vector3 n = triangulator.mesh.m_gradients[i];

			addNormalPointLine(v, count);
			addNormalPointLine(v+lineLength*n.normalize(),count);
		}

		glGenVertexArrays(1, &LVAO);
		glGenBuffers(1, &LVBO);
		glGenBuffers(1, &LEBO);

		glBindVertexArray(LVAO);
		// load data into vertex buffers
		glBindBuffer(GL_ARRAY_BUFFER, LVBO);

		//vector<float> aVertices = getVertices();
		//cout << vertexes.size() << endl;

		glBufferData(GL_ARRAY_BUFFER, normalPointVertexes.size() * sizeof(Models::Vertex), &normalPointVertexes[0], GL_STATIC_DRAW);
		//glBufferData(GL_ARRAY_BUFFER, sizeof(aVertices.data()), aVertices.data(), GL_STATIC_DRAW);

		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, LEBO);
		glBufferData(GL_ELEMENT_ARRAY_BUFFER, normalPointIndices.size() * sizeof(unsigned int), &normalPointIndices[0], GL_STATIC_DRAW);


		// set the vertex attribute pointers
		// Vector3 Positions
		glEnableVertexAttribArray(0);
		glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(Models::Vertex), (void*)0);
		//glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 9*sizeof(float), (void*)0);

		// Vector3 Normals
		glEnableVertexAttribArray(1);
		glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, sizeof(Models::Vertex), (void*)offsetof(Models::Vertex, Models::Vertex::Normal));
		//glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 9*sizeof(float), (void*)(3 * sizeof(float)));

		// Vector3 (x,y,0) Texture Coords
		glEnableVertexAttribArray(2);
		glVertexAttribPointer(2, 2, GL_FLOAT, GL_FALSE, sizeof(Models::Vertex), (void*)offsetof(Models::Vertex, Models::Vertex::TexCoords));
		//glVertexAttribPointer(2, 3, GL_FLOAT, GL_FALSE, 9*sizeof(float), (void*)(6 * sizeof(float)));


		//(void*)(3 * sizeof(float))

		glBindVertexArray(0);
	}

	void addNormalPointLine(Vector3 p, int &count) {

		normalPointVertices.push_back(p);

		Models::Vertex v = Models::Vertex();
		v.Position = p;
		v.Normal = -p;
		normalPointVertexes.push_back(v);

		normalPointIndices.push_back(count);

		count++;

	}

	void renderNormalPointLines() {


		setupNormalPointLines();

		glBindVertexArray(LVAO);


		glDrawArrays(GL_LINES, 0, normalPointIndices.size());
		glBindVertexArray(0);

	}

	void renderGridCells(){

		if (linerendered == false) {
			//setupGridCells();
			linerendered = true;
		}
		// draw mesh
		glBindVertexArray(LVAO);

		//glActiveTexture(GL_TEXTURE0);
		//glBindTexture(GL_TEXTURE_2D, texture);
		//glDrawElements(GL_TRIANGLES, triangulator.mesh.m_indices.size(), GL_UNSIGNED_INT, 0);
		//glDrawElements(GL_TRIANGLES, triangulator.mesh.m_vertices.size()/9, GL_UNSIGNED_INT, 0);

		//works
		glDrawArrays(GL_LINES, 0, gridIndices.size());
		glBindVertexArray(0);


	}

	void initCornerArrays(float value=-1) {
		A.clear();
		A.resize((grid->vars.kMcNx + 1)*(grid->vars.kMcNy + 1)*(grid->vars.kMcNz + 1), value);

		B.clear();
		B.resize((grid->vars.kMcNx + 1)*(grid->vars.kMcNy + 1)*(grid->vars.kMcNz + 1), value);

		C.clear();
		C.resize((grid->vars.kMcNx + 1)*(grid->vars.kMcNy + 1)*(grid->vars.kMcNz + 1), value);

		D.clear();
		D.resize((grid->vars.kMcNx + 1)*(grid->vars.kMcNy + 1)*(grid->vars.kMcNz + 1), value);

		E.clear();
		E.resize((grid->vars.kMcNx + 1)*(grid->vars.kMcNy + 1)*(grid->vars.kMcNz + 1), value);

		F.clear();
		F.resize((grid->vars.kMcNx + 1)*(grid->vars.kMcNy + 1)*(grid->vars.kMcNz + 1), value);

		G.clear();
		G.resize((grid->vars.kMcNx + 1)*(grid->vars.kMcNy + 1)*(grid->vars.kMcNz + 1), value);

		H.clear();
		H.resize((grid->vars.kMcNx + 1)*(grid->vars.kMcNy + 1)*(grid->vars.kMcNz + 1), value);
	}

	void prepareMCKernelCL() {
		mcKernel = cl::Kernel(program, "marchingCubes");

		mcKernel.setArg(0, grid->vars.kMcNx);
		mcKernel.setArg(1, grid->vars.kMcNy);
		mcKernel.setArg(2, grid->vars.kMcNz);
		mcKernel.setArg(3, grid->vars.kMcStep);
		mcKernel.setArg(4, grid->vars.kMcMargin);

		mcKernel.setArg(5, voxelGradXBuf);
		mcKernel.setArg(6, voxelGradYBuf);
		mcKernel.setArg(7, voxelGradZBuf);

		mcKernel.setArg(8, voxelBuf);

		clVertXBuf = cl::Buffer(*context, CL_MEM_READ_WRITE | CL_MEM_HOST_READ_ONLY | CL_MEM_COPY_HOST_PTR, sizeof(float)*((gridNum) * 15+16), clVertX.data());
		mcKernel.setArg(9, clVertXBuf);

		clVertYBuf = cl::Buffer(*context, CL_MEM_READ_WRITE | CL_MEM_HOST_READ_ONLY | CL_MEM_COPY_HOST_PTR, sizeof(float)*((gridNum) * 15+16), clVertY.data());
		mcKernel.setArg(10, clVertYBuf);

		clVertZBuf = cl::Buffer(*context, CL_MEM_READ_WRITE | CL_MEM_HOST_READ_ONLY | CL_MEM_COPY_HOST_PTR, sizeof(float)*((gridNum) * 15+16), clVertZ.data());
		mcKernel.setArg(11, clVertZBuf);

		clGradXBuf = cl::Buffer(*context, CL_MEM_READ_WRITE | CL_MEM_HOST_READ_ONLY | CL_MEM_COPY_HOST_PTR, sizeof(float)*((gridNum) * 15+16), clGradX.data());
		mcKernel.setArg(12, clGradXBuf);

		clGradYBuf = cl::Buffer(*context, CL_MEM_READ_WRITE | CL_MEM_HOST_READ_ONLY | CL_MEM_COPY_HOST_PTR, sizeof(float)*((gridNum) * 15+16), clGradY.data());
		mcKernel.setArg(13, clGradYBuf);

		clGradZBuf = cl::Buffer(*context, CL_MEM_READ_WRITE | CL_MEM_HOST_READ_ONLY | CL_MEM_COPY_HOST_PTR, sizeof(float)*((gridNum) * 15+16), clGradZ.data());
		mcKernel.setArg(14, clGradZBuf);

	}

	void executeMCKernelCL() {

		mcKernel.setArg(0, grid->vars.kMcNx);
		mcKernel.setArg(1, grid->vars.kMcNy);
		mcKernel.setArg(2, grid->vars.kMcNz);
		mcKernel.setArg(3, grid->vars.kMcStep);
		mcKernel.setArg(4, grid->vars.kMcMargin);

		mcKernel.setArg(5, voxelGradXBuf);
		mcKernel.setArg(6, voxelGradYBuf);
		mcKernel.setArg(7, voxelGradZBuf);

		mcKernel.setArg(8, voxelBuf);

		cl::Buffer clVertXBuf_ = cl::Buffer(*context, CL_MEM_READ_WRITE | CL_MEM_HOST_READ_ONLY | CL_MEM_USE_HOST_PTR, sizeof(float)*((gridNum) * 15+16), clVertX.data());
		cl::Buffer clVertYBuf_ = cl::Buffer(*context, CL_MEM_READ_WRITE | CL_MEM_HOST_READ_ONLY | CL_MEM_USE_HOST_PTR, sizeof(float)*((gridNum) * 15+16), clVertY.data());
		cl::Buffer clVertZBuf_ = cl::Buffer(*context, CL_MEM_READ_WRITE | CL_MEM_HOST_READ_ONLY | CL_MEM_USE_HOST_PTR, sizeof(float)*((gridNum) * 15+16), clVertZ.data());
		cl::Buffer clGradXBuf_ = cl::Buffer(*context, CL_MEM_READ_WRITE | CL_MEM_HOST_READ_ONLY | CL_MEM_USE_HOST_PTR, sizeof(float)*((gridNum) * 15+16), clGradX.data());
		cl::Buffer clGradYBuf_ = cl::Buffer(*context, CL_MEM_READ_WRITE | CL_MEM_HOST_READ_ONLY | CL_MEM_USE_HOST_PTR, sizeof(float)*((gridNum) * 15+16), clGradY.data());
		cl::Buffer clGradZBuf_ = cl::Buffer(*context, CL_MEM_READ_WRITE | CL_MEM_HOST_READ_ONLY | CL_MEM_USE_HOST_PTR, sizeof(float)*((gridNum) * 15+16), clGradZ.data());

		mcKernel.setArg(9, clVertXBuf_);
		mcKernel.setArg(10, clVertYBuf_);
		mcKernel.setArg(11, clVertZBuf_);
		mcKernel.setArg(12, clGradXBuf_);
		mcKernel.setArg(13, clGradYBuf_);
		mcKernel.setArg(14, clGradZBuf_);

		cl::CommandQueue queueMetaball(*context, device);
		
		//queueMetaball.enqueueNDRangeKernel(mcKernel, cl::NullRange, cl::NDRange(grid->vars.kMcNx, grid->vars.kMcNy, grid->vars.kMcNz));
		queueMetaball.enqueueNDRangeKernel(mcKernel, cl::NDRange(1,1,1), cl::NDRange(grid->vars.kMcNx-1, grid->vars.kMcNy-1, grid->vars.kMcNz-1));

		queueMetaball.enqueueReadBuffer(clVertXBuf_, GL_TRUE, 0, sizeof(float)*((gridNum) * 15+16), clVertX.data());
		queueMetaball.enqueueReadBuffer(clVertYBuf_, GL_TRUE, 0, sizeof(float)*((gridNum) * 15+16), clVertY.data());
		queueMetaball.enqueueReadBuffer(clVertZBuf_, GL_TRUE, 0, sizeof(float)*((gridNum) * 15+16), clVertZ.data());

		queueMetaball.enqueueReadBuffer(clGradXBuf_, GL_TRUE, 0, sizeof(float)*((gridNum) * 15+16), clGradX.data());
		queueMetaball.enqueueReadBuffer(clGradYBuf_, GL_TRUE, 0, sizeof(float)*((gridNum) * 15+16), clGradY.data());
		queueMetaball.enqueueReadBuffer(clGradZBuf_, GL_TRUE, 0, sizeof(float)*((gridNum) * 15+16), clGradZ.data());

		//#pragma omp parallel for
		for (int i = 0; i < ((gridNum) * 15+16); i += 3) {

			if (clVertX[i] == NULL || clVertX[i + 1] == NULL || clVertX[i + 2] == NULL) continue;
			if(clVertY[i] == NULL || clVertY[i + 1] == NULL || clVertY[i + 2] == NULL) continue;
			if(clVertZ[i] == NULL || clVertZ[i + 1] == NULL || clVertZ[i + 2] == NULL) continue;

			if (clVertX[i] == -10000 || clVertX[i + 1] == -10000 || clVertX[i + 2] == -10000) continue;
			if (clVertY[i] == -10000 || clVertY[i + 1] == -10000 || clVertY[i + 2] == -10000) continue;
			if (clVertZ[i] == -10000 || clVertZ[i + 1] == -10000 || clVertZ[i + 2] == -10000) continue;

			if (clGradX[i] == -10000 || clGradX[i + 1] == -10000 || clGradX[i + 2] == -10000) continue;
			if (clGradY[i] == -10000 || clGradY[i + 1] == -10000 || clGradY[i + 2] == -10000) continue;
			if (clGradZ[i] == -10000 || clGradZ[i + 1] == -10000 || clGradZ[i + 2] == -10000) continue;

			
			int last = triangulator.mesh.m_indices.size();

			triangulator.mesh.m_indices.push_back(last);
			triangulator.mesh.m_indices.push_back(last + 1);
			triangulator.mesh.m_indices.push_back(last + 2);

			triangulator.mesh.m_vertices.push_back(Vector3(clVertX[i], clVertY[i], clVertZ[i]));
			triangulator.mesh.m_vertices.push_back(Vector3(clVertX[i + 1], clVertY[i + 1], clVertZ[i + 1]));
			triangulator.mesh.m_vertices.push_back(Vector3(clVertX[i + 2], clVertY[i + 2], clVertZ[i + 2]));
			triangulator.mesh.m_gradients.push_back(Vector3(clGradX[i], clGradY[i], clGradZ[i]));
			triangulator.mesh.m_gradients.push_back(Vector3(clGradX[i + 1], clGradY[i + 1], clGradZ[i + 1]));
			triangulator.mesh.m_gradients.push_back(Vector3(clGradX[i + 2], clGradY[i + 2], clGradZ[i + 2]));

		}

		triangulator.polishTriangulation();
	}

	void marchingCubesCL() {
		cl::Kernel kernelMetaBall(program, "marchingCubes");

		kernelMetaBall.setArg(0, grid->vars.kMcNx);
		kernelMetaBall.setArg(1, grid->vars.kMcNy);
		kernelMetaBall.setArg(2, grid->vars.kMcNz);
		kernelMetaBall.setArg(3, grid->vars.kMcStep);
		kernelMetaBall.setArg(4, grid->vars.kMcMargin);

		//voxelGradXBuf = cl::Buffer(*context, CL_MEM_READ_WRITE | CL_MEM_USE_HOST_PTR, sizeof(float)*((grid->vars.kMcNx + 1)*(grid->vars.kMcNy + 1)*(grid->vars.kMcNz + 1)), voxelGradX.data());
		//voxelGradYBuf = cl::Buffer(*context, CL_MEM_READ_WRITE | CL_MEM_USE_HOST_PTR, sizeof(float)*((grid->vars.kMcNx + 1)*(grid->vars.kMcNy + 1)*(grid->vars.kMcNz + 1)), voxelGradY.data());
		//voxelGradZBuf = cl::Buffer(*context, CL_MEM_READ_WRITE | CL_MEM_USE_HOST_PTR, sizeof(float)*((grid->vars.kMcNx + 1)*(grid->vars.kMcNy + 1)*(grid->vars.kMcNz + 1)), voxelGradZ.data());

		kernelMetaBall.setArg(5, voxelGradXBuf);
		kernelMetaBall.setArg(6, voxelGradYBuf);
		kernelMetaBall.setArg(7, voxelGradZBuf);

		//cl::Buffer voxelBuf(*context, CL_MEM_READ_WRITE | CL_MEM_HOST_READ_ONLY | CL_MEM_USE_HOST_PTR, sizeof(float)*((grid->vars.kMcNx + 1)*(grid->vars.kMcNy + 1)*(grid->vars.kMcNz + 1)), voxels.data());
		kernelMetaBall.setArg(8, voxelBuf);

		cl::Buffer clVertXBuf(*context, CL_MEM_READ_WRITE | CL_MEM_USE_HOST_PTR, sizeof(float)*(gridNum) * 16, clVertX.data());
		kernelMetaBall.setArg(9, clVertXBuf);

		cl::Buffer clVertYBuf(*context, CL_MEM_READ_WRITE | CL_MEM_USE_HOST_PTR, sizeof(float)*(gridNum) * 16, clVertY.data());
		kernelMetaBall.setArg(10, clVertYBuf);

		cl::Buffer clVertZBuf(*context, CL_MEM_READ_WRITE | CL_MEM_USE_HOST_PTR, sizeof(float)*(gridNum) * 16, clVertZ.data());
		kernelMetaBall.setArg(11, clVertZBuf);

		cl::Buffer clGradXBuf(*context, CL_MEM_READ_WRITE | CL_MEM_USE_HOST_PTR, sizeof(float)*(gridNum) * 16, clGradX.data());
		kernelMetaBall.setArg(12, clGradXBuf);

		cl::Buffer clGradYBuf(*context, CL_MEM_READ_WRITE | CL_MEM_USE_HOST_PTR, sizeof(float)*(gridNum) * 16, clGradY.data());
		kernelMetaBall.setArg(13, clGradYBuf);

		cl::Buffer clGradZBuf(*context, CL_MEM_READ_WRITE | CL_MEM_USE_HOST_PTR, sizeof(float)*(gridNum) * 16, clGradZ.data());
		kernelMetaBall.setArg(14, clGradZBuf);

		cl::CommandQueue queueMetaball(*context, device);
		queueMetaball.enqueueNDRangeKernel(kernelMetaBall, cl::NullRange, cl::NDRange(grid->vars.kMcNz, grid->vars.kMcNy, grid->vars.kMcNx));

		queueMetaball.enqueueReadBuffer(clVertXBuf, GL_TRUE, 0, sizeof(float)*(gridNum) * 16, clVertX.data());
		queueMetaball.enqueueReadBuffer(clVertYBuf, GL_TRUE, 0, sizeof(float)*(gridNum) * 16, clVertY.data());
		queueMetaball.enqueueReadBuffer(clVertZBuf, GL_TRUE, 0, sizeof(float)*(gridNum) * 16, clVertZ.data());

		queueMetaball.enqueueReadBuffer(clGradXBuf, GL_TRUE, 0, sizeof(float)*(gridNum) * 16, clGradX.data());
		queueMetaball.enqueueReadBuffer(clGradYBuf, GL_TRUE, 0, sizeof(float)*(gridNum) * 16, clGradY.data());
		queueMetaball.enqueueReadBuffer(clGradZBuf, GL_TRUE, 0, sizeof(float)*(gridNum) * 16, clGradZ.data());

	
		//#pragma omp parallel for
		for (int i = 0; i < ((gridNum) * 15);i+=3) {
			
			if (clVertX[i] != -10000 && clVertX[i+1]!= -10000 && clVertX[i+2]!=-10000) {
				//cout << clVertX[i] << endl;

				int last = triangulator.mesh.m_indices.size();

				triangulator.mesh.m_indices.push_back(last);
				triangulator.mesh.m_indices.push_back(last+1);
				triangulator.mesh.m_indices.push_back(last+2);
				
				triangulator.mesh.m_vertices.push_back(Vector3(clVertX[i], clVertY[i], clVertZ[i]));
				triangulator.mesh.m_vertices.push_back(Vector3(clVertX[i+1], clVertY[i+1], clVertZ[i+1]));
				triangulator.mesh.m_vertices.push_back(Vector3(clVertX[i+2], clVertY[i+2], clVertZ[i+2]));
				triangulator.mesh.m_gradients.push_back(Vector3(clGradX[i], clGradY[i], clGradZ[i]));
				triangulator.mesh.m_gradients.push_back(Vector3(clGradX[i+1], clGradY[i+1], clGradZ[i+1]));
				triangulator.mesh.m_gradients.push_back(Vector3(clGradX[i+2], clGradY[i+2], clGradZ[i+2]));

			}
			
		}
		

		triangulator.polishTriangulation();
	}

	void mapArraysCL() {
		queueInitializeArrays = cl::CommandQueue(*context, device);
		voxelBuf = cl::Buffer(*context, CL_MEM_READ_WRITE | CL_MEM_HOST_READ_ONLY | CL_MEM_ALLOC_HOST_PTR, sizeof(float)*((grid->vars.kMcNx + 1)*(grid->vars.kMcNy + 1)*(grid->vars.kMcNz + 1)), voxels.data());
		voxelGradXBuf = cl::Buffer(*context, CL_MEM_READ_WRITE | CL_MEM_HOST_READ_ONLY | CL_MEM_ALLOC_HOST_PTR, sizeof(float)*((grid->vars.kMcNx + 1)*(grid->vars.kMcNy + 1)*(grid->vars.kMcNz + 1)), voxelGradX.data());
		voxelGradYBuf = cl::Buffer(*context, CL_MEM_READ_WRITE | CL_MEM_HOST_READ_ONLY | CL_MEM_ALLOC_HOST_PTR, sizeof(float)*((grid->vars.kMcNx + 1)*(grid->vars.kMcNy + 1)*(grid->vars.kMcNz + 1)), voxelGradY.data());
		voxelGradZBuf = cl::Buffer(*context, CL_MEM_READ_WRITE | CL_MEM_HOST_READ_ONLY | CL_MEM_ALLOC_HOST_PTR, sizeof(float)*((grid->vars.kMcNx + 1)*(grid->vars.kMcNy + 1)*(grid->vars.kMcNz + 1)), voxelGradZ.data());

		/*clVertXBuf = cl::Buffer(*context, CL_MEM_READ_WRITE | CL_MEM_HOST_READ_ONLY | CL_MEM_ALLOC_HOST_PTR, sizeof(float)*((grid->vars.kMcNx + 1)*(grid->vars.kMcNy + 1)*(grid->vars.kMcNz + 1)) * 16, clVertX.data());
		clVertYBuf = cl::Buffer(*context, CL_MEM_READ_WRITE | CL_MEM_HOST_READ_ONLY | CL_MEM_ALLOC_HOST_PTR, sizeof(float)*((grid->vars.kMcNx + 1)*(grid->vars.kMcNy + 1)*(grid->vars.kMcNz + 1)) * 16, clVertY.data());
		clVertZBuf = cl::Buffer(*context, CL_MEM_READ_WRITE | CL_MEM_HOST_READ_ONLY | CL_MEM_ALLOC_HOST_PTR, sizeof(float)*((grid->vars.kMcNx + 1)*(grid->vars.kMcNy + 1)*(grid->vars.kMcNz + 1)) * 16, clVertZ.data());

		clGradXBuf = cl::Buffer(*context, CL_MEM_READ_WRITE | CL_MEM_HOST_READ_ONLY | CL_MEM_ALLOC_HOST_PTR, sizeof(float)*((grid->vars.kMcNx + 1)*(grid->vars.kMcNy + 1)*(grid->vars.kMcNz + 1)) * 16, clGradX.data());
		clGradYBuf = cl::Buffer(*context, CL_MEM_READ_WRITE | CL_MEM_HOST_READ_ONLY | CL_MEM_ALLOC_HOST_PTR, sizeof(float)*((grid->vars.kMcNx + 1)*(grid->vars.kMcNy + 1)*(grid->vars.kMcNz + 1)) * 16, clGradY.data());
		clGradZBuf = cl::Buffer(*context, CL_MEM_READ_WRITE | CL_MEM_HOST_READ_ONLY | CL_MEM_ALLOC_HOST_PTR, sizeof(float)*((grid->vars.kMcNx + 1)*(grid->vars.kMcNy + 1)*(grid->vars.kMcNz + 1)) * 16, clGradZ.data());
		*/
	}

	void initArraysCL() {

		//cl::Kernel kernelInitializeArrays(program, "initArrays");

		voxelBuf = cl::Buffer(*context, CL_MEM_READ_WRITE | CL_MEM_HOST_READ_ONLY | CL_MEM_USE_HOST_PTR, sizeof(float)*((grid->vars.kMcNx + 1)*(grid->vars.kMcNy + 1)*(grid->vars.kMcNz + 1)), voxels.data());
		voxelGradXBuf = cl::Buffer(*context, CL_MEM_READ_WRITE | CL_MEM_HOST_READ_ONLY | CL_MEM_USE_HOST_PTR, sizeof(float)*((grid->vars.kMcNx + 1)*(grid->vars.kMcNy + 1)*(grid->vars.kMcNz + 1)), voxelGradX.data());
		voxelGradYBuf = cl::Buffer(*context, CL_MEM_READ_WRITE | CL_MEM_HOST_READ_ONLY | CL_MEM_USE_HOST_PTR, sizeof(float)*((grid->vars.kMcNx + 1)*(grid->vars.kMcNy + 1)*(grid->vars.kMcNz + 1)), voxelGradY.data());
		voxelGradZBuf = cl::Buffer(*context, CL_MEM_READ_WRITE | CL_MEM_HOST_READ_ONLY | CL_MEM_USE_HOST_PTR, sizeof(float)*((grid->vars.kMcNx + 1)*(grid->vars.kMcNy + 1)*(grid->vars.kMcNz + 1)), voxelGradZ.data());
		
		/*clVertXBuf = cl::Buffer(*context, CL_MEM_READ_WRITE | CL_MEM_HOST_READ_ONLY | CL_MEM_USE_HOST_PTR, sizeof(float)*((grid->vars.kMcNx + 1)*(grid->vars.kMcNy + 1)*(grid->vars.kMcNz + 1))*16, clVertX.data());
		clVertYBuf = cl::Buffer(*context, CL_MEM_READ_WRITE | CL_MEM_HOST_READ_ONLY | CL_MEM_USE_HOST_PTR, sizeof(float)*((grid->vars.kMcNx + 1)*(grid->vars.kMcNy + 1)*(grid->vars.kMcNz + 1))*16, clVertY.data());
		clVertZBuf = cl::Buffer(*context, CL_MEM_READ_WRITE | CL_MEM_HOST_READ_ONLY | CL_MEM_USE_HOST_PTR, sizeof(float)*((grid->vars.kMcNx + 1)*(grid->vars.kMcNy + 1)*(grid->vars.kMcNz + 1))*16, clVertZ.data());

		clGradXBuf = cl::Buffer(*context, CL_MEM_READ_WRITE | CL_MEM_HOST_READ_ONLY | CL_MEM_USE_HOST_PTR, sizeof(float)*((grid->vars.kMcNx + 1)*(grid->vars.kMcNy + 1)*(grid->vars.kMcNz + 1)) * 16, clGradX.data());
		clGradYBuf = cl::Buffer(*context, CL_MEM_READ_WRITE | CL_MEM_HOST_READ_ONLY | CL_MEM_USE_HOST_PTR, sizeof(float)*((grid->vars.kMcNx + 1)*(grid->vars.kMcNy + 1)*(grid->vars.kMcNz + 1)) * 16, clGradY.data());
		clGradZBuf = cl::Buffer(*context, CL_MEM_READ_WRITE | CL_MEM_HOST_READ_ONLY | CL_MEM_USE_HOST_PTR, sizeof(float)*((grid->vars.kMcNx + 1)*(grid->vars.kMcNy + 1)*(grid->vars.kMcNz + 1)) * 16, clGradZ.data());


		
		//kernelInitializeArrays.setArg(0, voxelBuf);
		queueInitializeArrays = cl::CommandQueue(*context, device);
		//queueInitializeArrays.enqueueNDRangeKernel(kernelInitializeArrays, cl::NullRange, cl::NDRange(((grid->vars.kMcNx + 1)*(grid->vars.kMcNy + 1)*(grid->vars.kMcNz + 1))));
		queueInitializeArrays.enqueueFillBuffer(voxelBuf,-1,0, sizeof(float)*((grid->vars.kMcNx + 1)*(grid->vars.kMcNy + 1)*(grid->vars.kMcNz + 1)));	
		queueInitializeArrays.enqueueFillBuffer(voxelGradXBuf, 0, 0, sizeof(float)*((grid->vars.kMcNx + 1)*(grid->vars.kMcNy + 1)*(grid->vars.kMcNz + 1)));
		queueInitializeArrays.enqueueFillBuffer(voxelGradYBuf, 0, 0, sizeof(float)*((grid->vars.kMcNx + 1)*(grid->vars.kMcNy + 1)*(grid->vars.kMcNz + 1)));
		queueInitializeArrays.enqueueFillBuffer(voxelGradZBuf, 0, 0, sizeof(float)*((grid->vars.kMcNx + 1)*(grid->vars.kMcNy + 1)*(grid->vars.kMcNz + 1)));
		
		queueInitializeArrays.enqueueFillBuffer(clVertXBuf, -10000, 0, sizeof(float)*((grid->vars.kMcNx + 1)*(grid->vars.kMcNy + 1)*(grid->vars.kMcNz + 1))*16);
		queueInitializeArrays.enqueueFillBuffer(clVertYBuf, -10000, 0, sizeof(float)*((grid->vars.kMcNx + 1)*(grid->vars.kMcNy + 1)*(grid->vars.kMcNz + 1))*16);
		queueInitializeArrays.enqueueFillBuffer(clVertZBuf, -10000, 0, sizeof(float)*((grid->vars.kMcNx + 1)*(grid->vars.kMcNy + 1)*(grid->vars.kMcNz + 1))*16);

		queueInitializeArrays.enqueueFillBuffer(clGradXBuf, -10000, 0, sizeof(float)*((grid->vars.kMcNx + 1)*(grid->vars.kMcNy + 1)*(grid->vars.kMcNz + 1)) * 16);
		queueInitializeArrays.enqueueFillBuffer(clGradYBuf, -10000, 0, sizeof(float)*((grid->vars.kMcNx + 1)*(grid->vars.kMcNy + 1)*(grid->vars.kMcNz + 1)) * 16);
		queueInitializeArrays.enqueueFillBuffer(clGradZBuf, -10000, 0, sizeof(float)*((grid->vars.kMcNx + 1)*(grid->vars.kMcNy + 1)*(grid->vars.kMcNz + 1)) * 16);

		
		queueInitializeArrays.enqueueReadBuffer(voxelBuf, GL_TRUE, 0, sizeof(float)*((grid->vars.kMcNx + 1)*(grid->vars.kMcNy + 1)*(grid->vars.kMcNz + 1)), voxels.data());
		queueInitializeArrays.enqueueReadBuffer(voxelGradXBuf, GL_TRUE, 0, sizeof(float)*((grid->vars.kMcNx + 1)*(grid->vars.kMcNy + 1)*(grid->vars.kMcNz + 1)), voxelGradX.data());
		queueInitializeArrays.enqueueReadBuffer(voxelGradYBuf, GL_TRUE, 0, sizeof(float)*((grid->vars.kMcNx + 1)*(grid->vars.kMcNy + 1)*(grid->vars.kMcNz + 1)), voxelGradY.data());
		queueInitializeArrays.enqueueReadBuffer(voxelGradZBuf, GL_TRUE, 0, sizeof(float)*((grid->vars.kMcNx + 1)*(grid->vars.kMcNy + 1)*(grid->vars.kMcNz + 1)), voxelGradZ.data());
		
		queueInitializeArrays.enqueueReadBuffer(clVertXBuf, GL_TRUE, 0, sizeof(float)*((grid->vars.kMcNx + 1)*(grid->vars.kMcNy + 1)*(grid->vars.kMcNz + 1))*16,clVertX.data());
		queueInitializeArrays.enqueueReadBuffer(clVertYBuf, GL_TRUE, 0, sizeof(float)*((grid->vars.kMcNx + 1)*(grid->vars.kMcNy + 1)*(grid->vars.kMcNz + 1))*16, clVertY.data());
		queueInitializeArrays.enqueueReadBuffer(clVertZBuf, GL_TRUE, 0, sizeof(float)*((grid->vars.kMcNx + 1)*(grid->vars.kMcNy + 1)*(grid->vars.kMcNz + 1))*16, clVertZ.data());

		queueInitializeArrays.enqueueReadBuffer(clGradXBuf, GL_TRUE, 0, sizeof(float)*((grid->vars.kMcNx + 1)*(grid->vars.kMcNy + 1)*(grid->vars.kMcNz + 1)) * 16, clGradX.data());
		queueInitializeArrays.enqueueReadBuffer(clGradYBuf, GL_TRUE, 0, sizeof(float)*((grid->vars.kMcNx + 1)*(grid->vars.kMcNy + 1)*(grid->vars.kMcNz + 1)) * 16, clGradY.data());
		queueInitializeArrays.enqueueReadBuffer(clGradZBuf, GL_TRUE, 0, sizeof(float)*((grid->vars.kMcNx + 1)*(grid->vars.kMcNy + 1)*(grid->vars.kMcNz + 1)) * 16, clGradZ.data());
		*/
	
	}

	void fillArraysCL() {
		
		//queueInitializeArrays.enqueueNDRangeKernel(kernelInitializeArrays, cl::NullRange, cl::NDRange(((grid->vars.kMcNx + 1)*(grid->vars.kMcNy + 1)*(grid->vars.kMcNz + 1))));
		/*queueInitializeArrays.enqueueFillBuffer(voxelBuf, -1, 0, sizeof(float)*((grid->vars.kMcNx + 1)*(grid->vars.kMcNy + 1)*(grid->vars.kMcNz + 1)));
		queueInitializeArrays.enqueueFillBuffer(voxelGradXBuf, 0, 0, sizeof(float)*((grid->vars.kMcNx + 1)*(grid->vars.kMcNy + 1)*(grid->vars.kMcNz + 1)));
		queueInitializeArrays.enqueueFillBuffer(voxelGradYBuf, 0, 0, sizeof(float)*((grid->vars.kMcNx + 1)*(grid->vars.kMcNy + 1)*(grid->vars.kMcNz + 1)));
		queueInitializeArrays.enqueueFillBuffer(voxelGradZBuf, 0, 0, sizeof(float)*((grid->vars.kMcNx + 1)*(grid->vars.kMcNy + 1)*(grid->vars.kMcNz + 1)));

		queueInitializeArrays.enqueueFillBuffer(clVertXBuf, -10000, 0, sizeof(float)*((grid->vars.kMcNx + 1)*(grid->vars.kMcNy + 1)*(grid->vars.kMcNz + 1)) * 16);
		queueInitializeArrays.enqueueFillBuffer(clVertYBuf, -10000, 0, sizeof(float)*((grid->vars.kMcNx + 1)*(grid->vars.kMcNy + 1)*(grid->vars.kMcNz + 1)) * 16);
		queueInitializeArrays.enqueueFillBuffer(clVertZBuf, -10000, 0, sizeof(float)*((grid->vars.kMcNx + 1)*(grid->vars.kMcNy + 1)*(grid->vars.kMcNz + 1)) * 16);

		queueInitializeArrays.enqueueFillBuffer(clGradXBuf, -10000, 0, sizeof(float)*((grid->vars.kMcNx + 1)*(grid->vars.kMcNy + 1)*(grid->vars.kMcNz + 1)) * 16);
		queueInitializeArrays.enqueueFillBuffer(clGradYBuf, -10000, 0, sizeof(float)*((grid->vars.kMcNx + 1)*(grid->vars.kMcNy + 1)*(grid->vars.kMcNz + 1)) * 16);
		queueInitializeArrays.enqueueFillBuffer(clGradZBuf, -10000, 0, sizeof(float)*((grid->vars.kMcNx + 1)*(grid->vars.kMcNy + 1)*(grid->vars.kMcNz + 1)) * 16);
		*/
		//queueInitializeArrays.enqueueWriteBuffer(voxelBuf, GL_TRUE, 0, sizeof(float)*((grid->vars.kMcNx + 1)*(grid->vars.kMcNy + 1)*(grid->vars.kMcNz + 1)),voxels.data());

		queueInitializeArrays.enqueueWriteBuffer(voxelBuf, GL_TRUE, 0, sizeof(float)*((grid->vars.kMcNx + 1)*(grid->vars.kMcNy + 1)*(grid->vars.kMcNz + 1)), voxels.data());
		queueInitializeArrays.enqueueWriteBuffer(voxelGradXBuf, GL_TRUE, 0, sizeof(float)*((grid->vars.kMcNx + 1)*(grid->vars.kMcNy + 1)*(grid->vars.kMcNz + 1)), voxelGradX.data());
		queueInitializeArrays.enqueueWriteBuffer(voxelGradYBuf, GL_TRUE, 0, sizeof(float)*((grid->vars.kMcNx + 1)*(grid->vars.kMcNy + 1)*(grid->vars.kMcNz + 1)), voxelGradY.data());
		queueInitializeArrays.enqueueWriteBuffer(voxelGradZBuf, GL_TRUE, 0, sizeof(float)*((grid->vars.kMcNx + 1)*(grid->vars.kMcNy + 1)*(grid->vars.kMcNz + 1)), voxelGradZ.data());

		/*queueInitializeArrays.enqueueWriteBuffer(clVertXBuf, GL_TRUE, 0, sizeof(float)*((grid->vars.kMcNx + 1)*(grid->vars.kMcNy + 1)*(grid->vars.kMcNz + 1)) * 16,clVertX.data());
		queueInitializeArrays.enqueueWriteBuffer(clVertYBuf, GL_TRUE, 0, sizeof(float)*((grid->vars.kMcNx + 1)*(grid->vars.kMcNy + 1)*(grid->vars.kMcNz + 1)) * 16,clVertY.data());
		queueInitializeArrays.enqueueWriteBuffer(clVertZBuf, GL_TRUE, 0, sizeof(float)*((grid->vars.kMcNx + 1)*(grid->vars.kMcNy + 1)*(grid->vars.kMcNz + 1)) * 16,clVertZ.data());

		queueInitializeArrays.enqueueWriteBuffer(clGradXBuf, GL_TRUE, 0, sizeof(float)*((grid->vars.kMcNx + 1)*(grid->vars.kMcNy + 1)*(grid->vars.kMcNz + 1)) * 16,clGradX.data());
		queueInitializeArrays.enqueueWriteBuffer(clGradYBuf, GL_TRUE, 0, sizeof(float)*((grid->vars.kMcNx + 1)*(grid->vars.kMcNy + 1)*(grid->vars.kMcNz + 1)) * 16,clGradY.data());
		queueInitializeArrays.enqueueWriteBuffer(clGradZBuf, GL_TRUE, 0, sizeof(float)*((grid->vars.kMcNx + 1)*(grid->vars.kMcNy + 1)*(grid->vars.kMcNz + 1)) * 16,clGradZ.data());
		*/


		/*queueInitializeArrays.enqueueReadBuffer(voxelBuf, GL_TRUE, 0, sizeof(float)*((grid->vars.kMcNx + 1)*(grid->vars.kMcNy + 1)*(grid->vars.kMcNz + 1)), voxels.data());
		queueInitializeArrays.enqueueReadBuffer(voxelGradXBuf, GL_TRUE, 0, sizeof(float)*((grid->vars.kMcNx + 1)*(grid->vars.kMcNy + 1)*(grid->vars.kMcNz + 1)), voxelGradX.data());
		queueInitializeArrays.enqueueReadBuffer(voxelGradYBuf, GL_TRUE, 0, sizeof(float)*((grid->vars.kMcNx + 1)*(grid->vars.kMcNy + 1)*(grid->vars.kMcNz + 1)), voxelGradY.data());
		queueInitializeArrays.enqueueReadBuffer(voxelGradZBuf, GL_TRUE, 0, sizeof(float)*((grid->vars.kMcNx + 1)*(grid->vars.kMcNy + 1)*(grid->vars.kMcNz + 1)), voxelGradZ.data());

		queueInitializeArrays.enqueueReadBuffer(clVertXBuf, GL_TRUE, 0, sizeof(float)*((grid->vars.kMcNx + 1)*(grid->vars.kMcNy + 1)*(grid->vars.kMcNz + 1)) * 15, clVertX.data());
		queueInitializeArrays.enqueueReadBuffer(clVertYBuf, GL_TRUE, 0, sizeof(float)*((grid->vars.kMcNx + 1)*(grid->vars.kMcNy + 1)*(grid->vars.kMcNz + 1)) * 15, clVertY.data());
		queueInitializeArrays.enqueueReadBuffer(clVertZBuf, GL_TRUE, 0, sizeof(float)*((grid->vars.kMcNx + 1)*(grid->vars.kMcNy + 1)*(grid->vars.kMcNz + 1)) * 15, clVertZ.data());

		queueInitializeArrays.enqueueReadBuffer(clGradXBuf, GL_TRUE, 0, sizeof(float)*((grid->vars.kMcNx + 1)*(grid->vars.kMcNy + 1)*(grid->vars.kMcNz + 1)) * 15, clGradX.data());
		queueInitializeArrays.enqueueReadBuffer(clGradYBuf, GL_TRUE, 0, sizeof(float)*((grid->vars.kMcNx + 1)*(grid->vars.kMcNy + 1)*(grid->vars.kMcNz + 1)) * 15, clGradY.data());
		queueInitializeArrays.enqueueReadBuffer(clGradZBuf, GL_TRUE, 0, sizeof(float)*((grid->vars.kMcNx + 1)*(grid->vars.kMcNy + 1)*(grid->vars.kMcNz + 1)) * 15, clGradZ.data());
		*/
	}

	void allocArraysCL() {


		/*clVertXBuf = cl::Buffer(*context, CL_MEM_READ_WRITE | CL_MEM_ALLOC_HOST_PTR, sizeof(float)*((grid->vars.kMcNx + 1)*(grid->vars.kMcNy + 1)*(grid->vars.kMcNz + 1)) * 16, clVertX.data());
		clVertYBuf = cl::Buffer(*context, CL_MEM_READ_WRITE | CL_MEM_ALLOC_HOST_PTR, sizeof(float)*((grid->vars.kMcNx + 1)*(grid->vars.kMcNy + 1)*(grid->vars.kMcNz + 1)) * 16, clVertY.data());
		clVertZBuf = cl::Buffer(*context, CL_MEM_READ_WRITE | CL_MEM_ALLOC_HOST_PTR, sizeof(float)*((grid->vars.kMcNx + 1)*(grid->vars.kMcNy + 1)*(grid->vars.kMcNz + 1)) * 16, clVertZ.data());

		clGradXBuf = cl::Buffer(*context, CL_MEM_READ_WRITE | CL_MEM_ALLOC_HOST_PTR, sizeof(float)*((grid->vars.kMcNx + 1)*(grid->vars.kMcNy + 1)*(grid->vars.kMcNz + 1)) * 16, clGradX.data());
		clGradYBuf = cl::Buffer(*context, CL_MEM_READ_WRITE | CL_MEM_ALLOC_HOST_PTR, sizeof(float)*((grid->vars.kMcNx + 1)*(grid->vars.kMcNy + 1)*(grid->vars.kMcNz + 1)) * 16, clGradY.data());
		clGradZBuf = cl::Buffer(*context, CL_MEM_READ_WRITE | CL_MEM_ALLOC_HOST_PTR, sizeof(float)*((grid->vars.kMcNx + 1)*(grid->vars.kMcNy + 1)*(grid->vars.kMcNz + 1)) * 16, clGradZ.data());
		*/
	}

	void generateSurfaceCL() {
		
		
		
		triangulator.clearTriangles();
		
		fill(voxels.begin(), voxels.end(), -1);
		fill(voxelGradX.begin(), voxelGradX.end(), 0);
		fill(voxelGradY.begin(), voxelGradY.end(), 0);
		fill(voxelGradZ.begin(), voxelGradZ.end(), 0);
		fill(clVertX.begin(), clVertX.end(), -10000);
		fill(clVertY.begin(), clVertY.end(), -10000);
		fill(clVertZ.begin(), clVertZ.end(), -10000);
		fill(clGradX.begin(), clGradX.end(), -10000);
		fill(clGradY.begin(), clGradY.end(), -10000);
		fill(clGradZ.begin(), clGradZ.end(), -10000);
		//fillArraysCL();
		//calcMetaballField();

		//GPU computation
		float time0 = glfwGetTime();
		//calcMetaballFieldCL();
		executeMetaballFieldKernelCL();
		float time1 = glfwGetTime();

		//cout << "Time to complete: " << (time1 - time0) << endl;

		//marchingCubesCL();
		executeMCKernelCL();

		//cout << (grid->vars.kMcNx + 1)*(grid->vars.kMcNy + 1)*(grid->vars.kMcNz + 1)*15 << endl;
	}

	void generateSurface() {
		voxels.clear();
		voxels.resize(gridNum, -1);
		triangulator.clearTriangles();

		prevoxels.clear();
		prevoxels.resize(gridNum, -1);

		voxelGradX.clear();
		voxelGradX.resize(gridNum, 0);

		voxelGradY.clear();
		voxelGradY.resize(gridNum, 0);

		voxelGradZ.clear();
		voxelGradZ.resize(gridNum, 0);

		voxelGrads.clear();
		voxelGrads.resize(gridNum, Vector3());

		float time0 = glfwGetTime();
		calcMetaballField();
		float time1 = glfwGetTime();

		//cout << "Time to complete: " << (time1 - time0) << endl;

		//GPU computation
		//calcMetaballFieldCL();


		int delta[][3] = {
			{ 0, 0, 1 },
			{ 1, 0, 1 },
			{ 1, 0, 0 },
			{ 0, 0, 0 },
			{ 0, 1, 1 },
			{ 1, 1, 1 },
			{ 1, 1, 0 },
			{ 0, 1, 0 }
		};

		isoLevel = 1;
		
		//#pragma omp parallel for
		for (int z = 0; z < grid->vars.kMcNz; ++z) {
			for (int y = 0; y < grid->vars.kMcNy; ++y) {
				for (int x = 0; x < grid->vars.kMcNx; ++x) {

					int cubeindex = 0;
					float vs[8];
					Vector3 ps[8];
					Vector3 gs[8];
					Vector3 pBase = Vector3(x, y, z);
					const float ro = 1.5;
					for (int i = 0; i < 8; i++) {
						Vector3 vVector = Vector3(pBase.x + delta[i][0], pBase.y + delta[i][1], pBase.z + delta[i][2]);
						vs[i] = voxels[offset_3d(vVector.x, vVector.y, vVector.z)];
						ps[i] = vVector*grid->vars.kMcStep - grid->vars.kMcMargin;

						//float gX = (voxels[offset_3d(vVector.x - 1, vVector.y, vVector.z)] - voxels[offset_3d(vVector.x + 1, vVector.y, vVector.z)]) / grid->vars.kMcStep;
						//float gY = (voxels[offset_3d(vVector.x, vVector.y - 1, vVector.z)] - voxels[offset_3d(vVector.x, vVector.y + 1, vVector.z)]) / grid->vars.kMcStep;
						//float gZ = (voxels[offset_3d(vVector.x, vVector.y, vVector.z - 1)] - voxels[offset_3d(vVector.x, vVector.y, vVector.z + 1)]) / grid->vars.kMcStep;
						
						//gs[i] = Vector3(gX,gY,gZ)*0.5;
						
						gs[i] = voxelGrads[offset_3d(vVector.x, vVector.y, vVector.z)];
						
						if (vs[i] < isoLevel) cubeindex |= 1 << i;
					}

					if (triangulator.kEdgeTable[cubeindex] == 0)
						continue;
					Vector3 vertlist[12];
					Vector3 gradlist[12];
					if (triangulator.kEdgeTable[cubeindex] & 1) {
						vertlist[0] = triangulator.VertexInterp1(isoLevel, ps[0], ps[1], vs[0], vs[1]);
						gradlist[0] = triangulator.VertexInterp1(isoLevel, gs[0], gs[1], vs[0], vs[1]);
					}
					if (triangulator.kEdgeTable[cubeindex] & 2) {
						vertlist[1] = triangulator.VertexInterp1(isoLevel, ps[1], ps[2], vs[1], vs[2]);
						gradlist[1] = triangulator.VertexInterp1(isoLevel, gs[1], gs[2], vs[1], vs[2]);
					}
					if (triangulator.kEdgeTable[cubeindex] & 4) {
						vertlist[2] = triangulator.VertexInterp1(isoLevel, ps[2], ps[3], vs[2], vs[3]);
						gradlist[2] = triangulator.VertexInterp1(isoLevel, gs[2], gs[3], vs[2], vs[3]);
					}
					if (triangulator.kEdgeTable[cubeindex] & 8) {
						vertlist[3] = triangulator.VertexInterp1(isoLevel, ps[3], ps[0], vs[3], vs[0]);
						gradlist[3] = triangulator.VertexInterp1(isoLevel, gs[3], gs[0], vs[3], vs[0]);
					}
					if (triangulator.kEdgeTable[cubeindex] & 16) {
						vertlist[4] = triangulator.VertexInterp1(isoLevel, ps[4], ps[5], vs[4], vs[5]);
						gradlist[4] = triangulator.VertexInterp1(isoLevel, gs[4], gs[5], vs[4], vs[5]);
					}
					if (triangulator.kEdgeTable[cubeindex] & 32) {
						vertlist[5] = triangulator.VertexInterp1(isoLevel, ps[5], ps[6], vs[5], vs[6]);
						gradlist[5] = triangulator.VertexInterp1(isoLevel, gs[5], gs[6], vs[5], vs[6]);
					}
					if (triangulator.kEdgeTable[cubeindex] & 64) {
						vertlist[6] = triangulator.VertexInterp1(isoLevel, ps[6], ps[7], vs[6], vs[7]);
						gradlist[6] = triangulator.VertexInterp1(isoLevel, gs[6], gs[7], vs[6], vs[7]);
					}
					if (triangulator.kEdgeTable[cubeindex] & 128) {
						vertlist[7] = triangulator.VertexInterp1(isoLevel, ps[7], ps[4], vs[7], vs[4]);
						gradlist[7] = triangulator.VertexInterp1(isoLevel, gs[7], gs[4], vs[7], vs[4]);
					}
					if (triangulator.kEdgeTable[cubeindex] & 256) {
						vertlist[8] = triangulator.VertexInterp1(isoLevel, ps[0], ps[4], vs[0], vs[4]);
						gradlist[8] = triangulator.VertexInterp1(isoLevel, gs[0], gs[4], vs[0], vs[4]);
					}
					if (triangulator.kEdgeTable[cubeindex] & 512) {
						vertlist[9] = triangulator.VertexInterp1(isoLevel, ps[1], ps[5], vs[1], vs[5]);
						gradlist[9] = triangulator.VertexInterp1(isoLevel, gs[1], gs[5], vs[1], vs[5]);
					}
					if (triangulator.kEdgeTable[cubeindex] & 1024) {
						vertlist[10] = triangulator.VertexInterp1(isoLevel, ps[2], ps[6], vs[2], vs[6]);
						gradlist[10] = triangulator.VertexInterp1(isoLevel, gs[2], gs[6], vs[2], vs[6]);
					}
					if (triangulator.kEdgeTable[cubeindex] & 2048) {
						vertlist[11] = triangulator.VertexInterp1(isoLevel, ps[3], ps[7], vs[3], vs[7]);
						gradlist[11] = triangulator.VertexInterp1(isoLevel, gs[3], gs[7], vs[3], vs[7]);
					}

					/* Create the triangle */
					
					const float th = 0.5;

					for (int i = 0; triangulator.kTriTable[cubeindex][i] != -1; i += 3) {
						int i0 = triangulator.kTriTable[cubeindex][i];
						int i1 = triangulator.kTriTable[cubeindex][i+1];
						int i2 = triangulator.kTriTable[cubeindex][i+2];

						Vector3 v0 = vertlist[i0];
						Vector3 v1 = vertlist[i1];
						Vector3 v2 = vertlist[i2];

						Vector3 g0 = gradlist[i0];
						Vector3 g1 = gradlist[i1];
						Vector3 g2 = gradlist[i2];

						/*g0.x = fmax(g0.x, th);
						g0.y = fmax(g0.y, th);
						g0.z = fmax(g0.z, th);

						g1.x = fmax(g1.x, th);
						g1.y = fmax(g1.y, th);
						g1.z = fmax(g1.z, th);

						g2.x = fmax(g2.x, th);
						g2.y = fmax(g2.y, th);
						g2.z = fmax(g2.z, th);*/

						int last = triangulator.mesh.m_indices.size();

						triangulator.mesh.m_indices.push_back(last);
						triangulator.mesh.m_indices.push_back(last+1);
						triangulator.mesh.m_indices.push_back(last+2);

						triangulator.mesh.m_vertices.push_back(v0);
						triangulator.mesh.m_vertices.push_back(v1);
						triangulator.mesh.m_vertices.push_back(v2);

						triangulator.mesh.m_gradients.push_back(g0);
						triangulator.mesh.m_gradients.push_back(g1);
						triangulator.mesh.m_gradients.push_back(g2);
					}

				}
			}
		}

		triangulator.polishTriangulation();
	}

	void generateGeometry() {

		voxels.clear();
		voxels.resize((grid->vars.kMcNx +1)*(grid->vars.kMcNy +1)*(grid->vars.kMcNz +1),-1);
		triangulator.clearTriangles();

		prevoxels.clear();
		prevoxels.resize((grid->vars.kMcNx + 1)*(grid->vars.kMcNy + 1)*(grid->vars.kMcNz + 1), -1);

		//default -1
		//initCornerArrays(-1);
		
		//default: 0.0
		//grid->r = 1;

		//calcDensityField();

		float time0 = glfwGetTime();
		calcMetaballField();
		float time1 = glfwGetTime();

		cout << "Time to complete: " << (time1 - time0) << endl;
		//normalizeVoxels(-0.25,0.25);
		//generateVoxels();
		
		//for (int i = 0; i < 12; i++)edge_indices[i] = 0;
		static const int delta[][3] = {
			{ 0, 0, 0 },
			{ 1, 0, 0 },
			{ 0, 1, 0 },
			{ 1, 1, 0 },
			{ 0, 0, 1 },
			{ 1, 0, 1 },
			{ 0, 1, 1 },
			{ 1, 1, 1 }
		};
		
		isoLevel = 1;// grid->particleSize * 3 - 1;// 2.5 / kPoly6;
		for (int z = 0; z < grid->vars.kMcNz;++z) {
			for (int y = 0; y < grid->vars.kMcNy; ++y) {
				for (int x = 0; x < grid->vars.kMcNx; ++x) {

					
					float vs[8];/* = {
						voxels[offset_3d(Vector3(x,   y,   z ))],
						voxels[offset_3d(Vector3( x + 1, y,   z ))],
						voxels[offset_3d(Vector3(x,   y + 1, z ))],
						voxels[offset_3d(Vector3( x + 1, y + 1, z ))],
						voxels[offset_3d(Vector3( x,   y,   z + 1))],
						voxels[offset_3d(Vector3( x + 1, y,   z + 1))],
						voxels[offset_3d(Vector3( x,   y + 1, z + 1 ))],
						voxels[offset_3d(Vector3( x + 1, y + 1, z + 1 ))],
					};*/
					Vector3 pBase = Vector3(x, y, z);// *grid->vars.kMcStep - kMcMargin;
					for (int i = 0; i < 8; i++) {
						Vector3 vVector = Vector3(pBase.x + delta[i][0], pBase.y + delta[i][1], pBase.z + delta[i][2]);
						vs[i] = voxels[offset_3d(vVector.x,vVector.y,vVector.z)];
					}

					int config_n =
						((vs[0] < isoLevel) << 0) |
						((vs[1] < isoLevel) << 1) |
						((vs[2] < isoLevel) << 2) |
						((vs[3] < isoLevel) << 3) |
						((vs[4] < isoLevel) << 4) |
						((vs[5] < isoLevel) << 5) |
						((vs[6] < isoLevel) << 6) |
						((vs[7] < isoLevel) << 7);

					if (config_n == 0 || config_n == 255)
						continue;

					
					float gX = (voxels[offset_3d(x - 1, y, z)] - voxels[offset_3d(x + 1, y, z)]) / grid->vars.kMcStep;
					float gY = (voxels[offset_3d(x, y - 1, z)] - voxels[offset_3d(x, y + 1, z)]) / grid->vars.kMcStep;
					float gZ = (voxels[offset_3d(x, y, z - 1)] - voxels[offset_3d(x, y, z + 1)]) / grid->vars.kMcStep;

					Vector3 g = Vector3(gX, gY, gZ)*0.5;


					auto do_edge = [&](int n_edge, float va, float vb, int axis, Vector3 base) {
						if ((va < isoLevel) == (vb < isoLevel))
							return;

						Vector3 v = base;// *grid->vars.kMcStep;// -kMcMargin;

						v[axis] += (va-isoLevel) / (va - vb);

						v *= grid->vars.kMcStep;
						v = v - grid->vars.kMcMargin;
						
						

						Vector3 vg = g+ (va - isoLevel) / (va - vb);

						vg *= grid->vars.kMcStep;
						vg = vg - grid->vars.kMcMargin;

						triangulator.mesh.m_gradients.push_back(vg);
						edge_indices[n_edge] = triangulator.mesh.m_vertices.size();
						triangulator.mesh.m_vertices.push_back(v);
					};

					do_edge(0, vs[0], vs[1], 0, Vector3(x, y, z));
					do_edge(1, vs[2], vs[3], 0, Vector3(x, y + 1, z));
					do_edge(2, vs[4], vs[5], 0, Vector3(x, y, z + 1));
					do_edge(3, vs[6], vs[7], 0, Vector3(x, y + 1, z + 1));

					do_edge(4, vs[0], vs[2], 1, Vector3(x, y, z));
					do_edge(5, vs[1], vs[3], 1, Vector3(x + 1, y, z));
					do_edge(6, vs[4], vs[6], 1, Vector3(x, y, z + 1));
					do_edge(7, vs[5], vs[7], 1, Vector3(x + 1, y, z + 1));

					do_edge(8, vs[0], vs[4], 2, Vector3(x, y, z));
					do_edge(9, vs[1], vs[5], 2, Vector3(x + 1, y, z));
					do_edge(10, vs[2], vs[6], 2, Vector3(x, y + 1, z));
					do_edge(11, vs[3], vs[7], 2, Vector3(x + 1, y + 1, z));

					const uint64_t config = marching_cube_tris[config_n];
					const int n_triangles = config & 0xF;
					const int n_indices = n_triangles * 3;
					const int index_base = triangulator.mesh.m_indices.size();

					
					int offset = 4;
					for (int i = 0; i < n_indices; i++) {

						const int edge = (config >> offset) & 0xF;
						triangulator.mesh.m_indices.push_back(edge_indices[edge]);
						offset += 4;
					}
					
				}
			}
		}
		
		triangulator.polishTriangulation();
	}

	float getDensityAt(Vector3 v) {

		return getDensityAt(v.x,v.y,v.z);
	}

	Vector4 Metaball(Vector3 pos, Vector3 center) const{

		Vector4 o;

		float radius = grid->particleSize;
		float radiusSq = radius*radius;

		Vector3 dist = pos - center;
		float invDistSq = 1 / (dist.dot(dist));

		Vector3 grad = -2 * radiusSq*invDistSq*invDistSq*dist;

		o.x = grad.x;
		o.y = grad.y;
		o.z = grad.z;

		o.w = radiusSq*invDistSq;

		return o;
		
	}

	Vector4 simpleMetaball(Vector3 pos, Vector3 center) {

		Vector4 o;

		float denominator = (pos.x - center.x)*(pos.x - center.x) + (pos.y - center.y)*(pos.y - center.y) + (pos.z - center.z)*(pos.z - center.z);
		o.w = grid->particleSize / denominator;

		Vector3 oxyz = ((grid->particleSize) / (denominator*denominator))*(Vector3(2 * (center.x - pos.x), 2 * (center.y - pos.y), 2 * (center.z - pos.z)));

		o.x = oxyz.x;
		o.y = oxyz.y;
		o.z = oxyz.z;

		return o;
	}

	float SoftObjects(Vector3 pos, Vector3 center) const {

		float r2 = (pos - center).lengthSq();
		float r4 = r2*r2;
		float r6 = r4*r2;

		float a = 1;
		float b = 3;

		float b6 = pow(b, 6);
		float b4 = pow(b, 4);

		return a*(1-((4*r6)/(9*b6))+((17*r4)/(9*b4))-((22*r2)/(9*b*b)));

	}

	void calcMetaballField() {


		int delta[][3] = {
			{ 0, 0, 1 },
			{ 1, 0, 1 },
			{ 1, 0, 0 },
			{ 0, 0, 0 },
			{ 0, 1, 1 },
			{ 1, 1, 1 },
			{ 1, 1, 0 },
			{ 0, 1, 0 }
		};

		int kdelta[][3] = {
			{-1,  0, 0},
			{ 0, -1, 0},
			{ 0,  0,-1},
			{ 1, 0,  1},
			{ 0, 1, 0 },
			{ 0, 0, 1 },
		};


		//#pragma omp parallel for
		for (int z = 0; z < (grid->vars.kMcNz + 1); z++) {
			for (int y = 0; y < (grid->vars.kMcNy + 1); y++) {
				for (int x = 0; x < (grid->vars.kMcNx + 1); x++) {

					Vector3 vPos = Vector3(x, y, z);
					int i = offset_3d(vPos.x,vPos.y,vPos.z);

					//works !
					/*float sum = 0;

					for (int l = 0; l < grid->nParticles; l++) {

						Vector3 pos = grid->particles[l].pos;

						if ((pos - vPos).length() > kSmoothRadius)continue;

						sum += grid->particleSize / ((vPos.x-pos.x)*(vPos.x - pos.x)+ (vPos.y - pos.y)*(vPos.y - pos.y)+ (vPos.z - pos.z)*(vPos.z - pos.z));
					}

					voxels[i] = sum;*/

					Vector4 density = calcDensityAt(x, y, z);
					
					voxels[i] = density.w;
					voxelGrads[i] = -Vector3(density.x, density.y, density.z);
				}
			}
		}

	}

	void prepareMetaballFieldKernelCL() {
		metaballKernel = cl::Kernel(program,"calcMetaballField");

		metaballKernel.setArg(0, grid->vars.kMcNx);
		metaballKernel.setArg(1, grid->vars.kMcNy);
		metaballKernel.setArg(2, grid->vars.kMcNz);
		metaballKernel.setArg(3, grid->vars.kMcStep);
		metaballKernel.setArg(4, grid->vars.kMcMargin);
		metaballKernel.setArg(5, grid->vars.kSmoothRadius);

		posXBuf = cl::Buffer(*context, CL_MEM_READ_WRITE | CL_MEM_HOST_READ_ONLY | CL_MEM_COPY_HOST_PTR, sizeof(float)*grid->posX.size(), grid->posX.data());
		posYBuf = cl::Buffer(*context, CL_MEM_READ_WRITE | CL_MEM_HOST_READ_ONLY | CL_MEM_COPY_HOST_PTR, sizeof(float)*grid->posY.size(), grid->posY.data());
		posZBuf = cl::Buffer(*context, CL_MEM_READ_WRITE | CL_MEM_HOST_READ_ONLY | CL_MEM_COPY_HOST_PTR, sizeof(float)*grid->posZ.size(), grid->posZ.data());

		metaballKernel.setArg(6, posXBuf);
		metaballKernel.setArg(7, posYBuf);
		metaballKernel.setArg(8, posZBuf);

		startsBuf = cl::Buffer(*context, CL_MEM_READ_WRITE | CL_MEM_HOST_READ_ONLY | CL_MEM_COPY_HOST_PTR, sizeof(int)*grid->starts.size(), grid->starts.data());
		endsBuf = cl::Buffer(*context, CL_MEM_READ_WRITE | CL_MEM_HOST_READ_ONLY | CL_MEM_COPY_HOST_PTR, sizeof(int)*grid->ends.size(), grid->ends.data());

		metaballKernel.setArg(9, startsBuf);
		metaballKernel.setArg(10, endsBuf);

		metaballKernel.setArg(11, grid->vars.kNX);
		metaballKernel.setArg(12, grid->vars.kNY);
		metaballKernel.setArg(13, grid->vars.kNZ);

		voxelBuf = cl::Buffer(*context, CL_MEM_READ_WRITE | CL_MEM_HOST_READ_ONLY | CL_MEM_USE_HOST_PTR, sizeof(float)*(gridNum), voxels.data());


		metaballKernel.setArg(14, voxelBuf);

		voxelGradXBuf = cl::Buffer(*context, CL_MEM_READ_WRITE | CL_MEM_HOST_READ_ONLY | CL_MEM_USE_HOST_PTR, sizeof(float)*(gridNum), voxelGradX.data());
		voxelGradYBuf = cl::Buffer(*context, CL_MEM_READ_WRITE | CL_MEM_HOST_READ_ONLY | CL_MEM_USE_HOST_PTR, sizeof(float)*(gridNum), voxelGradY.data());
		voxelGradZBuf = cl::Buffer(*context, CL_MEM_READ_WRITE | CL_MEM_HOST_READ_ONLY | CL_MEM_USE_HOST_PTR, sizeof(float)*(gridNum), voxelGradZ.data());

		metaballKernel.setArg(15, voxelGradXBuf);
		metaballKernel.setArg(16, voxelGradYBuf);
		metaballKernel.setArg(17, voxelGradZBuf);

		metaballKernel.setArg(18, grid->vars.kPolyThreshold);
	}

	void executeMetaballFieldKernelCL() {

		metaballKernel.setArg(0, grid->vars.kMcNx);
		metaballKernel.setArg(1, grid->vars.kMcNy);
		metaballKernel.setArg(2, grid->vars.kMcNz);
		metaballKernel.setArg(3, grid->vars.kMcStep);
		metaballKernel.setArg(4, grid->vars.kMcMargin);
		metaballKernel.setArg(5, grid->vars.kSmoothRadius);

		posXBuf = cl::Buffer(*context, CL_MEM_READ_WRITE | CL_MEM_HOST_READ_ONLY | CL_MEM_USE_HOST_PTR, sizeof(float)*grid->posX.size(), grid->posX.data());
		posYBuf = cl::Buffer(*context, CL_MEM_READ_WRITE | CL_MEM_HOST_READ_ONLY | CL_MEM_USE_HOST_PTR, sizeof(float)*grid->posY.size(), grid->posY.data());
		posZBuf = cl::Buffer(*context, CL_MEM_READ_WRITE | CL_MEM_HOST_READ_ONLY | CL_MEM_USE_HOST_PTR, sizeof(float)*grid->posZ.size(), grid->posZ.data());

		metaballKernel.setArg(6, posXBuf);
		metaballKernel.setArg(7, posYBuf);
		metaballKernel.setArg(8, posZBuf);

		startsBuf = cl::Buffer(*context, CL_MEM_READ_WRITE | CL_MEM_HOST_READ_ONLY | CL_MEM_USE_HOST_PTR, sizeof(int)*grid->starts.size(), grid->starts.data());
		endsBuf = cl::Buffer(*context, CL_MEM_READ_WRITE | CL_MEM_HOST_READ_ONLY | CL_MEM_USE_HOST_PTR, sizeof(int)*grid->ends.size(), grid->ends.data());

		metaballKernel.setArg(9, startsBuf);
		metaballKernel.setArg(10, endsBuf);

		metaballKernel.setArg(11, grid->vars.kNX);
		metaballKernel.setArg(12, grid->vars.kNY);
		metaballKernel.setArg(13, grid->vars.kNZ);

		voxelBuf = cl::Buffer(*context, CL_MEM_READ_WRITE | CL_MEM_HOST_READ_ONLY | CL_MEM_USE_HOST_PTR, sizeof(float)*(gridNum), voxels.data());


		metaballKernel.setArg(14, voxelBuf);

		voxelGradXBuf = cl::Buffer(*context, CL_MEM_READ_WRITE | CL_MEM_HOST_READ_ONLY | CL_MEM_USE_HOST_PTR, sizeof(float)*(gridNum), voxelGradX.data());
		voxelGradYBuf = cl::Buffer(*context, CL_MEM_READ_WRITE | CL_MEM_HOST_READ_ONLY | CL_MEM_USE_HOST_PTR, sizeof(float)*(gridNum), voxelGradY.data());
		voxelGradZBuf = cl::Buffer(*context, CL_MEM_READ_WRITE | CL_MEM_HOST_READ_ONLY | CL_MEM_USE_HOST_PTR, sizeof(float)*(gridNum), voxelGradZ.data());

		metaballKernel.setArg(15, voxelGradXBuf);
		metaballKernel.setArg(16, voxelGradYBuf);
		metaballKernel.setArg(17, voxelGradZBuf);

		metaballKernel.setArg(18, grid->vars.kPolyThreshold);

		cl::CommandQueue queueMetaball(*context, device);

		//void * add = queueMetaball.enqueueMapBuffer(voxelBuf, GL_TRUE, CL_MAP_READ | CL_MAP_WRITE, 0, sizeof(float)*((grid->vars.kMcNx + 1)*(grid->vars.kMcNy + 1)*(grid->vars.kMcNz + 1)));

		queueMetaball.enqueueNDRangeKernel(metaballKernel, cl::NullRange, cl::NDRange(grid->vars.kMcNx, grid->vars.kMcNy, grid->vars.kMcNy));

		//queueMetaball.enqueueUnmapMemObject(voxelBuf, add);

		queueMetaball.enqueueReadBuffer(voxelBuf, GL_TRUE, 0, sizeof(float)*(gridNum), voxels.data());
		queueMetaball.enqueueReadBuffer(voxelGradXBuf, GL_TRUE, 0, sizeof(float)*(gridNum), voxelGradX.data());
		queueMetaball.enqueueReadBuffer(voxelGradYBuf, GL_TRUE, 0, sizeof(float)*(gridNum), voxelGradY.data());
		queueMetaball.enqueueReadBuffer(voxelGradZBuf, GL_TRUE, 0, sizeof(float)*(gridNum), voxelGradZ.data());
	}

	void calcMetaballFieldCL() {


		cl::Kernel kernelMetaBall(program, "calcMetaballField");
		//set input args
		kernelMetaBall.setArg(0, grid->vars.kMcNx);
		kernelMetaBall.setArg(1, grid->vars.kMcNy);
		kernelMetaBall.setArg(2, grid->vars.kMcNz);
		kernelMetaBall.setArg(3, grid->vars.kMcStep);
		kernelMetaBall.setArg(4, grid->vars.kMcMargin);
		kernelMetaBall.setArg(5, grid->vars.kSmoothRadius);

		cl::Buffer posXBuf(*context, CL_MEM_READ_WRITE | CL_MEM_HOST_READ_ONLY | CL_MEM_USE_HOST_PTR, sizeof(float)*grid->posX.size(), grid->posX.data()); 
		cl::Buffer posYBuf(*context, CL_MEM_READ_WRITE | CL_MEM_HOST_READ_ONLY | CL_MEM_USE_HOST_PTR, sizeof(float)*grid->posY.size(), grid->posY.data());
		cl::Buffer posZBuf(*context, CL_MEM_READ_WRITE | CL_MEM_HOST_READ_ONLY | CL_MEM_USE_HOST_PTR, sizeof(float)*grid->posZ.size(), grid->posZ.data());

		kernelMetaBall.setArg(6, posXBuf);
		kernelMetaBall.setArg(7, posYBuf);
		kernelMetaBall.setArg(8, posZBuf);

		cl::Buffer startsBuf(*context, CL_MEM_READ_WRITE | CL_MEM_HOST_READ_ONLY | CL_MEM_USE_HOST_PTR, sizeof(int)*grid->starts.size(), grid->starts.data());
		cl::Buffer endsBuf(*context, CL_MEM_READ_WRITE | CL_MEM_HOST_READ_ONLY | CL_MEM_USE_HOST_PTR, sizeof(int)*grid->ends.size(), grid->ends.data());

		kernelMetaBall.setArg(9, startsBuf);
		kernelMetaBall.setArg(10, endsBuf);

		kernelMetaBall.setArg(11, grid->vars.kNX);
		kernelMetaBall.setArg(12, grid->vars.kNY);
		kernelMetaBall.setArg(13, grid->vars.kNZ);

		//(grid->vars.kMcNx + 1)*(grid->vars.kMcNy + 1)*(grid->vars.kMcNz + 1)

		voxelBuf = cl::Buffer(*context, CL_MEM_READ_WRITE | CL_MEM_HOST_READ_ONLY | CL_MEM_USE_HOST_PTR, sizeof(float)*(gridNum), voxels.data());
		
		
		kernelMetaBall.setArg(14, voxelBuf);

		voxelGradXBuf = cl::Buffer(*context, CL_MEM_READ_WRITE | CL_MEM_HOST_READ_ONLY | CL_MEM_USE_HOST_PTR, sizeof(float)*(gridNum), voxelGradX.data());
		voxelGradYBuf = cl::Buffer(*context, CL_MEM_READ_WRITE | CL_MEM_HOST_READ_ONLY | CL_MEM_USE_HOST_PTR, sizeof(float)*(gridNum), voxelGradY.data());
		voxelGradZBuf = cl::Buffer(*context, CL_MEM_READ_WRITE | CL_MEM_HOST_READ_ONLY | CL_MEM_USE_HOST_PTR, sizeof(float)*(gridNum), voxelGradZ.data());

		kernelMetaBall.setArg(15, voxelGradXBuf);
		kernelMetaBall.setArg(16, voxelGradYBuf);
		kernelMetaBall.setArg(17, voxelGradZBuf);

		//
		cl::CommandQueue queueMetaball(*context,device);

		//void * add = queueMetaball.enqueueMapBuffer(voxelBuf, GL_TRUE, CL_MAP_READ | CL_MAP_WRITE, 0, sizeof(float)*((grid->vars.kMcNx + 1)*(grid->vars.kMcNy + 1)*(grid->vars.kMcNz + 1)));

		queueMetaball.enqueueNDRangeKernel(kernelMetaBall,cl::NullRange, cl::NDRange(grid->vars.kMcNz,grid->vars.kMcNy,grid->vars.kMcNx));

		//queueMetaball.enqueueUnmapMemObject(voxelBuf, add);

		queueMetaball.enqueueReadBuffer(voxelBuf, GL_TRUE, 0, sizeof(float)*(gridNum), voxels.data());
		queueMetaball.enqueueReadBuffer(voxelGradXBuf, GL_TRUE, 0, sizeof(float)*(gridNum), voxelGradX.data());
		queueMetaball.enqueueReadBuffer(voxelGradYBuf, GL_TRUE, 0, sizeof(float)*(gridNum), voxelGradY.data());
		queueMetaball.enqueueReadBuffer(voxelGradZBuf, GL_TRUE, 0, sizeof(float)*(gridNum), voxelGradZ.data());
		
		//queueMetaball.finish();
		//queueMetaball.flush();

		/*for (int i = 0; i < ((grid->vars.kMcNx + 1)*(grid->vars.kMcNy + 1)*(grid->vars.kMcNz + 1));i++) {
			voxelGrads[i].x = voxelGradX[i];
			voxelGrads[i].y = voxelGradY[i];
			voxelGrads[i].z = voxelGradZ[i];
		}*/

		
	}

	float meta(float r) {

		return grid->particleSize / r;
		//return (1 - r*r)*(1 - r*r);
	}

	void calcDensityField() {
		/*scalarField->clear();
		#pragma omp parallel for
		for (int i = 0; i < ScalarField<float>::kSize; i++) {
			int x = i / ((kMcNx + 1) * (kMcNz + 1));
			int y = (i - x * ((kMcNx + 1) * (kMcNz + 1))) / (kMcNz + 1);
			int z = i - x * ((kMcNx + 1) * (kMcNz + 1)) - y * (kMcNz + 1);

			scalarField->at(x, y, z) = calcDensityAt(x, y, z);
			voxels[offset_3d(Vector3(x, y, z))] = scalarField->at(x, y, z);
			//scalarField->at(x, y, z) = calcMetaballsAt(x, y, z);
		}*/

		/*int N = voxels.size();

		
		#pragma omp parallel for
		for (int z = 0; z < (kMcNz + 1); z++) {
			for (int y = 0; y < (kMcNy+1); y++) {
				for (int x = 0; x < (kMcNx+1); x++) {

					int i = offset_3d(Vector3(x,y,z));
					voxels[i] = calcDensityAt(x,y,z);
					
				}
			}
		}*/

		Vector3 minV = Vector3(0, 0, 0);
		Vector3 maxV = Vector3(grid->vars.kMcNx, grid->vars.kMcNy, grid->vars.kMcNz);

		/*int delta[][3] = {
			{ 0, 0, 1 },
			{ 1, 0, 1 },
			{ 1, 0, 0 },
			{ 0, 0, 0 },
			{ 0, 1, 1 },
			{ 1, 1, 1 },
			{ 1, 1, 0 },
			{ 0, 1, 0 }
		};*/

		int delta[][3] = {
			{ 0, 0, 0 },
			{ 1, 0, 0 },
			{ 0, 1, 0 },
			{ 1, 1, 0 },
			{ 0, 0, 1 },
			{ 1, 0, 1 },
			{ 0, 1, 1 },
			{ 1, 1, 1 }
		};

		//Vector3 cm = getCM();
		//precompute values in/out surface
		for (int i = 0; i < grid->nParticles; i++) {

			
			Vector3 pos = grid->particles[i].pos;

			//find grid cell and adjust values of neighbourhood
			
			//Vector3 numCells = (maxV - minV)/kMcStep;
			Vector3 cell = ceilV((pos - minV) / grid->vars.kMcStep-1);
			cell.x = (int)cell.x;
			cell.y = (int)cell.y;
			cell.z = (int)cell.z;

			Vector3 ACell = (cell + Vector3(delta[0][0], delta[0][1], delta[0][2]));
			Vector3 BCell = (cell + Vector3(delta[1][0], delta[1][1], delta[1][2]));
			Vector3 CCell = (cell + Vector3(delta[2][0], delta[2][1], delta[2][2]));
			Vector3 DCell = (cell + Vector3(delta[3][0], delta[3][1], delta[3][2]));
			Vector3 ECell = (cell + Vector3(delta[4][0], delta[4][1], delta[4][2]));
			Vector3 FCell = (cell + Vector3(delta[5][0], delta[5][1], delta[5][2]));
			Vector3 GCell = (cell + Vector3(delta[6][0], delta[6][1], delta[6][2]));
			Vector3 HCell = (cell + Vector3(delta[7][0], delta[7][1], delta[7][2]));

			/*for (int l = 0; l < 8; l++) {
				Vector3 kCell = (cell + Vector3(delta[l][0], delta[l][1], delta[l][2]));

				currentCell.v[l] = 1;
				//voxels[offset_3d(kCell.x, kCell.y, kCell.z)] += implicitSurface(0);
				prevoxels[offset_3d(kCell.x, kCell.y, kCell.z)] = 1;

			}*/
			A[offset_3d(ACell.x, ACell.y, ACell.z)] = 1;
			B[offset_3d(BCell.x, BCell.y, BCell.z)] = 1;
			C[offset_3d(CCell.x, CCell.y, CCell.z)] = 1;
			D[offset_3d(DCell.x, DCell.y, DCell.z)] = 1;
			E[offset_3d(ECell.x, ECell.y, ECell.z)] = 1;
			F[offset_3d(FCell.x, FCell.y, FCell.z)] = 1;
			G[offset_3d(GCell.x, GCell.y, GCell.z)] = 1;
			H[offset_3d(HCell.x, HCell.y, HCell.z)] = 1;

		}

		//calculate values from 8 neighboring cubes
		for (int z = 0; z < (grid->vars.kMcNz+1); z++) {
			for (int y = 0; y < (grid->vars.kMcNy+1); y++) {
				for (int x = 0; x < (grid->vars.kMcNx+1); x++) {

					int i = offset_3d(Vector3(x, y, z));
					
					/*//up four group
					int i_0 = offset_3d(Vector3(x - 1, y + 1, z + 1));
					int i_1 = offset_3d(Vector3(x, y + 1, z + 1));
					int i_2 = offset_3d(Vector3(x - 1, y + 1, z));
					int i_3 = offset_3d(Vector3(x, y + 1, z));

					//down four group
					int i_4 = offset_3d(Vector3(x - 1, y, z + 1));
					int i_5 = offset_3d(Vector3(x, y, z + 1));
					int i_6 = offset_3d(Vector3(x - 1, y, z));
					int i_7 = offset_3d(Vector3(x, y, z));

					int indexes[] = {i_0,i_1,i_2,i_3,i_4,i_5,i_6,i_7};*/
					float sum = 0;

					/*for (int l = 0; l < 8; l++) {
						//sum += prevoxels[indexes[l]];

						MCCube cube = mcCubeCells[indexes[l]];
						//cout << cube.v[0] << endl;;
						sum+= cube.v[0];
					}*/

					//boundary conditions
					if (x == 0) {
						if (A[i] == 1)A[i] = -1;
						if (D[i] == 1)D[i] = -1;
						if (E[i] == 1)E[i] = -1;
						if (H[i] == 1)H[i] = -1;
					}
					if (y == 0) {
						if (A[i] == 1)A[i] = -1;
						if (B[i] == 1)B[i] = -1;
						if (C[i] == 1)C[i] = -1;
						if (D[i] == 1)D[i] = -1;
					}
					if (z == 0) {
						if (C[i] == 1)C[i] = -1;
						if (D[i] == 1)D[i] = -1;
						if (G[i] == 1)G[i] = -1;
						if (H[i] == 1)H[i] = -1;
					}

					if (x == grid->vars.kMcNx - 1) {
						if (B[i] == 1)B[i] = -1;
						if (C[i] == 1)C[i] = -1;
						if (F[i] == 1)F[i] = -1;
						if (G[i] == 1)G[i] = -1;
					}

					if (y == grid->vars.kMcNy - 1) {
						if (E[i] == 1)E[i] = -1;
						if (F[i] == 1)F[i] = -1;
						if (G[i] == 1)G[i] = -1;
						if (H[i] == 1)H[i] = -1;
					}

					if (z == grid->vars.kMcNz - 1) {
						if (A[i] == 1)A[i] = -1;
						if (B[i] == 1)B[i] = -1;
						if (F[i] == 1)F[i] = -1;
						if (E[i] == 1)E[i] = -1;
					}

					sum = A[i] + B[i] + C[i] + D[i] + E[i] + F[i] + G[i] + H[i];

					if (abs(sum) < 8) {
						//inside surface case
						voxels[i] = -1;
					}
					else {
						//in between
						voxels[i] = 1;
					}
				}
			}
		}
		
	}

	float cubeFun(float x) {
		return x*x*x;
	}

	Vector3 ceilV(Vector3 v) {

		return Vector3(ceilf(v.x),ceilf(v.y),ceilf(v.z));
	}

	float implicitSurface(float r) {

		float a = 1;
		float b = 3;


		if (r >= 0 && r <= b / 3) {
			return a*(1-3*(r*r)/(b*b));
		}
		else if (r>b/3 && r<=b) {
			return (3 / 2 * a)*(1-r/b)*(1-r/b);
		}
		else {
			return 0;
		}
	}

	Vector4 calcDensityAt(int x, int y, int z) {


		//float g = 2.0f;
		Vector4 field = Vector4(0,0,0,0);
		float density = 0.0;
		Vector3 vPos = Vector3(x, y, z) *grid->vars.kMcStep - grid->vars.kMcMargin;// +Vector3(-kMcMargin, -kMcMargin, -kMcMargin);
		int gridX = (int)(vPos.x / grid->vars.kSmoothRadius), gridY = (int)(vPos.y /grid->vars.kSmoothRadius), gridZ = (int)(vPos.z / grid->vars.kSmoothRadius);

		int count=0;
		float sum = 0;
		Vector3 denV = Vector3(0,0,0);
		for (int i = std::max(gridX - 1, 0); i <= std::min(gridX + 1, grid->vars.kNX - 1); i++) {
			for (int j = std::max(gridY - 1, 0); j <= std::min(gridY + 1, grid->vars.kNY - 1); j++) {
				for (int k = std::max(gridZ - 1, 0); k <= std::min(gridZ + 1, grid->vars.kNZ - 1); k++) {
					int gridIndex = i * grid->vars.kNY * grid->vars.kNZ + j * grid->vars.kNZ + k;
					
					for (int particle = grid->starts[gridIndex];
						//field.w <= grid->r & particle < grid->ends[gridIndex]; particle++) {
						//density<grid->r & particle < grid->ends[gridIndex]; particle++) {
						//density<kPolyThreshold & particle < grid->ends[gridIndex]; particle++) {
						particle < grid->ends[gridIndex]; particle++) {
						
						Vector3 pos = grid->particles[particle].pos;
						//assert(particle < kNumParticles);
						Vector3 r = (pos - vPos);

						if (r.length() > grid->vars.kSmoothRadius)continue;

						
						/*float sqrDiff = r.dot(r)-kSmoothRadius2;// fmax(0.f, kSmoothRadius2 - r.dot(r));
						//density += grid->particles[particle].nearDensity * sqrDiff * sqrDiff * sqrDiff;
						float deltaDensity = (grid->particles[particle].density) * sqrDiff * sqrDiff * sqrDiff; // m_j * W_{ij}
						//density += grid->particles[particle].nearDensity;

						//field += Metaball(point, grid->particles[particle].pos);
						//float den = grid->particles[particle].density;
						
						//density += grid->particles[particle].nearDensity;
						//denV += grid->particles[particle].pos;
						//density += (grid->r ) / ((grid->particles[particle].pos -point).length());
						//count++;
						//density += grid->particles[particle].nearDensity;

						//density += implicitSurface(r.length());
						density += grid->W(r, kSmoothRadius);*/

						//sum += grid->particleSize / ((vPos.x - pos.x)*(vPos.x - pos.x) + (vPos.y - pos.y)*(vPos.y - pos.y) + (vPos.z - pos.z)*(vPos.z - pos.z));
						//sum += meta(r.length());
						
						field += simpleMetaball(vPos,pos);
						density += SoftObjects(vPos,pos);
						//field.w += implicitSurface(r.length());//SoftObjects(vPos, pos);
					}
				}
			}
		}
		//density = field.w;
	
		field.w = density;

		return field;
	}

	float getDensityAt(int x, int y, int z) {
		return voxels[offset_3d(Vector3(x,y,z))];
	}

	bool compareDistances(const int &i1, const int &i2) {

		Vector3 v1 = triangulator.mesh.m_vertices[i1];
		Vector3 v2 = triangulator.mesh.m_vertices[i2];

		return ((v1-camPos).length() > (v2-camPos).length()) ? true : false;
	}

	void sortTriangles(Vector3 camPos) {
		
		this->camPos = camPos;

		//sort triangles from back to front for corrent alpha blending

		/*sort(triangulator.mesh.m_vertices.begin(), triangulator.mesh.m_vertices.end(),
			[&](const Vector3 & a, const Vector3 & b) -> bool
		{
			

			return ((a - camPos).length() <= (b - camPos).length()) ? true : false;
		});*/

		

		std::map<float, Triangle> sorted;
		int N = triangulator.mesh.m_triangles.size();
		for (int i = 0; i < N;i++) {

			Triangle t = triangulator.mesh.m_triangles[i];
			t.computePlane();
			Vector3 middlePoint = t.getCenter();

			float distance = abs(camPos.x*t.a + camPos.y*t.b + camPos.z*t.c + t.d);

			sorted[distance] = t;
		}

		triangulator.mesh.m_indices.clear();
		for (std::map<float, Triangle>::reverse_iterator it = sorted.rbegin(); it != sorted.rend(); ++it)
		{

			Triangle t = it->second;

			int i1 = t.i1;
			int i2 = t.i2;
			int i3 = t.i3;

			triangulator.mesh.m_indices.push_back(i1);
			triangulator.mesh.m_indices.push_back(i2);
			triangulator.mesh.m_indices.push_back(i3);
			
		}

	}

	void render() {

		//cout << triangulator.mesh.m_triangles.size() << endl;

		if (rendered == false) {
			//setupRenderBuffers();
			rendered = true;
		}

		

		glDeleteVertexArrays(1, &VAO);
		glDeleteVertexArrays(1, &VBO);
		glDeleteVertexArrays(1, &EBO);
		setupRenderBuffers();
		
		glActiveTexture(GL_TEXTURE3);
		glBindTexture(GL_TEXTURE_CUBE_MAP, texture);
		

		//calcTangentSpace();
		// draw mesh
		glBindVertexArray(VAO);

		//glActiveTexture(GL_TEXTURE0);
		//glBindTexture(GL_TEXTURE_2D, texture);
		//glDrawElements(GL_TRIANGLES, triangulator.mesh.m_indices.size(), GL_UNSIGNED_SHORT, (void*)0);
		//glDrawElements(GL_TRIANGLES, triangulator.mesh.m_vertices.size()/9, GL_UNSIGNED_INT, 0);
		
		//todo: instanced rendering
		//glDrawArraysInstanced(GL_TRIANGLES, 0, triangulator.mesh.m_indices.size(),0);

		//works
		glDrawArrays(GL_TRIANGLES, 0, triangulator.mesh.m_indices.size());
		
		//glDrawElements(GL_TRIANGLES, triangulator.mesh.m_indices.size(), GL_UNSIGNED_INT, 0);
		
		glBindVertexArray(0);

		
	}

	
};