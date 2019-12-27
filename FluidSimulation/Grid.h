#pragma once

#define _USE_MATH_DEFINES
#include "Particle.h"
#include "Vectors.h"
#include "Models.h"
#include "Clock.h"

#include <math.h>
#include <vector>
#include <memory>
#include <omp.h>
#include<algorithm> // for copy() and assign() 
#include<iterator> // for back_inserter 


using namespace std;

class Grid {

public :
	int width, height;

	//parameters
	//default : 0.02
	float mass = 1;// 0.02;

	//varied from 1.5 to 6
	//float volume = 1.8183930978059426;
	float particleSize = 0.035;// 0.05;// 0.035;// 0.15;// 0.035;
	int nParticles;

	//kernel radius: varied from 0.15 to 1
	float h = 0.35;

	//for water: 998.2
	float rest_density = 998.2;

	//pressure constant
	const float K = 5;

	//default: water 3.5
	float viscosity = 3.5;

	float surfaceTension = 0.0728;
	float surfaceTensionThreshold = 6;

	//for water:0
	float buoyancyConstant = 0;

	//surface threshold
	float atmPressure = 14.696;//101325;
	float gasConstant = 8.31446261815324;

	Clock clk;

	vector<vector<int>> hashmap;

	vector<Particle> particles;

	//large prime numbers
	int p1 = 73856093;
	int p2 = 19349663;
	int p3 = 83492791;
	int nH=0;

	float delay = 0.01;
	bool isSimulating = false;
	bool isRendering = false;

	float R;
	/*vector<vector<int>> mcGridIdx;
	vector<vector<int>> mcGridIdy;
	vector<vector<int>> mcGridIdz;*/
	
	vector<vector<int>> mcGridIdx;

	int nX = 0;
	int nY = 0;
	int nZ = 0;

	bool initMC = false;

	vector<float> cellDensities;

	bool simulateFirst = true;


public:

	Grid() {

		clk = Clock();

	}

	Grid(int width,int height) {
		this->width = width;
		this->height = height;

		//nParticles = (int) (volume / (4 * M_PI*pow(particleSize, 3)));

		float particlePerDim = 5;
		nParticles = pow(particlePerDim, 3);

		clk = Clock();
		nH = clk.NextPrime(2*nParticles);

		//initialize hash table
		for (int i = 0; i < nH; i++) {
			
			vector<int> emptyRow;
			emptyRow.push_back(0);
			//hashmap.push_back(indexes);
			hashmap.push_back(emptyRow);
		}
		clearHashMap();
		
		//radius : 1.0

		//cout << nParticles << endl;

		int cnt = 0;
		float length = particlePerDim*particleSize;
		//default gap: 0.005
		float gap = particleSize;//0.09;//0.015;
		length += gap*particlePerDim;
		#pragma omp parallel for
		for (float ix = 0; ix < length; ix+=(particleSize+gap)) {
			for (float j = 0; j < length; j += (particleSize+gap)) {
				for (float k = 0; k < length; k += (particleSize+gap)) {

					Particle pIJK = Particle(ix, j, k,cnt);

					particles.push_back(pIJK);
					cnt++;
				}
			}
		}
		//cout << particles.size() << endl;

		//for (int i = 0; i < 10000;i++)simulate();
		nParticles = cnt;

		
		//for (int i = 0; i < 1; i++)simulate();
	}


	/***
	 Kriging Method 3D
	 ***/
	void attachParticlesToGrid(Vector3 min, Vector3 size, Vector3 cubeSize) {

		
		
		nX = size.x / cubeSize.x;
		nY = size.y / cubeSize.y;
		nZ = size.z / cubeSize.z;


		int size3D = nX*nY*nZ;

		int neighbours = 10;
		cellDensities.clear();
		#pragma omp parallel for
		for (int i = 0; i < size3D;i++) {

			int iX = i / (nY*nZ);
			int iY = (i-iX*nY*nZ) / nZ;
			int iZ = i-iX*nY*nZ-iY*nZ;

			float cellX = min.x + iX*cubeSize.x;
			float cellY = min.y + iY*cubeSize.y;
			float cellZ = min.z + iZ*cubeSize.z;

			Vector3 cellPos = Vector3(cellX,cellY,cellZ);

			vector<float> distances;

			for (int j = 0; j < nParticles; j++) {
				Vector3 pos = particles[j].pos;

				float dist = (cellPos - pos).length();

				distances.push_back(dist);
			}

			sort(distances.begin(), distances.end());

			vector<float> nDistances;
			vector<float> nDensities;
			for (int k = 0; k < neighbours; k++) {

				for (int j = 0; j < nParticles; j++) {

					Vector3 pos = particles[j].pos;

					float dist = (cellPos - pos).length();

					if (distances[k] == dist) {
						nDistances.push_back(dist);
						nDensities.push_back(particles[j].density);
						break;
					}
				}
			}

			float distSum = 0;
			for (int j = 0; j < nDistances.size();j++) {
				distSum += nDistances[j];
			}

			float sumDensity = 0;

			//inverse weighting function
			float exponent = 0.002;
			float weightSum = 0;
			for (int j = 0; j < neighbours;j++) {
				float dist = nDistances[j];

				float weight = 1. / pow(dist,exponent);
				weightSum += weight;
			}

			for (int j = 0; j < neighbours;j++) {
				float dist = nDistances[j];

				float weight = 1. / pow(dist, exponent);

				sumDensity += weight*nDensities[j];
			}

			sumDensity /= weightSum;
			cellDensities.push_back(sumDensity);
		}
		
	}

	void normalizeCells() {

		float minValue = FLT_MAX;
		float maxValue = -FLT_MAX;

		int N = cellDensities.size();
		for (int i = 0; i < N;i++) {
			float v = cellDensities[i];
			if (v < minValue)minValue = v;
			if (v > maxValue)maxValue = v;
		}

		//normalize to 0-1
		float newMax = 1;
		float newMin = 0;

		for (int i = 0; i < N; i++) {
			float v = cellDensities[i];
			
			cellDensities[i] = (v - minValue)*(newMax - newMin) / (maxValue-minValue) + newMin;
		}
	}

	float getCellDensity(int x, int y, int z) {

		if (x < 0)x = 0;
		if (y < 0)y = 0;
		if (z < 0)z = 0;

		if (x > nX - 1)x = nX - 1;
		if (y > nY - 1)y = nY - 1;
		if (z > nZ - 1)z = nZ - 1;

		return cellDensities[x*nY*nZ + y*nZ + z];
	}

	void showme(Vector3 v, string s) {

		cout << s << "|" <<v.x << "|" << v.y << "|" << v.z << endl;
	}

	void clearMCGrid() {
		

		for (int i = 0; i < mcGridIdx.size();i++) {
			mcGridIdx[i].clear();
		}
		mcGridIdx.clear();

	}

	void simulate() {

		//if (isRendering)return;
		isSimulating = true;
		
		
		#pragma omp parallel for
		for (int i = 0; i < nParticles;i++) {
			//particles[i].isRendering = false;
			updatePosition(i);
			addToHashMap(i);
		}

		//int cnt = 0;
		#pragma omp parallel for
		for (int i = 0; i < nParticles; i++) {
			
			getNeighbours(i);

		}

		
		//cout << "num of cnts: " << cnt << endl;
		#pragma omp parallel for
		for (int i = 0; i < nParticles; i++) {

			particles[i].density =  getParticleDensity(i);
			particles[i].pressure = getParticlePressure(i);

		}
		
		//default: 0.2
		float damp = 1;// 0.2;

		float damp2 = 0.2;

		#pragma omp parallel for
		for (int i = 0; i < nParticles; i++) {

			
			
			Vector3 internalForces = getPressureForce(i) + getViscosityForce(i);
			Vector3 externalForces = getGravityForce(i)*damp + getSurfaceTension(i) + getBuoyancyForce(i);

			particles[i].totalForce = internalForces + externalForces;
			
			particles[i].totalForce *= damp2;

			//cout << "(" << particles[i].totalForce.x << "," << particles[i].totalForce.y << "," << particles[i].totalForce.z << ")" << endl;
		}
		
		#pragma omp parallel for
		for (int i = 0; i < nParticles;i++) {
			
			updateAcceleration(i);
			updateVelocity(i);
			XSPH(i);
			
			//cout << "accel" << particles[i].acceleration<<endl;
			
			if (particles[i].pos.y < (-2.5+particleSize/2)) {

				particles[i].pos.y = (-2.5 + particleSize / 2);
				Vector3 n = Vector3(0,1,0);
				
				float k = n.dot(particles[i].vel);
				//particles[i].vel -= 2 * k*n*clk.deltaTime*10;
				
				particles[i].totalForce -= getGravityForce(i)*damp;
				

				//particles[i].vel *= delay;
				//particles[i].vel = delay*particles[i].vel;

			}

			//updateVelocity(i);
		}
		
		clearHashMap();
		
		isSimulating = false;

	}

	float getSurfaceThreshold() {

		atmPressure = K*rest_density;
		return atmPressure / gasConstant;
	}

	void setTriangleData(vector<Triangle> triangles) {

		for (int i = 0; i < nParticles;i++) {

			particles[i].setData(triangles);
		}
	}

	void setData(vector<Vector3> vertices, vector<int> indices) {
		for (int i = 0; i < nParticles; i++) {

			particles[i].setData(vertices,indices);
		}
	}

	void setData(vector<Vector3> vertices, vector<int> indices, vector<Models::Vertex> vertexes) {
		for (int i = 0; i < nParticles; i++) {

			particles[i].setData(vertices, indices,vertexes);
		}
	}

	void setupParticles() {
		for (int i = 0; i < nParticles; i++) {

			particles[i].setupBuffers();
		}
	}

	float Wpoly6(Vector3 r, float h) {

		float rL = r.length();
		float result = 0;
		if (rL >= 0 && rL <= h) {

			float term1 = 315.f / (64.f*M_PI*pow(h,9.f));
			float term2 = pow(h*h-rL*rL,3);
			result = term1*term2;
		}
		return result;
	}

	Vector3 gradWpoly6(Vector3 r, float h) {

		float rL = r.length();
		float term1 = -945.f / (32.f*M_PI*pow(h,9.f));
		Vector3 term2 = r * (pow(h*h-rL*rL, 2.f));

		return term1 * term2;
	}

	float lapWpoly6(Vector3 r, float h) {

		float rL = r.length();
		float term1 = -945.f / (32.f * M_PI*pow(h, 9.f));
		float term2 = (h*h - rL * rL)*(3*h*h-7*rL*rL);

		return term1 * term2;
	}

	float Wspiky(Vector3 r, float h) {

		
		float rL = r.length();
		float result = 0;

		if (rL >= 0 && rL <= h) {
			float term1 = 15.f / (M_PI*pow(h, 6.f));
			float term2 = pow(h-rL,3);
			return term1 * term2;
		}

		return result;
	}

	Vector3 gradWspiky(Vector3 r, float h) {
		float term1 = -45.f / (M_PI*pow(h,6.f));
		float rL = r.length();
		//if (rL == 0)rL=1;
		Vector3 term2 = r * (pow(h-rL,2))/rL;

		return term1*term2;
	}

	float lapWSpiky(Vector3 r, float h) {

		float rL = r.length();
		float term1 = -90.f / (M_PI*pow(h,6.f));
		float term2 = (h-rL)*(h-2*rL);

		return term1 * term2;
	}

	float Wviscosity(Vector3 r, float h) {

		float rL = r.length();

		//if (rL == 0)return 0;
		float result = 0;

		if (rL >= 0 && rL <= h) {
			float term1 = 15 / (2*M_PI*pow(h, 3));
			float term2 = -(pow(rL,3)/(2*pow(h,3)))+(rL*rL)/(h*h)+h/(2*rL)-1;
			return term1 * term2;
		}

		return result;
	}

	Vector3 gradWviscosity(Vector3 r, float h) {

		float rL = r.length();

		float term1 = 15 / (2 * M_PI*pow(h, 3));
		Vector3 term2 = r * (-3*rL/(2*pow(h,3))+2/(h*h)-h/(2*pow(rL,3)));

		return term1 * term2;
	}

	float lapWviscosity(Vector3 r, float h) {

		float rL = r.length();

		float term1 = 45.f / (M_PI*pow(h, 6.f));
		float term2 = h - rL;

		return term1 * term2;
	}

	float getParticleDensity(int index) {
		
		float sum = 0;

		Vector3 rI = particles[index].pos;
		
		Particle pI = particles[index];

		int count = pI.n_indices.size();
		
		for (int i = 0; i < count;i++) {

			if (isNeighbour(particles[pI.n_indices[i]], index) == false)continue;

			Vector3 rJ = particles[pI.n_indices[i]].pos;
			sum += mass * Wpoly6(rI-rJ,h);

		}
		

		sum += rest_density;

		return sum;
	}

	float getCornerDensity(vector<int> pIdx, Vector3 pos) {

		float sum = 0;

		for (int i = 0; i < pIdx.size(); i++) {
			Vector3 rI = particles[pIdx[i]].pos;

			sum += mass/particles[pIdx[i]].density* Wpoly6(pos-rI,h);

		}

		sum /= (float)pIdx.size();
		
		sum += rest_density;
		
		/*float ri = R * 3 / 2;
		vector<float> wi_s;

		float den = 0;
		for (int i = 0; i < pIdx.size(); i++) {
			Vector3 xi = particles[pIdx[i]].pos;
			den += k2(((pos - xi).dot(pos - xi)) / (R*R));
			
		}

		if (den == 0)den = 1;

		for (int i = 0; i < pIdx.size();i++) {

			float wi = 0;
			Vector3 xi = particles[pIdx[i]].pos;
			
			float nom = k2(((pos - xi).dot(pos - xi)) / (R*R));

			wi = nom / den;
			wi_s.push_back(wi);
		}



		Vector3 term1 = Vector3(0,0,0);
		float term2 = 0;
		for (int i = 0; i < pIdx.size(); i++) {
			Vector3 xi = particles[pIdx[i]].pos;

			term1 += wi_s[i] * xi;

			term2 += wi_s[i] * ri;
		}

		sum = (pos - term1).length() - term2+rest_density;
		//cout << sum-1 << endl;
		*/
		return sum;
	}

	float k2(float s) {

		return fmax(0, pow(1 - s, 3));
	}

	float getParticlePressure(int index) {

		float Po = K * rest_density;
		
		float density = particles[index].density;

		if (density == 0)density = rest_density;
		float volume_per_unit = 1. / density;
		return K/volume_per_unit - Po;
	}

	//internal force
	Vector3 getPressureForce(int index) {

		Vector3 sumForce = Vector3(0,0,0);

		Vector3 rI = particles[index].pos;
		float pI = particles[index].pressure;
		
		Particle p_I = particles[index];

		int count = p_I.n_indices.size();

		for (int i = 0; i < count;i++) {

			if (isNeighbour(particles[p_I.n_indices[i]], index) == false)continue;

			Vector3 rJ = particles[p_I.n_indices[i]].pos;
			float pJ = particles[p_I.n_indices[i]].pressure;
			float dJ = particles[p_I.n_indices[i]].density;

			if (dJ == 0)dJ = rest_density;

			sumForce = sumForce+((pI+pJ)/2.)*(mass/dJ)*gradWspiky(rI-rJ,h);
		}

		sumForce = -sumForce;

		return sumForce;
	}

	//internal force
	Vector3 getViscosityForce(int index) {

		Vector3 sumForce = Vector3(0, 0, 0);

		

		Vector3 rI = particles[index].pos;
		Vector3 velI = particles[index].vel;

		Particle pI = particles[index];

		int count = pI.n_indices.size();
		for (int i = 0; i < count; i++) {

			if (isNeighbour(particles[pI.n_indices[i]], index) == false)continue;

			Vector3 rJ = particles[pI.n_indices[i]].pos;
			//float pJ = particles[pI.n_indices[i]].pressure;
			float dJ = particles[pI.n_indices[i]].density;
			Vector3 velJ = particles[pI.n_indices[i]].vel;

			if (dJ == 0)dJ = rest_density;

			sumForce = sumForce + (velJ-velI)*((mass/dJ)*lapWviscosity(rI-rJ,h));
		}

		sumForce = sumForce * viscosity;

		return sumForce;
	}
	
	//external force
	Vector3 getGravityForce(int index) {

		Vector3 gravity = Vector3(0,-9.8,0);

		float density = particles[index].density;

		return density * gravity;
		
	}

	//external force
	Vector3 getSurfaceTension(int index) {


		float sum = 0;
		float curvature = 0;

		//colour field
		Vector3 rI = particles[index].pos;
		
		Particle pI = particles[index];

		int count = pI.n_indices.size();
		for (int i = 0; i < count; i++) {

			if (isNeighbour(particles[pI.n_indices[i]], index) == false)continue;

			Vector3 rJ = particles[pI.n_indices[i]].pos;
			float dJ = particles[pI.n_indices[i]].density;

			if (dJ == 0)dJ = rest_density;

			sum += (mass / dJ)*Wpoly6(rI-rJ,h);
			curvature += (mass / dJ)*lapWpoly6(rI - rJ, h);
		}


		Vector3 sumForce = Vector3(0, 0, 0);

		Vector3 n = Vector3(0, 0, 0);
		
		for (int i = 0; i < count; i++) {

			if (isNeighbour(particles[pI.n_indices[i]], index) == false)continue;

			Vector3 rJ = particles[pI.n_indices[i]].pos;
			float dJ = particles[pI.n_indices[i]].density;

			if (dJ == 0)dJ = rest_density;
			n += (mass / dJ)*gradWpoly6(rI - rJ, h);
		}

		if (n.length() < surfaceTensionThreshold)return Vector3(0,0,0);

		if(n.length()!=0)curvature = -curvature/n.length();

		sumForce = surfaceTension*curvature * n;

		

		return sumForce;
	}

	//external force
	Vector3 getBuoyancyForce(int index) {

		Vector3 gravity = Vector3(0,-9.8,0);

		Vector3 sumForce = buoyancyConstant * (particles[index].density - rest_density)*gravity;

		return sumForce;
	}

	//TODO force by interaction

	void updateAcceleration(int index) {

		Particle pI = particles[index];

		if (pI.density != 0)particles[index].acceleration = pI.totalForce / pI.density;
		else particles[index].acceleration = pI.totalForce/rest_density;
		//particles[index].acceleration = pI.totalForce / mass;

	}

	void updateVelocity(int index) {
		
		Particle pI = particles[index];
		particles[index].vel = pI.vel + clk.deltaTime*pI.acceleration;
		
		//Vector3 n = Vector3(0, 1, 0);
		//float distance = particles[index].pos.y - (-2.5 + particleSize / 2);
		//float restitution = 0.5*distance / (clk.deltaTime*particles[index].vel.length());
		 
		//particles[index].vel = particles[index].vel - (1.0 + restitution)*(particles[index].vel.dot(n))*n;

	}

	void updatePosition(int index) {

		Particle pI = particles[index];
		particles[index].pos = pI.pos + clk.deltaTime*pI.vel;

		
	}

	void XSPH(int index) {

		float e = 0.1;

		Vector3 velI = particles[index].vel;
		float denI = particles[index].density;
		Vector3 posI = particles[index].pos;
		Vector3 sumVector = Vector3(0,0,0);

		Particle pI = particles[index];

		int count = pI.n_indices.size();
		for (int i = 0; i < count; i++) {

			if (isNeighbour(particles[pI.n_indices[i]], index) == false)continue;

			Vector3 velJ = particles[pI.n_indices[i]].vel;
			float denJ = particles[pI.n_indices[i]].density;
			Vector3 posJ = particles[pI.n_indices[i]].pos;
			
			sumVector = sumVector + (2 * mass) / (denI + denJ)*Wpoly6(posI-posJ,h);

		}
		velI = velI + e*sumVector;

		particles[index].vel = velI;
	}

	int hash(Vector3 pos) {

		Vector3 discX = discretizeX(pos);
		return ((((int)discX.x*p1) ^ ((int)discX.y*p2) ^ ((int)discX.z*p3)));
		//return (((int)discX.x*p1) ^ ((int)discX.y*p2) ^ ((int)discX.z*p3)) % nParticles;
	}

	Vector3 discretizeX(Vector3 floatVector) {

		return Vector3(((int)floor(floatVector.x/h)), (int)(floor(floatVector.y / h)), (int)(floor(floatVector.z / h)));
	}

	void addToHashMap(int index) {

		//hashmap[hash(particles[index].pos)] = index;
		
		int key = hash(particles[index].pos);

		if (key < 0) {
			key = (nH - key)%nH;
		}
		else key %= nH;

		hashmap[key].push_back((int)(particles[index].id));

	}

	void getNeighbours(int index) {

		
		particles[index].clearNeighbourIndices();
	
		/*int commonKey = hash(particles[index].pos);

		particles[index].count = 0;
		
		for (int i = 0; i < nParticles;i++) {
			if (i == index)continue;
			//Particle* p = &particles[i];
			int key = hash(particles[i].pos);
			if (key == commonKey) {
				if (isNeighbour(particles[i], index)) {
					particles[index].count++;
					//particles[index].neighbours.push_back(particlePointers[i]);
					particles[index].n_indices.push_back(i);
				}
			}
		}*/

		/*int size = hashmap[hash(particles[index].pos)].size();
		for (int i = 0; i < size;i++) {
			particles[index].n_indices.push_back(hashmap[hash(particles[index].pos)][i]);
		}*/
		
		int key = hash(particles[index].pos);
		if (key < 0) {
			key = (nH - key) % nH;
		}
		else key %= nH;

		particles[index].n_indices = hashmap[key];
		//copy(hashmap[key].begin(), hashmap[key].end(), back_inserter(particles[index].n_indices));
	}

	void clearHashMap() {
		for (int i = 0; i < nH;i++) {
			hashmap[i].clear();
			//hashmap[i] = -1;
		}
	}

	void terminate() {
		
		
	}

	//check if neighbours from hashmap
	//list are true neighbours
	/*bool isNeighbour(Particle pI, Particle pJ) {

		if ((pI.pos - pJ.pos).length() <= h)return true;

		return false;
	}*/

	bool isNeighbour(Particle p, int index) {

		if (p.id == index)return false;
		if ((particles[index].pos-p.pos).length() <= (2*h))return true;

		return false;
	}


	
};