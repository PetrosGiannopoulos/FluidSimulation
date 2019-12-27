#pragma once

#define _USE_MATH_DEFINES
#include "Particle.h"
#include "Vectors.h"
#include "Models.h"
#include "Clock.h"
//#include "constant.h"
#include <math.h>
#include <vector>
#include <memory>
#include <random>
#include <array>
#include <omp.h>
#include<algorithm> // for copy() and assign() 
#include<iterator> // for back_inserter 

#include "ConstVars.h"

using namespace std;

class SphGrid {

public:

	//std::array<int, kGridSize> start, end;// , prevStart, prevEnd;
	
	vector<int> start, end;
	
	//std::array<int, kNumParticles> next, prevNext;

	vector<int> next, prevNext;

	vector<int> starts, ends, starts_, ends_;

	/*std::array<int, kGridSize> starts;
	std::array<int, kGridSize> ends;
	std::array<int, kGridSize> starts_;
	std::array<int, kGridSize> ends_;*/
	//std::array<int, kNumParticles> nexts_;
	//std::array<Vector3, kNumParticles> positions_;
	//std::array<Vector3, kNumParticles> velocities_;

	vector<int> nexts_;
	vector<Vector3> positions_;
	vector<Vector3> velocities_;

	vector<float> posX;
	vector<float> posY;
	vector<float> posZ;

	Clock clk;

	vector<Particle> particles;
	vector<Vector3> offsets;

	//varied from 1.5 to 6
	//float volume = 1.8183930978059426;
	float particleSize = 0.5;// 0.5;
	int nParticles;

	//large prime numbers
	int p1 = 73856093;
	int p2 = 19349663;
	int p3 = 83492791;
	int nH = 0;

	float delay = 0.01;
	bool isSimulating = false;
	bool isRendering = false;

	bool simulateFirst = true;

	Vector3 axis, right;
	float externalForce;
	glm::mat3 reverseRotation;

	Vector3 bbmin, bbmax;


	//metaballs threshold;
	float r = 1;
	float multiplier = 1;
	float sigma;
	float simTime = 0;

	ConstVars vars;

	bool firstSim = true;
public:

	template<class F>
	void forEachNeighbor(int particle, F f) {
		Vector3 _index = particles[particle].pos / vars.kSmoothRadius;
		Vector3 &index = _index;
		#pragma omp parallel for
		for (int x = std::max(0, ((int)index.x) - 1);
			x <= std::min(vars.kNX - 1,((int)index.x) + 1); x++) {

			for (int y = std::max(0, ((int)index.y) - 1);
				y <= std::min(vars.kNY - 1, ((int)index.y) + 1); y++) {

				for (int z = std::max(0, ((int)index.z) - 1);
					z <= std::min(vars.kNZ - 1, ((int)index.z) + 1); z++) {

					for (int neighbor = starts[x * vars.kNY * vars.kNZ + y * vars.kNZ + z];
						neighbor < ends[x * vars.kNY * vars.kNZ + y * vars.kNY + z]; neighbor++) {

						//if(particle<neighbor){
						if (particle <= neighbor) {
							f(neighbor);
						}
					}
				}
			}
		}
	}

	SphGrid() {

		clk = Clock();
		clear();

		vars = ConstVars();
	}

	SphGrid(int width, int height) {

		vars = ConstVars();
		clk = Clock();
		clear();

		setReverseRotation(glm::mat3(1.f));

		auto gen = std::default_random_engine();
		auto distX = std::uniform_real_distribution<float>(vars.kXMax*0.7, vars.kXMax),
			distY = std::uniform_real_distribution<float>(vars.kYMax),
			distZ = std::uniform_real_distribution<float>(vars.kZMax);
		
		
		nParticles = vars.kNumParticles;
		int nCount = 0;
		//#pragma omp parallel for
		for (int i = 0; i < vars.kNumParticles; i++) {

			Particle p = Particle();
			p.pos = Vector3(distX(gen),distY(gen),distZ(gen));

			p.vel = Vector3(0, 0, 0);
			
			particles.push_back(p);
			nCount++;
			offsets.push_back(Vector3(0,0,0));

			posX.push_back(0);
			posY.push_back(0);
			posZ.push_back(0);
		}
		//nParticles = particles.size();
		//vars.kNumParticles = nParticles;
		
		updateGrid();

		sigma = 315.0 / 64.0 / vars.kPi;
		
	}

	void recalculate() {
		particles.clear();
		posX.clear();
		posY.clear();
		posZ.clear();
		offsets.clear();

		clear();

		setReverseRotation(glm::mat3(1.f));

		auto gen = std::default_random_engine();
		auto distX = std::uniform_real_distribution<float>(vars.kXMax*0.7, vars.kXMax),
			distY = std::uniform_real_distribution<float>(vars.kYMax),
			distZ = std::uniform_real_distribution<float>(vars.kZMax);


		nParticles = vars.kNumParticles;
		int nCount = 0;
		//#pragma omp parallel for
		for (int i = 0; i < vars.kNumParticles; i++) {

			Particle p = Particle();
			p.pos = Vector3(distX(gen), distY(gen), distZ(gen));

			p.vel = Vector3(0, 0, 0);

			particles.push_back(p);
			nCount++;
			offsets.push_back(Vector3(0, 0, 0));

			posX.push_back(0);
			posY.push_back(0);
			posZ.push_back(0);
		}
		//nParticles = particles.size();
		//vars.kNumParticles = nParticles;

		updateGrid();

		sigma = 315.0 / 64.0 / vars.kPi;
	}

	Vector3 getRandomSet(vector<float> &discrete) {

		int indX = rand()%discrete.size();
		int indY = rand() % discrete.size();
		int indZ = rand() % discrete.size();
		/*auto gen = std::default_random_engine();
		auto indexX = std::uniform_real_distribution<int>(discrete.size());
		auto indexY = std::uniform_real_distribution<int>(discrete.size());
		auto indexZ = std::uniform_real_distribution<int>(discrete.size());*/

		return Vector3(discrete[indX], discrete[indY], discrete[indZ]);
	}

	void simulate() {
		//#pragma omp parallel for
		for (int i = 0; i < nParticles; i++) {
			particles[i].density = 0;
			particles[i].field = Vector4(0,0,0,0);
			particles[i].acceleration = Vector3(0, 0, 0);
			particles[i].colourField = 0;
			
			//test
			//particles[i].vel = Vector3(0,0,0);
		}

		

		updateGrid();
		updateDensities();
		
		
		//updateAccelerations();
		updateAccelerations2();
		
		updateParticles();
		//updateParticlesRK4();
		
		updateXSPH();

	}


	
	void updateGrid() {

		clear();

		//#pragma omp parallel for
		for (int i = 0; i < vars.kNumParticles; i++) {
			auto &position = particles[i].pos;
			int g = (int)(position.x / vars.kSmoothRadius) * vars.kNY * vars.kNZ +
				(int)(position.y / vars.kSmoothRadius) * vars.kNZ + (int)(position.z / vars.kSmoothRadius);
			if (starts_[g] < 0) {
				starts_[g] = i;
				ends_[g] = i;
			}
			else {
				nexts_[ends_[g]] = i;
				ends_[g] = i;
			}
			nexts_[i] = -1;
		}

		int c = 0;
		//#pragma omp parallel for
		for (int g = 0; g < vars.kGridSize; ++g) {
			starts[g] = c;
			for (int p = starts_[g]; p >= 0; p = nexts_[p]) {
				positions_[c] = particles[p].pos;
				velocities_[c] = particles[p].vel;
				++c;
			}
			ends[g] = c;
		}

		#pragma omp parallel for
		for (int i = 0; i < vars.kNumParticles;i++) {
			particles[i].vel = velocities_[i];
			particles[i].pos = positions_[i];

			posX[i] = positions_[i].x;
			posY[i] = positions_[i].y;
			posZ[i] = positions_[i].z;
		}
		
	}

	void updateDensities() {
		//#pragma omp parallel for
		for (int i = 0; i < vars.kNumParticles; i++) {
			auto &particlePosition = particles[i].pos;
			forEachNeighbor(i, [&](int neighbor) {
				auto &neighborPosition = particles[neighbor].pos;
				auto r = particlePosition - neighborPosition;
				auto rl = r.length();
				auto dist2 = r.lengthSq();// r.dot(r); // dist = ||r||^2
				auto sqrDiff = vars.kSmoothRadius2 - dist2; // h^2 - r^2
				
				if (rl >= 0 && rl<=vars.kSmoothRadius) {
				//if (sqrDiff > 0) {
					auto deltaDensity = (vars.kParticleMass * vars.kPoly6) * sqrDiff * sqrDiff * sqrDiff; // m_j * W_{ij}
					
					particles[i].density += deltaDensity;
					particles[neighbor].density += deltaDensity;

				}
				
			});
		}

	}

	float W(Vector3 r, float h) {

		float rH = (r / h).length();

		if (rH > 0 && rH <= 1) {
			return (sigma / (pow(h, 3)))*(pow(2 - rH, 3) - 4 * pow(1 - rH, 3));
		}
		else if (rH>1 && rH <= 2) {
			return (sigma / (pow(h, 3)))*(pow(2 - rH, 3));
		}
		else {
			return 0;
		}
	}

	float Wpoly6(Vector3 r, float h) {

		float rL = r.length();
		float result = 0;
		if (rL<h) {

			float term1 = 315.f / (64.f*M_PI*pow(h, 9.f));
			float term2 = pow(h*h - rL*rL, 3);
			result = term1*term2;
		}

		//if (result != 0)cout << result*kSmoothRadius<< endl;
		return result;
	}

	void updateAccelerations() {
		//#pragma omp parallel for
		for (int i = 0; i < vars.kNumParticles; i++) {
			float particleDensity = particles[i].density;

			auto direction_ = reverseRotation*glm::normalize(glm(particles[i].pos) - glm::vec3(10.f));
			Vector3 direction = Vector3(direction_.x, direction_.y, direction_.z);
			particles[i].acceleration += axis * (vars.kGravity[1] + externalForce) +direction * 0.5f*externalForce;

			auto particlePosition = particles[i].pos;
			float particleDensityRecip = vars.kParticleMass / particleDensity;

			forEachNeighbor(i, [&](int neighbor) {

				
				auto r = particlePosition - particles[neighbor].pos;
				auto dist2 = r.lengthSq();// r.dot(r); // ||r||^2
				if (dist2 <= vars.kSmoothRadius2) {
					float dist = r.length(); // ||r||
					auto neighborDensity = particles[neighbor].density;
					auto neighborDensityRecip = vars.kParticleMass / neighborDensity;
					auto sqrDiff = vars.kSmoothRadius2 - dist2; // h^2 - r^2
					auto diff = vars.kSmoothRadius - dist; // h - r

					auto accPressure = dist == 0.f ? Vector3(0,0,0) :
						vars.kIdealGasCoeff * vars.kParticleMass *
						(particleDensity + neighborDensity -
							2 * vars.kStandardDensity) * vars.kSpiky * diff * diff /
						2.f * r / dist / particleDensity / neighborDensity;

					particles[i].acceleration -= accPressure;
					//particles[neighbor].acceleration += accPressure;

					auto accViscosity =
						vars.kViscosity * (particles[neighbor].vel - particles[i].vel)
						* vars.kLaplacianViscosity * diff * vars.kParticleMass;

					particles[i].acceleration += accViscosity / particleDensity;
					particles[neighbor].acceleration -= accViscosity / neighborDensity;

				}
			});

		}
		
		externalForce *= 0.3;
		
	}

	void updateAccelerations2() {
		//#pragma omp parallel for
		for (int i = 0; i < vars.kNumParticles; i++) {
			float particleDensity = particles[i].density;

			auto direction_ = reverseRotation*glm::normalize(glm(particles[i].pos) - glm::vec3(10.f));
			Vector3 direction = Vector3(direction_.x, direction_.y, direction_.z);
			particles[i].acceleration += axis * (vars.kGravity[1] + externalForce) + direction * 0.5f*externalForce;

			//particles[i].acceleration += (vars.kGravity + externalForce) + direction * 0.5f*externalForce;
			//particles[i].acceleration += vars.kGravity*particles[i].density;// +direction*0.5f*externalForce;

			auto particlePosition = particles[i].pos;
			float particleDensityRecip = vars.kParticleMass / particleDensity;

			forEachNeighbor(i, [&](int neighbor) {


				auto r = particlePosition - particles[neighbor].pos;
				auto dist2 = r.dot(r);// r.lengthSq();// r.dot(r); // ||r||^2
				float dist = r.length(); // ||r||
				if (dist2 <= vars.kSmoothRadius2 && dist>=0) {
					
					auto neighborDensity = particles[neighbor].density;
					auto neighborDensityRecip = vars.kParticleMass / neighborDensity;
					auto sqrDiff = vars.kSmoothRadius2 - dist2; // h^2 - r^2
					auto diff = vars.kSmoothRadius - dist; // h - r

					Vector3 pressureForce = Vector3();
					Vector3 pressureForceN = Vector3();
					//ideal gas law
					float particlePressure = vars.kIdealGasCoeff*(particleDensity - vars.kStandardDensity);
					float neighborPressure = vars.kIdealGasCoeff*(neighborDensity - vars.kStandardDensity);

					//tait equation
					//particlePressure = vars.kBConstant*(pow(particleDensity / vars.kStandardDensity, vars.kGamma) - 1);
					//neighborPressure = vars.kBConstant*(pow(neighborDensity / vars.kStandardDensity, vars.kGamma) - 1);

					if (dist == 0.f) {
						pressureForce = Vector3(0, 0, 0);
						pressureForceN = Vector3(0, 0, 0);
					}
					else {
						pressureForce = 0.5*(particlePressure + neighborPressure)*(particleDensityRecip)*vars.kSpiky*(r / dist)*(diff*diff);
						pressureForceN = 0.5*(particlePressure + neighborPressure)*(neighborDensityRecip)*vars.kSpiky*(r / dist)*(diff*diff);
					}
					/*auto accPressure = dist == 0.f ? Vector3(0, 0, 0) :
						vars.kIdealGasCoeff * vars.kParticleMass *
						(particleDensity + neighborDensity -
							2 * vars.kStandardDensity) * vars.kSpiky * diff * diff /
						2.f * r / dist / particleDensity / neighborDensity;
						*/
					particles[i].acceleration -= pressureForce;
					particles[neighbor].acceleration += pressureForceN;

					/*auto accViscosity =
						vars.kViscosity * (particles[neighbor].vel - particles[i].vel)
						* vars.kLaplacianViscosity * diff * vars.kParticleMass;*/
					Vector3 viscosityForce = vars.kViscosity*(particles[neighbor].vel - particles[i].vel)*vars.kLaplacianViscosity*diff*vars.kParticleMass;

					particles[i].acceleration += viscosityForce / neighborDensity;
					particles[neighbor].acceleration -= viscosityForce / particleDensity;

					//surface tension
					
					Vector3 surfaceTension = Vector3(0, 0, 0);
					Vector3 surfaceTensionN = Vector3(0, 0, 0);

					Vector3 surfNormal = neighborDensityRecip*vars.kDivPoly6*r*sqrDiff*sqrDiff;
					Vector3 surfNormalN = -particleDensityRecip*vars.kDivPoly6*r*sqrDiff*sqrDiff;

					float cs = neighborDensityRecip*vars.kLapPoly6*sqrDiff*(3 * vars.kSmoothRadius2 - 7 * dist2);
					float csN = particleDensityRecip*vars.kLapPoly6*sqrDiff*(3 * vars.kSmoothRadius2 - 7 * dist2);

					if (surfNormal.length() >= vars.kSurfaceTensionThreshold) {
						float kappa = -cs / surfNormal.length();
						float kappaN = -csN / surfNormalN.length();

						surfaceTension = kappa*vars.kParticleMass*vars.kPoly6*sqrDiff*sqrDiff*sqrDiff*r;
						surfaceTensionN = kappaN*vars.kParticleMass*vars.kPoly6*sqrDiff*sqrDiff*sqrDiff*r;



						particles[i].acceleration += kappa*vars.kSurfaceTensionSigma*surfNormal;
						particles[neighbor].acceleration -= kappaN*vars.kSurfaceTensionSigma*surfNormal;
					}
					

				
				}
			});
			
			particles[i].acceleration += particles[i].nextForce;
			particles[i].nextForce = Vector3(0, 0, 0);
		}

		externalForce *= 0.3;

	}

	const float eps = 0.5;

	void updateXSPH() {
		//#pragma omp parallel for
		for (int i = 0; i < vars.kNumParticles; i++) {
			auto particleDensity = particles[i].density;


			auto particlePosition = particles[i].pos;
			auto particleDensityRecip = vars.kParticleMass / particleDensity;

			Vector3 vSum = Vector3(0,0,0);
			forEachNeighbor(i, [&](int neighbor) {


				auto r = particlePosition - particles[neighbor].pos;
				auto dist2 = r.dot(r);// r.lengthSq();// r.dot(r); // ||r||^2
				auto dist = r.length(); // ||r||
				if (dist2 <= vars.kSmoothRadius2 && dist >= 0) {
					auto neighborDensity = particles[neighbor].density;
					auto sqrDiff = vars.kSmoothRadius2 - dist2; // h^2 - r^2
					Vector3 vSum = (eps*vars.kParticleMass / (0.5*(particleDensity + neighborDensity)))*(vars.kPoly6*sqrDiff*sqrDiff*sqrDiff*(particles[neighbor].vel-particles[i].vel));

				}
			});

			forEachNeighbor(i, [&](int neighbor) {


				auto r = particlePosition - particles[neighbor].pos;
				auto dist2 = r.dot(r);// r.lengthSq();// r.dot(r); // ||r||^2
				auto dist = r.length(); // ||r||
				if (dist2 <= vars.kSmoothRadius2 && dist >= 0) {
					auto neighborDensity = particles[neighbor].density;
					auto sqrDiff = vars.kSmoothRadius2 - dist2; // h^2 - r^2
					
					
					particles[neighbor].vel -= vSum;
				}
			});
			particles[i].vel += vSum;
		}
	}

	const float boundV = 30;

	void updateParticles() {
		//float delta = clk.deltaTime*0.1;
		//#pragma omp parallel for
		for (int p = 0; p < vars.kNumParticles; p++) {
			Vector3 &velocity = particles[p].vel;
			Vector3 &position = particles[p].pos;

			Vector3 &prePos = particles[p].pPos;
			Vector3 &preVel = particles[p].pVel;
			Vector3 &preAcc = particles[p].pAcc;
			

			if (particles[p].cX) {
				//velocity.x *= 1. / vars.kRebounceCoeff;
				particles[p].cX = false;
			}
			if (particles[p].cY) {
				//velocity.y *= 1. / vars.kRebounceCoeff;
				particles[p].cY = false;
			}
			if (particles[p].cZ) { 
				//velocity.z *= 1. / vars.kRebounceCoeff; 
				particles[p].cZ = false;
			}

			if (firstSim) {

				velocity += vars.kTimeStep * particles[p].acceleration;
				position += vars.kTimeStep * velocity;
			}
			else {

				//Leap-Frog integration - quite unstable here :(
				position += vars.kTimeStep * preVel + 0.5*preAcc*vars.kTimeStep*vars.kTimeStep;
				velocity += 0.5*vars.kTimeStep*(preAcc + particles[p].acceleration);
				
			}
			//velocity.x = min();

			
			prePos = position;
			preVel = velocity;
			preAcc = particles[p].acceleration;

			//default: 0
			if (position.x < 0) {
				position.x = 0;
				velocity.x *= -vars.kRebounceCoeff;
				particles[p].cX = true;
			}
			else if (position.x > vars.kXMax) {
				position.x = vars.kXMax;
				velocity.x *= -vars.kRebounceCoeff;
				particles[p].cX = true;
			}

			//default: 0
			if (position.y < particleSize) {
				position.y = particleSize;
				velocity.y *= -vars.kRebounceCoeff;
				particles[p].cY = true;
			}
			else if (position.y > vars.kYMax) {
				position.y = vars.kYMax;
				velocity.y *= -vars.kRebounceCoeff;
				particles[p].cY = true;
			}

			//default: 0
			if (position.z < 0) {
				position.z = 0;
				velocity.z *= -vars.kRebounceCoeff;
				particles[p].cZ = true;
			}
			else if (position.z > vars.kZMax) {
				position.z = vars.kZMax;
				velocity.z *= -vars.kRebounceCoeff;
				particles[p].cZ = true;
			}

			//clamp velocity for stability
			if (velocity.x < -boundV)velocity.x = -boundV;
			if (velocity.y < -boundV)velocity.y = -boundV;
			if (velocity.z < -boundV)velocity.z = -boundV;

			if (velocity.x > boundV)velocity.x = boundV;
			if (velocity.y > boundV)velocity.y = boundV;
			if (velocity.z > boundV)velocity.z = boundV;

		}

		//firstSim = false;
	}

	void updateParticlesRK4() {
		for (int p = 0; p < vars.kNumParticles; p++) {
			Vector3 &velocity = particles[p].vel;
			Vector3 &position = particles[p].pos;

			Vector3 v0 = velocity;
			float t0 = simTime;
			float t1 = t0 + vars.kTimeStep*0.5;
			Vector3 v1 = v0 + t1*particles[p].acceleration;
			float t2 = t0 + vars.kTimeStep*0.5;
			Vector3 v2 = v1 + t2*particles[p].acceleration;
			float t3 = t0 + vars.kTimeStep;
			Vector3 v3 = v2 + t3*particles[p].acceleration;

			velocity += (v0 + v1*2 + v2*2 + v3)*vars.kTimeStep;

			//velocity += vars.kTimeStep * particles[p].acceleration;
			//position += vars.kTimeStep * velocity;
			position += vars.kTimeStep*velocity;

			//default: 0
			if (position.x < 0) {
				position.x = 0;
				velocity.x *= -vars.kRebounceCoeff;
			}
			else if (position.x > vars.kXMax) {
				position.x = vars.kXMax;
				velocity.x *= -vars.kRebounceCoeff;
			}

			//default: 0
			if (position.y < particleSize) {
				position.y = particleSize;
				velocity.y *= -vars.kRebounceCoeff;
			}
			else if (position.y > vars.kYMax) {
				position.y = vars.kYMax;
				velocity.y *= -vars.kRebounceCoeff;
			}

			//default: 0
			if (position.z < 0) {
				position.z = 0;
				velocity.z *= -vars.kRebounceCoeff;
			}
			else if (position.z > vars.kZMax) {
				position.z = vars.kZMax;
				velocity.z *= -vars.kRebounceCoeff;

			}

		}
	}
	

	glm::vec3 glm(Vector3 v) {
		return glm::vec3(v.x,v.y,v.z);
	}

	void setReverseRotation(const glm::mat3 &model) {
		reverseRotation = model;

		setAxis(glm::normalize(reverseRotation * glm::vec3(0.f, 1.f, 0.f)));
	}

	void setAxis(const glm::vec3 &axis) {
		this->axis.x = axis.x;
		this->axis.y = axis.y;
		this->axis.z = axis.z;
	}

	void setExternalForce(float externalForce) {
		externalForce = std::min(externalForce, 350.f);
		this->externalForce = externalForce;
	}

	void clear() {
		//std::fill(start.begin(), start.end(), -1);
		//std::fill(end.begin(), end.end(), -1);
		//std::fill(next.begin(), next.end(), -1);

		starts.clear();
		starts.resize(vars.kGridSize, -1);

		ends.clear();
		ends.resize(vars.kGridSize, - 1);

		starts_.clear();
		starts_.resize(vars.kGridSize, - 1);

		ends_.clear();
		ends_.resize(vars.kGridSize,-1);

		start.clear();
		start.resize(vars.kGridSize,-1);

		end.clear();
		end.resize(vars.kGridSize,-1);

		next.clear();
		next.resize(vars.kNumParticles, -1);

		nexts_.clear();
		nexts_.resize(vars.kNumParticles, -1);

		velocities_.clear();
		velocities_.resize(vars.kNumParticles, Vector3(0,0,0));

		positions_.clear();
		positions_.resize(vars.kNumParticles, Vector3(0, 0, 0));
	}

	void setTriangleData(vector<Triangle> triangles) {

		for (int i = 0; i < nParticles; i++) {

			particles[i].setData(triangles);
			
		}
	}

	void setData(vector<Vector3> vertices, vector<int> indices) {
		for (int i = 0; i < nParticles; i++) {

			particles[i].setData(vertices, indices);
		}
	}

	void setData(vector<Vector3> vertices, vector<int> indices, vector<Models::Vertex> vertexes) {
		for (int i = 0; i < nParticles; i++) {

			particles[i].setData(vertices, indices, vertexes);
			//particles[i].updateMesh();
		}
	}

	void setupParticlesInstanced() {
		for (int i = 0; i < nParticles; i++) {

			particles[i].setupBuffersInstanced();
		}
	}

	void setupParticles() {
		for (int i = 0; i < nParticles; i++) {

			particles[i].setupBuffers();
		}
	}

	void terminate() {}

};
