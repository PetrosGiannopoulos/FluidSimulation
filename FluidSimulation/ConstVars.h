#pragma once

#include <cmath>
#include "Vectors.h"
#include <glm/glm.hpp>
class ConstVars {

public:

	double kPi = 3.14159265358979323836;
	int kNumParticles = 5000;
	float kSmoothRadius = 1.5;
	float kSmoothRadius2 = kSmoothRadius * kSmoothRadius;
	float kSmoothRadius3 = kSmoothRadius2 * kSmoothRadius;
	float kPoly6 = 315.0 / 64.0 / kPi / kSmoothRadius3 / kSmoothRadius3 / kSmoothRadius3;
	float kSpiky = -45.0 / kPi / kSmoothRadius3 / kSmoothRadius3;
	//0 - 1.5 
	float kViscosity = 0.25;// 0.25;
	float kParticleMass = 1.0;// 1.0;
	float kIdealGasCoeff = 1000.0;
	float kStandardDensity = 1.2;
	Vector3 kGravity = Vector3(0.0, -100, 0.0);
	float kLaplacianViscosity = 45.0 * kPi / kSmoothRadius3 / kSmoothRadius3;
	float kDivPoly6 = -945.0 / 32.0 / kPi / kSmoothRadius3 / kSmoothRadius3 / kSmoothRadius3;
	float kLapPoly6 = kDivPoly6;
	float kSurfaceTensionThreshold = 0.025;
	float kSurfaceTensionSigma = 0.0728;
	float kGamma = 7;
	float kBConstant = 1119;

	//default: 30
	float defaultSize = 15;
	float kXMax = 30;
	float kYMax = 30;
	float kZMax = 30;
	int kNX;
	int kNY;
	int kNZ;
	int kGridSize;

	//default: 0.02
	float kTimeStep = 0.02f;
	//default: 0.9
	float kRebounceCoeff = 0.9f;

	//default: 1.0f
	float kMcStep = 0.75;
	float kMcMargin = kSmoothRadius / kMcStep;
	int kMcNx;
	int kMcNy;
	int kMcNz;
	int kMcGridSize;

	float kPolyThreshold = kParticleMass * .6f / kPoly6;

public:

	ConstVars(){
		kNX = (kXMax / kSmoothRadius + 1);
		kNY = (kYMax / kSmoothRadius + 1);
		kNZ = (kZMax / kSmoothRadius + 1);

		kGridSize = kNX*kNY*kNZ;

		kMcNx = ((kXMax + 2 * kMcMargin) / kMcStep + 1);
		kMcNy = ((kYMax + 2 * kMcMargin) / kMcStep + 1);
		kMcNz = ((kZMax + 2 * kMcMargin) / kMcStep + 1);

		kMcGridSize = kMcNx*kMcNy*kMcNz;
	}

	void recalculate() {
		kNX = (kXMax / kSmoothRadius + 1);
		kNY = (kYMax / kSmoothRadius + 1);
		kNZ = (kZMax / kSmoothRadius + 1);

		kGridSize = kNX*kNY*kNZ;

		kMcNx = ((kXMax + 2 * kMcMargin) / kMcStep + 1);
		kMcNy = ((kYMax + 2 * kMcMargin) / kMcStep + 1);
		kMcNz = ((kZMax + 2 * kMcMargin) / kMcStep + 1);

		kMcGridSize = kMcNx*kMcNy*kMcNz;
	}
};