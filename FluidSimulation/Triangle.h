#pragma once

#include <vector>
#include "Vectors.h"

using namespace std;

class Triangle {


public:

	Vector3 v1, v2, v3;

	float a, b, c, d;

	int i1, i2, i3;

public:

	Triangle() {

	}

	Triangle(Vector3 v1, Vector3 v2, Vector3 v3) {
		this->v1 = v1;
		this->v2 = v2;
		this->v3 = v3;
	}

	void move(Vector3 v) {
		v1 = v1 + v;
		v2 = v2 + v;
		v3 = v3 + v;

	}

	Vector3 linear1(float amount) {
		Vector3 res = Vector3();
		res.x = v1.x*amount + v2.x*(1 - amount);
		res.y = v1.y*amount + v2.y*(1 - amount);
		res.z = v1.z*amount + v2.z*(1 - amount);
		return res;
	}

	Vector3 linear2(float amount) {
		Vector3 res = Vector3();
		res.x = v2.x*amount + v3.x*(1 - amount);
		res.y = v2.y*amount + v3.y*(1 - amount);
		res.z = v2.z*amount + v3.z*(1 - amount);
		return res;
	}

	Vector3 linear3(float amount) {
		Vector3 res = Vector3();
		res.x = v3.x*amount + v1.x*(1 - amount);
		res.y = v3.y*amount + v1.y*(1 - amount);
		res.z = v3.z*amount + v1.z*(1 - amount);
		return res;
	}

	Vector3 centroid() {
		Vector3 res = Vector3();
		res.x = (v1.x + v2.x + v3.x) / 3;
		res.y = (v1.y + v2.y + v3.y) / 3;
		res.z = (v1.z + v2.z + v3.z) / 3;

		res = (res*normal())*(sqrt(2. / 3.)*edgeDistance());
		return res;
	}

	Vector3 closest(Vector3 v) {

		float dist1 = (v - v1).length();
		float dist2 = (v - v2).length();
		float dist3 = (v - v3).length();

		if (dist1 < dist2 && dist1 < dist3)return v1;
		if (dist2 < dist1 && dist2 < dist3)return v2;

		return v3;
	}

	float getLongEdge() {
		
		return fmax(fmax((v1-v2).length(), (v1 - v3).length()),(v2-v3).length());
	}

	Vector3 getCenter() {

		return (v1 + v2 + v3) / 3;
	}

	bool intersects(Triangle t) {

		computePlane();

		bool isV1Left = (a*t.v1.x + b*t.v1.y + c*t.v1.z + d) < 0;
		bool isV2Left = (a*t.v2.x + b*t.v2.y + c*t.v2.z + d) < 0;
		bool isV3Left = (a*t.v3.x + b*t.v3.y + c*t.v3.z + d) < 0;

		bool leftV[3] = {isV1Left,isV2Left,isV3Left};

		int countLeft = 0;
		for (int i = 0; i < 3; i++) {
			if (leftV[i] < 0)countLeft++;
		}

		if (countLeft > 0 && countLeft < 3) return true;

		return false;
	}

	void computePlane() {
		this->a = v1.y * (v2.z - v3.z) + v2.y * (v3.z - v1.z) + v3.y * (v1.z - v2.z);
		this->b = v1.z * (v2.x - v3.x) + v2.z * (v3.x - v1.x) + v3.z * (v1.x - v2.x);
		this->c = v1.x * (v2.y - v3.y) + v2.x * (v3.y - v1.y) + v3.x * (v1.y - v2.y);
		this->d = -(v1.x * (v2.y * v3.z - v3.y * v2.z) + v2.x * (v3.y * v1.z - v1.y * v3.z) + v3.x * (v1.y * v2.z - v2.y * v1.z));
	}

	float edgeDistance() {
		return sqrt((this->v1.x - this->v2.x)*(this->v1.x - this->v2.x) + (this->v1.y - this->v2.y)*(this->v1.y - this->v2.y) + (this->v1.z - this->v2.z)*(this->v1.z - this->v2.z));
	}

	Vector3 normal() {
		Vector3 res = Vector3();

		
		
		//res = (v2 - v1).cross(v3 - v1);
		
		//res = (v3 - v2).cross(v1-v2);
		res = (v2 - v3).cross(v1 - v2);
		
		res = res.normalize();
		return res;
	}

	Vector3 normalU() {
		Vector3 res = Vector3();

		res = (v2 - v3).cross(v1 - v2);

		return res;
	}
};