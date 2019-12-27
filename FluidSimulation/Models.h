#pragma once

#define _USE_MATH_DEFINES
#include <vector>
#include "Vectors.h"
#include "Triangle.h"

#include <math.h>
#include <fstream>
#include <sstream>
#include <string>
#include <algorithm>

using namespace std;


class Models {

public:

	struct Vertex {
		// position
		Vector3 Position;
		// normal
		Vector3 Normal;
		// texCoords
		Vector2 TexCoords;

		Vector3 Tangent;
		Vector3 BiTangent;

		float wTangent;
		int faceCount;
		vector<Vector3> wnormals;
		Vector3 normal;
	};

	struct TriangleSSBO {

		glm::vec4 v1;
		glm::vec4 v2;
		glm::vec4 v3;
		glm::vec4 n1;
		glm::vec4 n2;
		glm::vec4 n3;
	};

	vector<TriangleSSBO> ssbo_triangles;
	vector<Vector3> m_vertices;
	vector<Vertex> m_vertexes;
	vector<Triangle> m_triangles;
	vector<int> m_indices;
	vector<Vector3> m_normals;
	vector<Vector3> m_gradients;
	vector<int> m_faceCount;
	vector<vector<int>> m_commonTris;

	glm::vec3 AABB_min;
	glm::vec3 AABB_max;
	//vector<Vector3> 

public :

	Models() {

	}

	void loadModel(char *s) {
		//read file
		std::ifstream infile(s);
		std::string line;

		//for splitting string into array
		string arr[4];
		int i = 0;

		//parse every line
		while (std::getline(infile, line)) {
			stringstream ssin(line);

			if (line.size() == 0)continue;
			if (line[0] == 'v') {

				while (ssin.good() && i < 4) {
					ssin >> arr[i]; ++i;
				}
				m_vertices.push_back(Vector3(stof(arr[1]), stof(arr[2]), stof(arr[3])));

				Vertex vertex = Vertex();
				vertex.Position = m_vertices[m_vertices.size() - 1];
				m_vertexes.push_back(vertex);
			}
			else if (line[0] == 'f') {
				while (ssin.good() && i < 4) {
					ssin >> arr[i]; ++i;
				}
				m_triangles.push_back(Triangle(m_vertices[stoi(arr[1]) - 1], m_vertices[stoi(arr[2]) - 1], m_vertices[stoi(arr[3]) - 1]));
				//showTriangles.push_back(true);
				m_indices.push_back(stoi(arr[1]) - 1);
				m_indices.push_back(stoi(arr[2]) - 1);
				m_indices.push_back(stoi(arr[3]) - 1);

			}
			i = 0;
		}

	}

	void clearAll() {
		m_gradients.clear();
		m_indices.clear();
		m_vertexes.clear();
		m_vertices.clear();
		m_normals.clear();
		m_triangles.clear();
		m_faceCount.clear();
		m_commonTris.clear();
		ssbo_triangles.clear();

	}

	void calcAABB() {

		int N = m_vertices.size();

		float minX = FLT_MAX;
		float minY = FLT_MAX;
		float minZ = FLT_MAX;

		float maxX = -FLT_MAX;
		float maxY = -FLT_MAX;
		float maxZ = -FLT_MAX;

		for (int i = 0; i < N;i++) {

			glm::vec3 v = glm::vec3(m_vertices[i].x, m_vertices[i].y, m_vertices[i].z);

			if (v.x < minX)minX = v.x;
			if (v.y < minY)minY = v.y;
			if (v.z < minZ)minZ = v.z;

			if (v.x > maxX)maxX = v.x;
			if (v.y > maxY)maxY = v.y;
			if (v.z > maxZ)maxZ = v.z;

		}

		AABB_min = glm::vec3(minX,minY,minZ);
		AABB_max = glm::vec3(maxX,maxY,maxZ);

	}

	void rescaleModel(float size) {


		Vector3 cm = Vector3(0,0,0);

		int N = m_vertices.size();

		for (int i = 0; i<N; i++) {
			Vector3 vertex = m_vertices[i];

			cm += vertex;

		}
		cm /= N;

		//cout << N << endl;
		
		for (int i = 0; i < N; i++) {
			Vector3 &vertex = m_vertices[i];

			vertex *= size;

			m_vertexes[i].Position = vertex;
		}
		
		updateTriangles();
	}

	void refineModel(int k_) {

		int N = m_triangles.size();

		for (int k = 0; k<k_; k++) {
			N = m_triangles.size();
			for (int i = 0; i<N; i++) {
				Triangle tri = m_triangles[i];

				Vector3 n1 = tri.linear1(0.5);
				Vector3 n2 = tri.linear2(0.5);
				Vector3 n3 = tri.linear3(0.5);

				//Vector3D norm = tri->normal();
				//n1 = (*n1) + (*((*norm)*(radius-n1->length())));
				//n2 = (*n2) + (*((*norm)*(radius-n2->length())));
				//n3 = (*n3) + (*((*norm)*(radius-n3->length())));
				n1 = n1.normalize();
				n2 = n2.normalize();
				n3 = n3.normalize();
				tri.v1 = tri.v1.normalize();
				tri.v2 = tri.v2.normalize();
				tri.v3 = tri.v3.normalize();

				//Vector3D *center = tri->centroid();
				Triangle t1 = Triangle(tri.v1, n1, n3);
				Triangle t2 = Triangle(tri.v2, n2, n1);
				Triangle t3 = Triangle(tri.v3, n3, n2);
				Triangle t4 = Triangle(n1, n2, n3);


				//m_triangles.push_back(t1);
				m_triangles[i] = t1;
				m_triangles.push_back(t2); m_triangles.push_back(t3); m_triangles.push_back(t4);
				//showTriangles.push_back(true); showTriangles.push_back(true); showTriangles.push_back(true); showTriangles.push_back(true);
				//showTriangles[i]=false;

			}
		}

		
	}

	void refineMCModel(int k_) {

		int N = m_triangles.size();

		for (int k = 0; k<k_; k++) {
			N = m_triangles.size();
			for (int i = 0; i<N; i++) {
				Triangle tri = m_triangles[i];

				Vector3 n1 = tri.linear1(0.5);
				Vector3 n2 = tri.linear2(0.5);
				Vector3 n3 = tri.linear3(0.5);

				n1 = n1.normalize();
				n2 = n2.normalize();
				n3 = n3.normalize();
				tri.v1 = tri.v1.normalize();
				tri.v2 = tri.v2.normalize();
				tri.v3 = tri.v3.normalize();

				int last = m_indices.size();

				int i11 = last, i12 = last + 1, i13 = last + 2;
				m_indices.push_back(i11); m_indices.push_back(i12); m_indices.push_back(i13);

				last = m_indices.size();

				int i21 = last, i22 = last + 1, i23 = last + 2;
				m_indices.push_back(i21); m_indices.push_back(i22); m_indices.push_back(i23);

				last = m_indices.size();

				int i31 = last, i32 = last + 1, i33 = last + 2;
				m_indices.push_back(i31); m_indices.push_back(i32); m_indices.push_back(i33);

				last = m_indices.size();

				int i41 = last, i42 = last + 1, i43 = last + 2;
				m_indices.push_back(i41); m_indices.push_back(i42); m_indices.push_back(i43);

				m_vertices.push_back(tri.v1); m_vertices.push_back(n1); m_vertices.push_back(n3);
				m_vertices.push_back(tri.v2); m_vertices.push_back(n2); m_vertices.push_back(n1);
				m_vertices.push_back(tri.v3); m_vertices.push_back(n3); m_vertices.push_back(n2);
				m_vertices.push_back(n1); m_vertices.push_back(n2); m_vertices.push_back(n3);
				
				Triangle t1 = Triangle(tri.v1, n1, n3);
				Triangle t2 = Triangle(tri.v2, n2, n1);
				Triangle t3 = Triangle(tri.v3, n3, n2);
				Triangle t4 = Triangle(n1, n2, n3);

				m_normals.push_back(-t1.normalU()); m_normals.push_back(-t1.normalU()); m_normals.push_back(-t1.normalU());
				m_normals.push_back(-t2.normalU()); m_normals.push_back(-t2.normalU()); m_normals.push_back(-t2.normalU());
				m_normals.push_back(-t3.normalU()); m_normals.push_back(-t3.normalU()); m_normals.push_back(-t3.normalU());
				m_normals.push_back(-t4.normalU()); m_normals.push_back(-t4.normalU()); m_normals.push_back(-t4.normalU());

				m_gradients.push_back(-t1.normal()); m_gradients.push_back(-t1.normal()); m_gradients.push_back(-t1.normal());
				m_gradients.push_back(-t2.normal()); m_gradients.push_back(-t2.normal()); m_gradients.push_back(-t2.normal());
				m_gradients.push_back(-t3.normal()); m_gradients.push_back(-t3.normal()); m_gradients.push_back(-t3.normal());
				m_gradients.push_back(-t4.normal()); m_gradients.push_back(-t4.normal()); m_gradients.push_back(-t4.normal());

				Vertex vertex11 = Vertex();vertex11.Position = tri.v1;vertex11.Normal = -t1.normal();vertex11.faceCount = 0;
				Vertex vertex12 = Vertex(); vertex12.Position = n1; vertex12.Normal = -t1.normal(); vertex12.faceCount = 0;
				Vertex vertex13 = Vertex(); vertex13.Position = n3; vertex13.Normal = -t1.normal(); vertex13.faceCount = 0;
				
				Vertex vertex21 = Vertex(); vertex21.Position = tri.v2; vertex21.Normal = -t2.normal(); vertex21.faceCount = 0;
				Vertex vertex22 = Vertex(); vertex22.Position = n2; vertex22.Normal = -t2.normal(); vertex22.faceCount = 0;
				Vertex vertex23 = Vertex(); vertex23.Position = n1; vertex23.Normal = -t2.normal(); vertex23.faceCount = 0;
				
				Vertex vertex31 = Vertex(); vertex31.Position = tri.v3; vertex31.Normal = -t3.normal(); vertex31.faceCount = 0;
				Vertex vertex32 = Vertex(); vertex32.Position = n3; vertex32.Normal = -t3.normal(); vertex32.faceCount = 0;
				Vertex vertex33 = Vertex(); vertex33.Position = n2; vertex33.Normal = -t3.normal(); vertex33.faceCount = 0;
				
				Vertex vertex41 = Vertex(); vertex41.Position = n1; vertex41.Normal = -t4.normal(); vertex41.faceCount = 0;
				Vertex vertex42 = Vertex(); vertex42.Position = n2; vertex42.Normal = -t4.normal(); vertex42.faceCount = 0;
				Vertex vertex43 = Vertex(); vertex43.Position = n3; vertex43.Normal = -t4.normal(); vertex43.faceCount = 0;
				
				m_triangles[i] = t1;
				m_triangles.push_back(t2); m_triangles.push_back(t3); m_triangles.push_back(t4);
				
			}
		}
	}

	Vector3 cm() {
		
		int N = m_vertices.size();
		Vector3 center = Vector3(0,0,0);
		Vector3 move = Vector3(0,-2.5, -30);
		for (int i = 0; i < N; i++) {
			Vector3 v = m_vertices[i];

			center += v;
		}
		center /= (float)N;
		return center;
	}

	float PI = 3.14159265359;

	void flattenMesh() {


		int N = m_gradients.size();
		
		for (int i = 0; i < N; i++) {

			Vector3 g0 = m_gradients[i];
		
			float r = g0.length();
			float phi = atan2(g0.y,g0.x)*(180/PI)+180;
			float theta = acos(g0.z/r)*(180/PI);

			//cout << phi << endl;

			float r2 = r;
			float phi2=0;
			float theta2=0;


			if (phi > 0 && phi < 90) {
				if (phi < abs(phi - 90))phi2 = 0;
				else phi2 = 90;
			}

			if (phi > 90 && phi < 180) {
				if (abs(phi-90) < abs(phi - 180))phi2 = 90;
				else phi2 = 180;
			}

			if (phi > 180 && phi < 270) {
				if (abs(phi - 180) < abs(phi - 270))phi2 = 180;
				else phi2 = 270;
			}

			if (phi > 270 && phi < 360) {
				if (abs(phi - 270) < abs(phi - 360))phi2 = 270;
				else phi2 = 360;
			}

			if (theta > 0 && theta < 90) {
				if (theta < abs(theta - 90))theta2 = 0;
				else theta2 = 90;
			}

			if (theta > 90 && theta < 180) {
				if (abs(theta - 90) < abs(theta - 180))theta2 = 90;
				else theta2 = 180;
			}

			phi2 -= 180;

			float x = r2*sin(theta2*(PI/180))*cos(phi2*(PI / 180));
			float y = r2*sin(theta2*(PI / 180))*sin(phi2*(PI / 180));
			float z = r2*cos(theta2*(PI / 180));
			
			m_gradients[i] = Vector3(x,y,z);

		}

		

	}

	void taubinSmoothing() {

	
		float lambda = 0.50;
		float mi = 0.53;

		int N = m_triangles.size();
		int vN = m_vertices.size();
		vector<int> p;
		m_commonTris.clear();
		m_commonTris.resize(vN, p);
		for (int i = 0; i < N; i++) {

			Triangle t = m_triangles[i];

			int i1 = t.i1;
			int i2 = t.i2;
			int i3 = t.i3;

			m_commonTris[i1].push_back(i);
			m_commonTris[i2].push_back(i);
			m_commonTris[i3].push_back(i);

		}

		vector<Vector3> n_gradients;

		for (int i = 0; i < vN;i++) {
			n_gradients.push_back(m_gradients[i]);
		}

		for (int k = 0; k < 10;k++) {

			for (int i = 0; i < vN;i++) {

				int cN = m_commonTris[i].size();

				Vector3 g0 = n_gradients[i];

				Vector3 sum = Vector3(0, 0, 0);
				float sumWeight = 0;
				for (int j = 0; j < cN;j++) {

					Triangle t = m_triangles[m_commonTris[i][j]];

					Vector3 vs[3] = {t.v1,t.v2,t.v3};
					int ids[3] = { t.i1 ,t.i2,t.i3 };
					
					for (int l = 0; l < 3;l++) {

						Vector3 vl = vs[l];
						int idl = ids[l];

						if(m_indices[i] !=idl){
						//if (!vl.equal(m_vertices[i], 0.0000005)) {
							
							float dist = vl.distance(m_vertices[i]);

							float weight = 1 / dist;
							sumWeight += weight;
							sum += weight*vl;
						}

					}

				}
				g0 = g0 + lambda * (sum / (sumWeight)-g0);
				n_gradients[i] = g0;
			}


			for (int i = 0; i < vN; i++) {

				int cN = m_commonTris[i].size();

				Vector3 g0 = n_gradients[i];

				Vector3 sum = Vector3(0, 0, 0);
				float sumWeight = 0;
				for (int j = 0; j < cN; j++) {

					Triangle t = m_triangles[m_commonTris[i][j]];

					Vector3 vs[3] = { t.v1,t.v2,t.v3 };
					int ids[3] = { t.i1 ,t.i2,t.i3 };

					for (int l = 0; l < 3; l++) {

						Vector3 vl = vs[l];
						int idl = ids[l];

						if(m_indices[i]!=idl){
						//if (!vl.equal(m_vertices[i], 0.0000005)) {

							float dist = vl.distance(m_vertices[i]);

							float weight = 1 / dist;
							sumWeight += weight;
							sum += weight*vl;
						}

					}

				}
				g0 = g0 - mi * (sum / (sumWeight)-g0);
				n_gradients[i] = g0;
			}

		}

		for (int i = 0; i < vN;i++) {
			m_gradients[i] = n_gradients[i];
		}

	}

	void buildMesh() {

		m_faceCount.clear();
		m_triangles.clear();
		m_vertexes.clear();
		m_normals.clear();
		m_normals.resize(m_vertices.size(),Vector3(0,0,0));

		m_faceCount.resize(m_vertices.size(),0);

		int N = m_indices.size();

		for (int i = 0; i < N-3; i+=3) {

			int i1 = m_indices[i];
			int i2 = m_indices[i+1];
			int i3 = m_indices[i+2];

			//showme(Vector3(i1, i2, i3), "ind");

			Vector3 v1 = m_vertices[i1];
			Vector3 v2 = m_vertices[i2];
			Vector3 v3 = m_vertices[i3];

			Triangle t = Triangle(v1,v2,v3);
			t.i1 = i1;
			t.i2 = i2;
			t.i3 = i3;
			m_triangles.push_back(t);
			Vector3 normal = -t.normalU();

			m_normals[i1] += m_gradients[i1];
			m_normals[i2] += m_gradients[i2];
			m_normals[i3] += m_gradients[i3];

			m_faceCount[i1] ++;
			m_faceCount[i2] ++;
			m_faceCount[i3] ++;
		}

		

	}

	void fillSSBO() {

		int N = m_triangles.size();
		for (int i = 0; i < N; i++) {
			
			Triangle t = m_triangles[i];

			Vector3 v1 = t.v1;
			Vector3 v2 = t.v2;
			Vector3 v3 = t.v3;

			Vector3 n1 = -t.normalU();
			Vector3 n2 = -t.normalU();
			Vector3 n3 = -t.normalU();

			TriangleSSBO ssbo_tri = TriangleSSBO();
			ssbo_tri.v1 = glm::vec4(v1.x, v1.y, v1.z,0);
			ssbo_tri.v2 = glm::vec4(v2.x, v2.y, v2.z, 0);
			ssbo_tri.v3 = glm::vec4(v3.x, v3.y, v3.z, 0);
			ssbo_tri.n1 = glm::vec4(n1.x, n1.y, n1.z, 0);
			ssbo_tri.n2 = glm::vec4(n2.x, n2.y, n2.z, 0);
			ssbo_tri.n3 = glm::vec4(n3.x, n3.y, n3.z, 0);

			ssbo_triangles.push_back(ssbo_tri);
		}
	}

	void fillSSBO_Gradients() {

		int N = m_triangles.size();
		for (int i = 0; i < N; i++) {

			Triangle t = m_triangles[i];

			Vector3 v1 = t.v1;
			Vector3 v2 = t.v2;
			Vector3 v3 = t.v3;

			Vector3 n1 = m_gradients[t.i1];
			Vector3 n2 = m_gradients[t.i2];
			Vector3 n3 = m_gradients[t.i3];

			TriangleSSBO ssbo_tri = TriangleSSBO();
			ssbo_tri.v1 = glm::vec4(v1.x, v1.y, v1.z, 0);
			ssbo_tri.v2 = glm::vec4(v2.x, v2.y, v2.z, 0);
			ssbo_tri.v3 = glm::vec4(v3.x, v3.y, v3.z, 0);
			ssbo_tri.n1 = glm::vec4(n1.x, n1.y, n1.z, 0);
			ssbo_tri.n2 = glm::vec4(n2.x, n2.y, n2.z, 0);
			ssbo_tri.n3 = glm::vec4(n3.x, n3.y, n3.z, 0);

			ssbo_triangles.push_back(ssbo_tri);
		}

	}

	void buildMeshTaubin() {
		
		int N = m_indices.size();

		int triN = 0;
		for (int i = 0; i < N - 2; i += 3) {

			int i1 = m_indices[i];
			int i2 = m_indices[i + 1];
			int i3 = m_indices[i + 2];

			//showme(Vector3(i1, i2, i3), "ind");

			Vector3 v1 = m_vertices[i1];
			Vector3 v2 = m_vertices[i2];
			Vector3 v3 = m_vertices[i3];

			Triangle t = Triangle(v1, v2, v3);
			t.i1 = i1;
			t.i2 = i2;
			t.i3 = i3;
			m_triangles[triN] = t;
			Vector3 normal = -t.normalU();

			m_normals[i1] += normal;// m_gradients[i1];
			m_normals[i2] += normal;// m_gradients[i2];
			m_normals[i3] += normal;// m_gradients[i3];

			triN++;
		}
	}

	void buildVertexes() {
		int N = m_indices.size();

		for (int i = 0; i < N; i++) {
			int i0 = m_indices[i];

			Vector3 v0 = m_vertices[i0];
			Vector3 n0 = m_normals[i0];
			Vector3 g0 = m_gradients[i0];

			Vertex vertex0 = Vertex();
			vertex0.Position = v0;
			vertex0.Normal = g0.normalize();
			vertex0.faceCount = 0;
			m_vertexes.push_back(vertex0);
		}
	}

	void vertexNormals() {

		int triN = m_triangles.size();
		int vertN = m_vertices.size();


		for (int i = 0; i < vertN; i++) {

			m_vertexes[i].normal = Vector3(0,0,0);
			m_vertexes[i].faceCount = 0;
		}
		for (int i = 0; i < triN; i++) {

			Triangle t = m_triangles[i];

			int i1 = t.i1;
			int i2 = t.i2;
			int i3 = t.i3;

			m_vertexes[i1].normal += m_gradients[i1];
			m_vertexes[i2].normal += m_gradients[i2];
			m_vertexes[i3].normal += m_gradients[i3];

			m_vertexes[i1].faceCount++;
			m_vertexes[i2].faceCount++;
			m_vertexes[i3].faceCount++;
		}

		
		for (int i = 0; i < vertN;i++) {

			m_vertexes[i].normal /= (float)m_vertexes[i].faceCount;

			m_vertexes[i].Normal = m_vertexes[i].normal;
			//m_vertexes[i].Normal = m_vertexes[i].normal.normalize();
		}

	}

	void smoothNormals() {

		int triN = m_triangles.size();

		for (int i = 0; i < triN; i++) {

			Triangle t = m_triangles[i];

			Vector3 p1 = t.v1;
			Vector3 p2 = t.v2;
			Vector3 p3 = t.v3;

			int i1 = t.i1;
			int i2 = t.i2;
			int i3 = t.i3;

			Vector3 n = -t.normalU(); //(p2 - p1).cross(p3 - p1);

			float a1 = angle(p2-p1, p3-p1);
			float a2 = angle(p3-p2, p1-p2);
			float a3 = angle(p1-p3, p2-p3);

			m_vertexes[i1].wnormals.push_back(n*a1);
			m_vertexes[i2].wnormals.push_back(n*a2);
			m_vertexes[i3].wnormals.push_back(n*a3);
		}

		int vertN = m_vertices.size();

		for (int i = 0; i < vertN;i++) {

			Vector3 n = Vector3(0,0,0);

			int wSize = m_vertexes[i].wnormals.size();

			for (int j = 0; j < wSize;j++) {
				n += m_vertexes[i].wnormals[j];
			}
			
			n = n.normalize();
			m_normals[i] = n;
			m_vertexes[i].Normal = n;
		}
	}

	float angle(Vector3 v1, Vector3 v2) {
		
		float l1 = v1.lengthSq();
		float l2 = v2.lengthSq();

		float dot = v1.dot(v2);

		return acos(dot/sqrt(l1*l2));

		//return acos(v1.normalize().dot(v2.normalize()));

	}

	void showme(Vector3 v, string s = "") {

		cout << s << "|" << v.x << "|" << v.y << "|" << v.z << endl;
	}

	void rebuildStructure() {

		m_vertices.clear();
		m_indices.clear();
		m_vertexes.clear();
		m_normals.clear();
		
		int N = m_triangles.size();

		for (int i = 0; i < N;i++) {

			Triangle tri = m_triangles[i];
			
			Vector3 v0 = tri.v1;
			Vector3 v1 = tri.v2;
			Vector3 v2 = tri.v3;

			Vector3 normal = -tri.normalU();// -tri.normalU();

			Vertex vertex0 = Vertex();
			vertex0.Position = v0;
			vertex0.Normal = normal;
			Vertex vertex1 = Vertex();
			vertex1.Position = v1;
			vertex1.Normal = normal;
			Vertex vertex2 = Vertex();
			vertex2.Position = v2;
			vertex2.Normal = normal;

			m_vertices.push_back(v0);
			m_vertices.push_back(v1);
			m_vertices.push_back(v2);

			int vN = m_vertices.size()-1;

			m_indices.push_back(vN - 2);
			m_indices.push_back(vN - 1);
			m_indices.push_back(vN);

			m_normals.push_back(normal);
			m_normals.push_back(normal);
			m_normals.push_back(normal);

			m_vertexes.push_back(vertex0);
			m_vertexes.push_back(vertex1);
			m_vertexes.push_back(vertex2);

		}

	}

	void redefineStructure() {

		float epsilon = 1.401298E-45;
		int N = m_triangles.size();

		m_vertices.clear();
		m_indices.clear();
		m_vertexes.clear();

		int V;

		Vector3 sphereCenter = cm();
		for (int i = 0; i < N;i++) {

			Triangle t = m_triangles[i];

			Vector3 v1 = t.v1;
			Vector3 v2 = t.v2;
			Vector3 v3 = t.v3;

			Vector3 n = t.normal();


			V = m_vertices.size();
			bool fV1 = false; bool fV2 = false; bool fV3 = false;
			for (int j = 0; j < V;j++) {

				//if (j % 3 == 1 || j %3 == 2)continue;

				if (v1.equal(m_vertices[j], epsilon)) {
					fV1 = true;
					break;
				}
				if (v2.equal(m_vertices[j], epsilon)) {
					fV2 = true;
					break;
				}
				if (v3.equal(m_vertices[j], epsilon)) {
					fV3 = true;
					break;
				}
				
			}
			
			if (fV1 == false) {

				

				m_vertices.push_back(v1);
				//m_vertices.push_back(n);

				//sphere center = (0,0,0)
				Vector3 sc = (v1-sphereCenter).normalize();
				float u = atan2(sc.z, sc.x) / (2 * M_PI) + 0.5;
				float v = 0.5 - asin(sc.y)/M_PI;
				Vector2 uv = Vector2(u, v);
				//m_vertices.push_back(uv);

				Vertex vertex = Vertex();
				vertex.Position = v1;
				vertex.Normal = n;
				vertex.TexCoords = uv;
				m_vertexes.push_back(vertex);
			}
			if (fV2 == false) {
				m_vertices.push_back(v2);
				//m_vertices.push_back(n);

				//sphere center = (0,0,0)
				Vector3 sc = (v2-sphereCenter).normalize();
				float u = atan2(sc.z, sc.x) / (2 * M_PI) + 0.5;
				float v = 0.5 - asin(sc.y) / M_PI;
				Vector2 uv = Vector2(u, v);
				//m_vertices.push_back(uv);

				Vertex vertex = Vertex();
				vertex.Position = v2;
				vertex.Normal = n;
				vertex.TexCoords = uv;
				m_vertexes.push_back(vertex);
				
			}
			if (fV3 == false) {
				m_vertices.push_back(v3);
				//m_vertices.push_back(n);

				//sphere center = (0,0,0)
				Vector3 sc = (v3-sphereCenter).normalize();
				float u = atan2(sc.z, sc.x) / (2 * M_PI) + 0.5;
				float v = 0.5 - asin(sc.y) / M_PI;
				Vector2 uv = Vector2(u, v);
				//m_vertices.push_back(uv);

				Vertex vertex = Vertex();
				vertex.Position = v3;
				vertex.Normal = n;
				vertex.TexCoords = uv;
				m_vertexes.push_back(vertex);
			}

		}

		for (int i = 0; i < N; i++) {


			Triangle t = m_triangles[i];

			Vector3 v1 = t.v1;
			Vector3 v2 = t.v2;
			Vector3 v3 = t.v3;

			V = m_vertices.size();
		
			int index1 = -1, index2 = -1, index3 = -1;
			for (int j = 0; j < V; j++) {

				//if (j % 3 == 1 || j % 3 == 2)continue;

				if (v1.equal(m_vertices[j], epsilon)) {
					index1 = j;
					continue;
				}
				if (v2.equal(m_vertices[j], epsilon)) {
					index2 = j;
					continue;
				}
				if (v3.equal(m_vertices[j], epsilon)) {
					index3 = j;
					continue;
				}

			}

			m_indices.push_back(index1);
			m_indices.push_back(index2);
			m_indices.push_back(index3);

		}

	}

	void updateTriangles() {

		m_triangles.clear();
		int N = m_indices.size();
		for (int i = 0; i < N - 2; i += 3) {
			int index1 = m_indices[i];
			int index2 = m_indices[i + 1];
			int index3 = m_indices[i + 2];

			Triangle tri = Triangle(m_vertices[index1], m_vertices[index2], m_vertices[index3]);

			m_triangles.push_back(tri);
		}
	}

	
};