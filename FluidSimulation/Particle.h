#pragma once

#include <glad\glad.h>
#include <GLFW\glfw3.h>

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

#include <vector>
#include <memory>
#include "Vectors.h"
#include "Triangle.h"
#include "Models.h"
#include <functional>

using namespace std;

class Particle {


public:
	Vector3 pos;
	Vector3 vel;
	Vector3 size;
	Vector3 colour;
	float alpha;
	int shape;
	float lifeSpan;

	int id;
	float pressure=1000;
	float density=1000;
	float nearDensity = 1000;
	float colourField = 0;

	Vector3 totalForce = Vector3(0,0,0);
	Vector3 acceleration = Vector3(0,0,0);
	Vector3 nextForce = Vector3(0, 0, 0);

	Vector3 pVel;
	Vector3 pPos;
	Vector3 pAcc;

	vector<Triangle> mesh;
	vector<int> indices;
	vector<Vector3> vertices;
	vector<Models::Vertex> vertexes;

	vector<int> n_indices;


	Models model;
	int count = 0;

	unsigned int VAO, VBO,EBO,OVBO;
	unsigned int texture;

	bool isRendering = false;
	bool rendered = false;

	bool cX = false;
	bool cY = false;
	bool cZ = false;
	
	Vector4 field;

public:

	Particle() {

		pos = Vector3(0,0,0);
		vel = Vector3(0,0,0);
		acceleration = Vector3(0,0,0);
		totalForce = Vector3(0,0,0);
		count = 0;

		field = Vector4(0,0,0,0);

		pPos = Vector3(0, 0, 0);
		pVel = Vector3(0, 0, 0);
		pAcc = Vector3(0, 0, 0);
		model = Models();
	}

	~Particle() {
	}

	/*Particle(const Particle &p) {
		pos = p.pos;
		vel = p.vel;
		acceleration = p.acceleration;
		totalForce = p.totalForce;
		pressure = p.pressure;
		density = p.density;
		mesh = p.mesh;
		neighbours = p.neighbours;
		id = p.id;
	}
	*/
	Particle(float x, float y, float z, int id) {

		this->pos = Vector3(x,y,z);
		this->id = id;

		vel = Vector3();
		acceleration = Vector3();
		totalForce = Vector3();

		pPos = Vector3(0, 0, 0);
		pVel = Vector3(0, 0, 0);
		pAcc = Vector3(0, 0, 0);
		count = 0;
	}

	/*void clearNeighbours() {
		
		
		neighbours.erase(neighbours.begin(), neighbours.end());
	}*/

	void clearNeighbourIndices() {
		n_indices.clear();
	}

	void setData(vector<Triangle> &mesh) {
		this->mesh = mesh;

	}

	void setData(vector<Vector3> &vertices, vector<int> &indices) {
		this->vertices = vertices;
		this->indices = indices;

	}

	void setData(vector<Vector3> &vertices, vector<int> &indices, vector<Models::Vertex> &vertexes) {
		this->vertices = vertices;
		this->indices = indices;
		this->vertexes = vertexes;

		/*model.m_indices = indices;
		model.m_vertices = vertices;

		model.buildMesh();

		this->vertices = model.m_vertices;
		this->indices = model.m_indices;
		this->vertexes = model.m_vertexes;*/
		
		

	}

	void updateMesh() {
		for (int i = 0; i < mesh.size(); i++) {
			this->mesh[i].move(pos);
		}

		int N = indices.size();
		for (int i = 0; i < N-2; i+=3) {
			vertices[indices[i]] += pos;
			vertices[indices[i+1]] += pos;
			vertices[indices[i+2]] += pos;

			vertexes[indices[i]].Position = vertices[indices[i]];
			vertexes[indices[i+1]].Position = vertices[indices[i+1]];
			vertexes[indices[i+2]].Position = vertices[indices[i+2]];
		}
	}

	void setupBuffersInstanced() {
		glGenVertexArrays(1, &VAO);
		glGenBuffers(1, &OVBO);
		glGenBuffers(1, &VBO);
		glGenBuffers(1, &EBO);

		glBindVertexArray(VAO);
		// load data into vertex buffers
		glBindBuffer(GL_ARRAY_BUFFER, VBO);

		//vector<float> aVertices = getVertices();
		//cout << vertexes.size() << endl;

		glBufferData(GL_ARRAY_BUFFER, vertexes.size() * sizeof(Models::Vertex), &vertexes[0], GL_STREAM_DRAW);
		//glBufferData(GL_ARRAY_BUFFER, sizeof(aVertices.data()), aVertices.data(), GL_STATIC_DRAW);

		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
		glBufferData(GL_ELEMENT_ARRAY_BUFFER, indices.size() * sizeof(unsigned int), &indices[0], GL_STREAM_DRAW);
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

		glEnableVertexAttribArray(3);
		glBindBuffer(GL_ARRAY_BUFFER, OVBO);
		glVertexAttribPointer(3, 3, GL_FLOAT, GL_FALSE, sizeof(Vector3), (void *)0);
		glBindBuffer(GL_ARRAY_BUFFER, 0);
		glVertexAttribDivisor(3, 1);

		//(void*)(3 * sizeof(float))

		glBindVertexArray(0);
	}

	void setupBuffers() {
		glGenVertexArrays(1, &VAO);
		glGenBuffers(1, &VBO);
		glGenBuffers(1, &EBO);

		glBindVertexArray(VAO);
		// load data into vertex buffers
		glBindBuffer(GL_ARRAY_BUFFER, VBO);

		//vector<float> aVertices = getVertices();
		//cout << vertexes.size() << endl;

		glBufferData(GL_ARRAY_BUFFER, vertexes.size()*sizeof(Models::Vertex), &vertexes[0], GL_STREAM_DRAW);
		//glBufferData(GL_ARRAY_BUFFER, sizeof(aVertices.data()), aVertices.data(), GL_STATIC_DRAW);

		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
		glBufferData(GL_ELEMENT_ARRAY_BUFFER, indices.size()*sizeof(unsigned int), &indices[0], GL_STREAM_DRAW);
		//glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(indices.data())*sizeof(unsigned int), &indices.data()[0], GL_STATIC_DRAW);

		// set the vertex attribute pointers
		// Vector3 Positions
		glEnableVertexAttribArray(0);
		glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(Models::Vertex), (void*)0);
		//glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 9*sizeof(float), (void*)0);
		
		// Vector3 Normals
		glEnableVertexAttribArray(1);
		glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, sizeof(Models::Vertex), (void*)offsetof(Models::Vertex,Models::Vertex::Normal));
		//glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 9*sizeof(float), (void*)(3 * sizeof(float)));
		
		// Vector3 (x,y,0) Texture Coords
		glEnableVertexAttribArray(2);
		glVertexAttribPointer(2, 2, GL_FLOAT, GL_FALSE, sizeof(Models::Vertex), (void*)offsetof(Models::Vertex,Models::Vertex::TexCoords));
		//glVertexAttribPointer(2, 3, GL_FLOAT, GL_FALSE, 9*sizeof(float), (void*)(6 * sizeof(float)));
		

		//(void*)(3 * sizeof(float))

		glBindVertexArray(0);
	}

	vector<float> getVertices() {

		int size = vertices.size();

		//float *vs = new float[size*3];

		vector<float> vs;

		for (int i = 0; i < size;i++) {

			vs.push_back(vertices[i].x);
			vs.push_back(vertices[i].y);
			vs.push_back(vertices[i].z);

		}
		
		//cout << vs.size() << endl;
		return vs;
	}

	void renderInstanced() {
		if (rendered == false) {
			setupBuffers();
			rendered = true;
		}

		// draw mesh
		glBindVertexArray(VAO);

		glActiveTexture(GL_TEXTURE0);
		glBindTexture(GL_TEXTURE_2D, texture);

		//works
		glDrawElementsInstanced(GL_TRIANGLES, indices.size(), GL_UNSIGNED_INT, 0, indices.size());
	}

	void render() {
		
		if (rendered == false) {
			setupBuffers();
			rendered = true;
		}

		// draw mesh
		glBindVertexArray(VAO);

		glActiveTexture(GL_TEXTURE0);
		glBindTexture(GL_TEXTURE_2D, texture);
		
		//works
		glDrawElements(GL_TRIANGLES, indices.size(), GL_UNSIGNED_INT, 0);

		//instanced
		//glDrawElementsInstanced(GL_TRIANGLES, indices.size(), GL_UNSIGNED_INT, 0,indices.size());


		//glDrawElements(GL_TRIANGLES, vertices.size()/9, GL_UNSIGNED_INT, 0);
		//glDrawArrays(GL_TRIANGLES, 0, indices.size());
		glBindVertexArray(0);
	}

	Particle(Vector3 pos, Vector3 vel, Vector3 size, Vector3 colour, float alpha, int shape, float lifeSpan) {
		this->pos = pos;
		this->vel = vel;
		this->size = size;
		this->colour = colour;
		this->alpha = alpha;
		this->lifeSpan = lifeSpan;
	}

};