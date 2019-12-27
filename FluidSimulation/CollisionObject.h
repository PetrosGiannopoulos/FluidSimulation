#pragma once
#define _USE_MATH_DEFINES

class CollisionObject {

public:

	vector<Vector3> m_vertices;
	vector<int> m_indices;
	vector<Models::Vertex> m_vertexes;
	vector<Triangle> mesh;
	vector<Models::TriangleSSBO> ssbo_triangles;

	unsigned int VAO, VBO, EBO, OVBO;
	unsigned int texture;

	float radius;
	Vector3 pos;
	Vector3 vel;
	Vector3 acceleration;
	Vector3 color;

	float mass = 1.5;
	float Volume;
	float density = 1;

	bool rendered = false;
	bool DRT;
	ConstVars vars;

	SphGrid *grid;

	float id;
	bool bounded = false;
	string type = "";

	Vector3 minV;
	Vector3 maxV;

	int collId;
	vector<CollisionObject> *collObjects;
	Vector3 nextForce;

	float cubeSize;
public:

	CollisionObject() {
		radius = 3;
		pos = Vector3(5, 20, 20);
		color = Vector3(1,0,0);
		vars = ConstVars();
		vars.kRebounceCoeff = pow(vars.kRebounceCoeff, 7);

		this->Volume = (4. / 3)*M_PI*pow(radius, 3);
	//	this->mass = density*this->Volume;
		
	}

	CollisionObject(SphGrid *grid) {
		radius = 3;
		pos = Vector3(5, 20, 20);
		color = Vector3(1, 0, 0);
		vars = ConstVars();
		vars.kRebounceCoeff = pow(vars.kRebounceCoeff, 7);

		this->grid = grid;

		this->Volume = (4. / 3)*M_PI*pow(radius,3);
		//this->mass = density*this->Volume;

		//cout << Volume << endl;
		//cout << mass << endl;
		calcAABB();
		nextForce = Vector3(0, 0, 0);

	}

	void setPos(glm::vec3 pos) {
		this->pos = Vector3(pos.x, pos.y, pos.z);
	}

	void simulate() {

		calcAABB();
		acceleration = Vector3(0, 0, 0);

		acceleration += nextForce;
		nextForce = Vector3(0, 0, 0);

		updateFluidCollision2();

		acceleration += vars.kGravity;

		acceleration /= mass;

		updatePosition();
	}

	void setCollObjects(vector<CollisionObject> *collObjects) {
		this->collObjects = collObjects;
	}

	void updateFluidCollision() {

		//logic: for each particle, for each triangle of the sphere , 
		//find the triangle that intersects with a triangle of the 
		//collision sphere and create a force based on that triangle's
		//normal

		
		float pSize = grid->particleSize;
		float collisionStrength = 1;

		float collTriSize = this->mesh.size();
		for (int i = 0; i < grid->nParticles;i++) {

			int triSize = grid->particles[i].mesh.size();
			Vector3 ppos = grid->particles[i].pos;
			
			for (int j = 0; j < triSize;j++) {

				Triangle t = grid->particles[i].mesh[j];

				Vector3 v1 = t.v1*pSize + ppos;
				Vector3 v2 = t.v2*pSize + ppos;
				Vector3 v3 = t.v3*pSize + ppos;

				Triangle postT = Triangle(v1, v2, v3);
				Vector3 pCenter = postT.getCenter();

				float longEdge = postT.getLongEdge();

				float macroDist = (pCenter - pos).length();
				if (macroDist>radius)continue;

				for (int k = 0; k < collTriSize;k++) {

					Triangle collT = mesh[k];

					Vector3 cv1 = collT.v1*pSize + pos;
					Vector3 cv2 = collT.v2*pSize + pos;
					Vector3 cv3 = collT.v3*pSize + pos;

					Triangle postCollT = Triangle(cv1,cv2,cv3);

					float collLongEdge = postCollT.getLongEdge();
					Vector3 collCenter = postCollT.getCenter();

					float dist = (pCenter - collCenter).length();

					if (dist > (longEdge + collLongEdge))continue;

					//check intersection

					bool intersects = postCollT.intersects(postT);

					if (intersects == false)continue;

					Vector3 collisionForce = postCollT.normal()*collisionStrength;

					grid->particles[i].acceleration += collisionForce;
				}
			}


		}

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

			Vector3 v = pos+m_vertices[i];

			if (v.x < minX)minX = v.x;
			if (v.y < minY)minY = v.y;
			if (v.z < minZ)minZ = v.z;

			if (v.x > maxX)maxX = v.x;
			if (v.y > maxY)maxY = v.y;
			if (v.z > maxZ)maxZ = v.z;

			
		}

		minV = Vector3(minX,minY,minZ);
		maxV = Vector3(maxX,maxY,maxZ);

	}

	bool checkCollision(Vector3 ppos, float radius) {

		if (type=="Sphere") {
			float dist = (ppos - pos).length();
			if (dist > (this->radius + radius)) return false;
		}
		else if (type == "Random") {
			//check AABB 
			bool xcondition = ((ppos.x + radius) >= minV.x) && ((ppos.x - radius) <= maxV.x)?true:false;
			bool ycondition = ((ppos.y + radius) >= minV.y) && ((ppos.y - radius) <= maxV.y)?true:false;
			bool zcondition = ((ppos.z + radius) >= minV.z) && ((ppos.z - radius) <= maxV.z)?true:false;
			
			if (xcondition && ycondition && zcondition) return true;
			else return false;
			
		}

		return true;
	}

	void handleCollision(Vector3 ppos, int i, int cType) {
		float collisionStrength = 10;
		

		if (type == "Sphere" && cType == 0) {
			Vector3 collisionDir = ppos - pos;
			Vector3 collisionForce = collisionDir * collisionStrength;

			//simple impluse force 
			grid->particles[i].nextForce += collisionForce * 50;
			if (bounded == false)acceleration -= collisionForce*0.5;
		}
		else if (type == "Random" && cType == 0) {
			
			int N = mesh.size();
			for (int j = 0; j < N;j++) {

				Triangle t = mesh[j];

				Vector3 n = -t.normal();

				t.computePlane();
				Vector3 middle = t.getCenter();
				Vector3 pOut = middle + 10 * n;

				float outDist = (t.a*pOut.x + t.b*pOut.y + t.c*pOut.z + t.d);
				float pDist = (t.a*ppos.x + t.b*ppos.y + t.c*ppos.z + t.d);
				
				if ((outDist < 0 && pDist > 0) || (outDist > 0 && pDist < 0)) {
				
					if (pDist < grid->particleSize*0.5) {

						
						Vector3 collisionDir = -n;
						Vector3 collisionForce = collisionDir * collisionStrength;

						//simple impluse force 
						grid->particles[i].nextForce += collisionForce * 0.05;
						if (bounded == false)acceleration -= collisionForce*0.5;
					}
				}
			}
		}
		else if (type == "Random" && cType == 1) {
			int N = mesh.size();
			for (int j = 0; j < N; j++) {

				Triangle t = mesh[j];

				Vector3 n = -t.normal();
				Vector3 un = -t.normalU();

				t.computePlane();
				Vector3 middle = t.getCenter();
				Vector3 pOut = middle + 10 * n;

				float outDist = (t.a*pOut.x + t.b*pOut.y + t.c*pOut.z + t.d);
				float pDist = (t.a*ppos.x + t.b*ppos.y + t.c*ppos.z + t.d);

				if ((outDist < 0 && pDist > 0) || (outDist > 0 && pDist < 0)) {

					if (pDist < (*collObjects)[i].radius) {
						

						Vector3 collisionDir = -un;
						Vector3 collisionForce = collisionDir * collisionStrength;

						//simple impluse force 
						(*collObjects)[i].nextForce += collisionForce*0.05;
						if (bounded == false)acceleration -= collisionForce*0.5;
					}
				}
			}
		}
	}

	void updateFluidCollision2() {

		for (int i = 0; i < grid->nParticles;i++) {
			Vector3 ppos = grid->particles[i].pos;

			//collision detection
			if (checkCollision(ppos, grid->particleSize*0.5) == false)continue;
			handleCollision(ppos, i, 0);
		}

		for (int i = 0; i < (*collObjects).size(); i++) {
			int cId = (*collObjects)[i].collId;

			if (cId == collId)continue;

			Vector3 cPos = (*collObjects)[i].pos;

			//collision detection
			if (checkCollision(cPos, (*collObjects)[i].radius) == false)continue;
			handleCollision(cPos, i, 1);
		}
	}

	float getCubeSize() {
		int N = m_vertices.size();

		float yMin = FLT_MAX;
		Vector3 cm = Vector3(0, 0, 0);
		for (int i = 0; i < N;i++) {

			Vector3 v = m_vertices[i];
			cm += v;

			if (v.y < yMin)yMin = v.y;
		}
		cm /= (float)N;

		return abs(cm.y-yMin);
	}

	void updatePosition() {

		if (type == "Plane")return;
		
		Vector3 &velocity = vel;
		Vector3 &position = pos;

		velocity += vars.kTimeStep * acceleration;
		position += vars.kTimeStep * velocity;


		
		//default: 0
		if (position.x < 0) {
			position.x = 0;
			velocity.x *= -vars.kRebounceCoeff;

		}
		else if (position.x > 30) {
			position.x = 30;
			velocity.x *= -vars.kRebounceCoeff;

		}

		//default: 0

		if (type == "Sphere" || type=="Random") {

			if (position.y < (radius)) {
				position.y = radius;
				velocity.y *= -vars.kRebounceCoeff;

			}
			else if (position.y > 30) {
				position.y = 30;
				velocity.y *= -vars.kRebounceCoeff;

			}
		}
		else if (type == "Cube") {

			float mul = DRT ? 0.5 : 1;
			if (position.y < (cubeSize*mul)) {
				position.y = cubeSize*mul;
				velocity.y *= -vars.kRebounceCoeff;
			}
			else if (position.y > 30) {
				position.y = 30;
				velocity.y *= -vars.kRebounceCoeff;

			}
		}

		//default: 0
		if (position.z < 0) {
			position.z = 0;
			velocity.z *= -vars.kRebounceCoeff;

		}
		else if (position.z > 30) {
			position.z = 30;
			velocity.z *= -vars.kRebounceCoeff;

		}

		
	}

	void setData(vector<Vector3> vertices, vector<int> indices, vector<Models::Vertex> vertexes) {
		this->m_vertices = vertices;
		this->m_indices = indices;
		this->m_vertexes = vertexes;

		//create UV coords
		if (type == "Sphere") {

			int N = m_vertexes.size();
			for (int i = 0; i < N;i++) {

				//m_vertexes[i].TexCoords = uvSphere((m_vertexes[i].Position-pos).normalize());
			}
		}
		else if (type == "Cube") {
			int N = m_vertexes.size();
			for (int i = 0; i < N; i++) {

				//m_vertexes[i].TexCoords = uvCube((m_vertexes[i].Position - pos).normalize());
			}
		}
	}

	Vector2 uvSphere(Vector3 p) {

		float u = 0.5 + atan2(p.z, p.x) / (2 * PI);
		float v = 0.5 - asin(p.y) / PI;

		return Vector2(u, v);
	}

	Vector2 uvCube(Vector3 p) {

		float absX = fabs(p.x);
		float absY = fabs(p.y);
		float absZ = fabs(p.z);

		int isXPositive = p.x > 0 ? 1 : 0;
		int isYPositive = p.y > 0 ? 1 : 0;
		int isZPositive = p.z > 0 ? 1 : 0;

		float maxAxis, uc, vc;

		if (isXPositive && absX >= absY && absX >= absZ) {
			maxAxis = absX;
			uc = -p.z;
			vc = p.y;
		}
		if (!isXPositive && absX >= absY && absX >= absZ) {
			maxAxis = absX;
			uc = p.z;
			vc = p.y;
		}
		if (isYPositive && absY >= absX && absY >= absZ) {
			maxAxis = absY;
			uc = p.x;
			vc = -p.z;
		}
		if (!isYPositive && absY >= absX && absY >= absZ) {
			maxAxis = absY;
			uc = p.x;
			vc = p.z;
		}
		if (isZPositive && absZ >= absX && absZ >= absY) {
			maxAxis = absZ;
			uc = p.x;
			vc = p.y;
		}
		if (!isZPositive && absZ >= absX && absZ >= absY) {
			maxAxis = absZ;
			uc = -p.x;
			vc = p.y;
		}

		float u = 0.5f * (uc / maxAxis + 1.0f);
		float v = 0.5f * (vc / maxAxis + 1.0f);

		return Vector2(u, v);
	}

	void setTriangleData(vector<Triangle> mesh, vector<Models::TriangleSSBO> ssbo_triangles) {
		this->mesh = mesh;
		this->ssbo_triangles = ssbo_triangles;
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

		glBufferData(GL_ARRAY_BUFFER, m_vertexes.size() * sizeof(Models::Vertex), &m_vertexes[0], GL_STATIC_DRAW);
		//glBufferData(GL_ARRAY_BUFFER, sizeof(aVertices.data()), aVertices.data(), GL_STATIC_DRAW);

		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
		glBufferData(GL_ELEMENT_ARRAY_BUFFER, m_indices.size() * sizeof(unsigned int), &m_indices[0], GL_STATIC_DRAW);
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


		//(void*)(3 * sizeof(float))

		glBindVertexArray(0);
	}

	glm::vec3 getPos() {
		return glm::vec3(pos.x,pos.y,pos.z);
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
		glDrawElements(GL_TRIANGLES, m_indices.size(), GL_UNSIGNED_INT, 0);
		//glDrawArrays(GL_TRIANGLES,0, m_indices.size());
		
		glBindVertexArray(0);
	}

};