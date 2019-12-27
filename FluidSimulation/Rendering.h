#pragma once

#include <GL/glut.h>
#include "Grid.h"
#include "Models.h"


class Rendering {


public:

	int width, height;
	Grid *grid;
	Models model;

	Vector3 angle;
	float angleMag;

	int mouseX, mouseY;
	int prevX, prevY;
	bool isPressed;

	Vector3 rot;

public:

	Rendering(int *argc,char *argv[], Grid *grid) {

		this->grid = grid;
		//cout << grid->nH << endl;

		glutInit(argc,argv);

		glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE);

		// Define the main window size and initial position
		// ( upper left corner, boundaries included )
		glutInitWindowSize(1000, 800);
		glutInitWindowPosition(50, 50);

		// Create and label the main window
		glutCreateWindow("Thesis");

		// Configure various properties of the OpenGL rendering context
		setup();
		
		
	}

	

	void render() {


		
		// Clean up the colour of the window and the depth buffer
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
		glMatrixMode(GL_MODELVIEW);
		glLoadIdentity();


		glPushMatrix();
		GLfloat lightDirection[] = { 0.0, -1.0,  0.0, 0.0 };
		GLfloat lightPosition[] = { 0.0, 1.0,   0.0,   0.0 };
		GLfloat lightAmbient[] = { 0.5,  0.5,   0.5, 0.0 };
		GLfloat lightDiffuse[] = { 0.8,  0.8,   0.8, 0.0 };
		GLfloat lightSpecular[] = { 0.6, 0.6, 0.6, 0.0 };

		glLightfv(GL_LIGHT0, GL_SPOT_DIRECTION, lightDirection);
		glLightfv(GL_LIGHT0, GL_POSITION, lightPosition);
		glLightfv(GL_LIGHT0, GL_AMBIENT, lightAmbient);
		glLightfv(GL_LIGHT0, GL_DIFFUSE, lightDiffuse);
		
		glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, lightDiffuse);
		glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, lightSpecular);
		glLightfv(GL_LIGHT0, GL_SPECULAR, lightSpecular);
		glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, 20);

		glPopMatrix();



		
		glTranslatef(0, 0, -10);
		//glRotatef(angleMag, angle.x, angle.y, angle.z);
		//rotate x-axis
		glRotatef(rot.x, 1, 0, 0);
		//rotate y-axis
		glRotatef(rot.y, 0, 1, 0);
		//rotate z-axis
		glRotatef(rot.z, 0, 0, 1);


		//floor
		renderFloor();

		//grid->simulate();
		
		grid->isRendering = true;
		glColor3f(0, 0, 1);
		glBegin(GL_TRIANGLES);
		{
			for (int i = 0; i < grid->nParticles; i++) {

				//draw particles as spheres
				Particle pI = grid->particles[i];
				Vector3 pos = pI.pos;

				vector<Triangle> meshI = pI.mesh;
				int meshIN = meshI.size();
				for (int j = 0; j < meshIN; j++) {

					Triangle t = meshI[j];

					Vector3 normT = t.normal();

					glNormal3f(normT.x, normT.y, normT.z);

					glVertex3f(t.v1.x+pos.x, t.v1.y + pos.y, t.v1.z + pos.z);
					glVertex3f(t.v2.x + pos.x, t.v2.y + pos.y, t.v2.z + pos.z);
					glVertex3f(t.v3.x + pos.x, t.v3.y + pos.y, t.v3.z + pos.z);

				}

			}
		}
		glEnd();
		grid->isRendering = false;
		
		
		
		

		glutSwapBuffers();

	}

	void resize(int w, int h) {
		// w and h are window sizes returned by glut
		// define the visible area of the window ( in pixels )

		width = w; height = h;
		if (h == 0) h = 1;

		// Aspect ratio
		float aspect = (float)w / (float)h;

		glViewport(0, 0, w, h);
		//glViewport(0, h/2, w/2, h/2);

		// Setup viewing volume
		glMatrixMode(GL_PROJECTION);
		glLoadIdentity();

		/*gluPerspective(45.0, aspect, 1.0, -200.0);

		glViewport(0, 0, w/2, h/2);
		glMatrixMode(GL_PROJECTION);
		glLoadIdentity();*/

		gluPerspective(45.0, aspect, 1.0, -200.0);
		//         L,      R,       B,      T,       N,       F
		//glOrtho(-100.0f, 100.0f, -100.0f, 100.0f, -200.0f, 200.0f);

		//gluPerspective(135.0, aspect, 1.0, -200.0);
	}

	void idle() {

		
		grid->simulate();

		glutPostRedisplay();

	}

	void mouseDragged(int x, int y) {


		rot.x -= prevY - y;
		rot.y -= prevX - x;

		prevX = x;
		prevY = y;

		
	}

	//Function for handling mouse events
	void mousePressed(int button, int state, int x, int y) {

		//cout << "x: " << x << " y: "<< y <<  endl;
		mouseX = x;
		mouseY = y;

		prevX = x;
		prevY = y;
		
	}


	void setup() {
		// Parameter handling

		// Make models with smooth edges
		glShadeModel(GL_SMOOTH);


		//-------------Task2b---------------
		//*
		// Points represented as circles
		glEnable(GL_POINT_SMOOTH);
		// Set point size
		glPointSize(3.0);

		glEnable(GL_LINE_SMOOTH);
		// Set Line width
		glLineWidth(2.0);
		//*/

		// Blending
		glEnable(GL_BLEND);
		//           incoming             stored
		glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

		//-------------Task7a---------------
		//*
		// Depth test
		glEnable(GL_DEPTH_TEST);
		// Renders a fragment if its z value is less or equal of the stored value
		glDepthFunc(GL_LEQUAL);
		glClearDepth(1.0f);
		//*/

		// Polygon rendering mode
		glEnable(GL_COLOR_MATERIAL);
		glColorMaterial(GL_FRONT, GL_AMBIENT_AND_DIFFUSE);

		// Set up light source
		GLfloat lightDirection[] = { 0, -1,  0, 0.0 };
		GLfloat lightPosition[] = { 0, 1,  0, 0.0 };
		GLfloat lightAmbient[] = { 0.2,  0.2,   0.2, 0.0 };
		GLfloat lightDiffuse[] = { 0.8,  0.8,   0.8, 0.0 };

		glLightfv(GL_LIGHT0, GL_POSITION, lightPosition);
		glLightfv(GL_LIGHT0, GL_AMBIENT, lightAmbient);
		glLightfv(GL_LIGHT0, GL_DIFFUSE, lightDiffuse);
		
		glLightf(GL_LIGHT0, GL_SPOT_CUTOFF, 120.0);
		glLightfv(GL_LIGHT0,GL_SPOT_DIRECTION,lightDirection);
		glLightf(GL_LIGHT0, GL_SPOT_EXPONENT, 1);


		glEnable(GL_LIGHT0);
		glEnable(GL_LIGHTING);

		// polygon rendering mode and material properties
		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

		glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);

		// Black background
		glClearColor(0.0f, 0.0f, 0.0f, 1.0f);

		angle = Vector3(0,1,0);
		angleMag = 0;

		rot = Vector3(0,0,0);

		model = Models();
		model.loadModel("icosahedron.obj");
		model.refineModel(2);
		model.rescaleModel(grid->particleSize);

		grid->setTriangleData(model.m_triangles);
	}

	void renderFloor() {
		glColor3f(0.3, 1, 0.4);

		float fH = -1.5;

		Vector3 v1 = Vector3(-20, fH, 20);
		Vector3 v2 = Vector3(-20,  fH, -20);
		Vector3 v3 = Vector3( 20, fH, 20);
		Triangle t1 = Triangle(v1,v2,v3);

		v1 = Vector3(20, fH, 20);
		v2 = Vector3(-20, fH, -20);
		v3 = Vector3(20, fH, -20);
		Triangle t2 = Triangle(v1,v2,v3);

		glBegin(GL_TRIANGLES);
		{

			Vector3 n1 = t1.normal();
			glNormal3f(n1.x,n1.y,n1.z);
			glVertex3f(t1.v1.x, t1.v1.y, t1.v1.z);
			glVertex3f(t1.v2.x, t1.v2.y, t1.v2.z);
			glVertex3f(t1.v3.x, t1.v3.y, t1.v3.z);

			Vector3 n2 = t2.normal();
			glNormal3f(n2.x, n2.y, n2.z);
			glVertex3f(t2.v1.x, t2.v1.y, t2.v1.z);
			glVertex3f(t2.v2.x, t2.v2.y, t2.v2.z);
			glVertex3f(t2.v3.x, t2.v3.y, t2.v3.z);
		}
		glEnd();
	}
};