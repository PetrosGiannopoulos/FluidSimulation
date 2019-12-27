// FluidSimulation.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include "pch.h"
#include <iostream>

#include "SphGrid.h"
#include "Vectors.h"
#include "Graphics.h"
//#include <GL/glut.h>

using namespace std;

void resizeWrapper(GLFWwindow* window, int width, int height);
void mouseWrapper(GLFWwindow* window, double xpos,double ypos);
void scrollWrapper(GLFWwindow* window, double xpos, double ypos);

//Rendering *rendering;

Graphics *graphics;

int main(int argc, char* argv[])
{
	
	SphGrid grid = SphGrid(1000,800);

	graphics = new Graphics(&grid);

	glfwSetFramebufferSizeCallback(graphics->window, resizeWrapper);
	glfwSetCursorPosCallback(graphics->window, mouseWrapper);
	glfwSetScrollCallback(graphics->window, scrollWrapper);

	graphics->mainLoop();

	graphics->terminate();
	
	/*rendering = new Rendering(&argc,argv, &grid);

	glutDisplayFunc(renderWrapper);
	glutReshapeFunc(resizeWrapper);
	glutIdleFunc(idleWrapper);
	glutMouseFunc(mousePressedWrapper);
	glutMotionFunc(mouseDraggedWrapper);
	glutMainLoop();*/

	//system("pause");

	//delete graphics;

	return 0;
}

void resizeWrapper(GLFWwindow* window,int width, int height) {
	graphics->framebuffer_resize_callback(window, width, height);
}

void mouseWrapper(GLFWwindow* window, double xpos, double ypos) {
	graphics->mouse_callback(window, xpos, ypos);
}
void scrollWrapper(GLFWwindow* window, double xoffset, double yoffset) {
	graphics->scroll_callback(window, xoffset, yoffset);
}

// Run program: Ctrl + F5 or Debug > Start Without Debugging menu
// Debug program: F5 or Debug > Start Debugging menu

// Tips for Getting Started: 
//   1. Use the Solution Explorer window to add/manage files
//   2. Use the Team Explorer window to connect to source control
//   3. Use the Output window to see build output and other messages
//   4. Use the Error List window to view errors
//   5. Go to Project > Add New Item to create new code files, or Project > Add Existing Item to add existing code files to the project
//   6. In the future, to open this project again, go to File > Open > Project and select the .sln file
