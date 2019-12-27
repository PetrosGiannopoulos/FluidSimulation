#pragma once

#include <glad/glad.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>

#include <vector>

// Defines several possible options for camera movement. Used as abstraction to stay away from window-system specific input methods
enum Camera_Movement {
	FORWARD,
	BACKWARD,
	LEFT,
	RIGHT
};

// Default camera values
const float YAW = -90.0f;
const float PITCH = 0.0f;
//default: 2.5f
const float SPEED = 10.f;
const float SENSITIVITY = 0.1f;
const float ZOOM = 45.0f;
const float PI = 3.14159265359;


// An abstract camera class that processes input and calculates the corresponding Euler Angles, Vectors and Matrices for use in OpenGL
class Camera
{
public:
	// Camera Attributes
	glm::vec3 Position;
	glm::vec3 Front;
	glm::vec3 Up;
	glm::vec3 Right;
	glm::vec3 WorldUp;
	// Euler Angles
	float Yaw;
	float Pitch;
	// Camera options
	float MovementSpeed;
	float MouseSensitivity;
	float Zoom;

	// Constructor with vectors
	Camera(glm::vec3 position = glm::vec3(0.0f, 0.0f, 0.0f), glm::vec3 up = glm::vec3(0.0f, 1.0f, 0.0f), float yaw = YAW, float pitch = PITCH) : Front(glm::vec3(0.0f, 0.0f, -1.0f)), MovementSpeed(SPEED), MouseSensitivity(SENSITIVITY), Zoom(ZOOM)
	{
		Position = position;
		WorldUp = up;
		Yaw = yaw;
		Pitch = pitch;
		updateCameraVectors();
	}
	// Constructor with scalar values
	Camera(float posX, float posY, float posZ, float upX, float upY, float upZ, float yaw, float pitch) : Front(glm::vec3(0.0f, 0.0f, -1.0f)), MovementSpeed(SPEED), MouseSensitivity(SENSITIVITY), Zoom(ZOOM)
	{
		Position = glm::vec3(posX, posY, posZ);
		WorldUp = glm::vec3(upX, upY, upZ);
		Yaw = yaw;
		Pitch = pitch;
		updateCameraVectors();
	}

	// Returns the view matrix calculated using Euler Angles and the LookAt Matrix
	glm::mat4 GetViewMatrix()
	{
		return glm::lookAt(Position, Position + Front, Up);
		//return calculate_lookAt_matrix(Position, Position + Front, Up);
	}

	glm::mat4 GetViewMatrix2()
	{
		//return glm::lookAt(Position, Position + Front, Up);
		return calculate_lookAt_matrix(Position, Front, Up);
	}

	glm::vec3 GetEyeRay(float x, float y, glm::mat4 invProjView) {
		glm::vec4 tmp = glm::vec4(x,y,0.0,1.0);
		tmp = invProjView*tmp;
		tmp /= tmp.w;
		glm::vec3 tmp0 = glm::vec3(tmp.x,tmp.y,tmp.z);

		return tmp0 - Position;
	}

	glm::mat4 getCameraToClipMatrix(float viewAngle, float frustumNear, float frustumFar)
	{
		//Set variables
		glm::mat4 cameraToClipMatrix = glm::mat4(0.0f);

		//Calc frustum scale
		const float degToRad = 3.14159f / 180.0f;
		float viewAngleRad = viewAngle * degToRad;
		float frustumScale = 1.0f / tan(viewAngleRad / 2.0f);

		//Create matrix
		cameraToClipMatrix[0].x = frustumScale;
		cameraToClipMatrix[1].y = frustumScale;
		cameraToClipMatrix[2].z = (frustumFar + frustumNear) / (frustumNear - frustumFar);
		cameraToClipMatrix[2].w = -1.0f;
		cameraToClipMatrix[3].z = (2.0f*frustumFar*frustumNear) / (frustumNear - frustumFar);

		return cameraToClipMatrix;
	}

	// Custom implementation of the LookAt function
	glm::mat4 calculate_lookAt_matrix(glm::vec3 position, glm::vec3 target, glm::vec3 worldUp)
	{
		// 1. Position = known
		// 2. Calculate cameraDirection
		glm::vec3 zaxis = glm::normalize(position - target);
		// 3. Get positive right axis vector
		glm::vec3 xaxis = glm::normalize(glm::cross(glm::normalize(worldUp), zaxis));
		// 4. Calculate camera up vector
		glm::vec3 yaxis = glm::cross(zaxis, xaxis);

		// Create translation and rotation matrix
		// In glm we access elements as mat[col][row] due to column-major layout
		glm::mat4 translation; // Identity matrix by default
		translation[3][0] = -position.x; // Third column, first row
		translation[3][1] = -position.y;
		translation[3][2] = -position.z;
		glm::mat4 rotation;
		rotation[0][0] = xaxis.x; // First column, first row
		rotation[1][0] = xaxis.y;
		rotation[2][0] = xaxis.z;
		rotation[0][1] = yaxis.x; // First column, second row
		rotation[1][1] = yaxis.y;
		rotation[2][1] = yaxis.z;
		rotation[0][2] = zaxis.x; // First column, third row
		rotation[1][2] = zaxis.y;
		rotation[2][2] = zaxis.z;

		// Return lookAt matrix as combination of translation and rotation matrix
		return rotation * translation; // Remember to read from right to left (first translation then rotation)
	}

	glm::mat4 lookAt(glm::vec3 position, glm::vec3 target, glm::vec3 worldUp)
	{
		// Camera Direction.
		glm::vec3 zAxis = glm::normalize(position - target);
		// Positive right axis vector
		glm::vec3 xAxis = glm::normalize(glm::cross(glm::normalize(worldUp), zAxis));
		// Camera up vector.
		glm::vec3 yAxis = glm::cross(zAxis, xAxis);

		glm::mat4 data = glm::mat4(0.0);

		glm::value_ptr(data)[0] = xAxis.x;
		glm::value_ptr(data)[1] = yAxis.x;
		glm::value_ptr(data)[2] = zAxis.x;

		glm::value_ptr(data)[4] = xAxis.y;
		glm::value_ptr(data)[5] = yAxis.y;
		glm::value_ptr(data)[6] = zAxis.y;

		glm::value_ptr(data)[8] = xAxis.z;
		glm::value_ptr(data)[9] = yAxis.z;
		glm::value_ptr(data)[10] = zAxis.z;

		glm::value_ptr(data)[12] = -glm::dot(xAxis,position);
		glm::value_ptr(data)[13] = -glm::dot(yAxis,position);
		glm::value_ptr(data)[14] = -glm::dot(zAxis,position);

		return data;
	}

	glm::mat4 setFrustrum(float l, float r, float b, float t, float n, float f)
	{
		glm::mat4 data = glm::mat4(0.0);
		glm::value_ptr(data)[0] = (2 * n) / (r - l);
		glm::value_ptr(data)[5] = (2 * n) / (t - b);
		glm::value_ptr(data)[8] = (r + l) / (r - l);
		glm::value_ptr(data)[9] = (t + b) / (t - b);
		glm::value_ptr(data)[10] = -((f + n) / (f - n));
		glm::value_ptr(data)[11] = -1;
		glm::value_ptr(data)[14] = -((2 * n * f) / (f - n));
		glm::value_ptr(data)[15] = 0;

		return data;
	}

	glm::mat4 perspectiveProjection(float fov, float aspectRatio, float nearPlane, float farPlane)
	{
		float tangent = tanf((fov * 0.5) * (PI / 180));
		float height = nearPlane * tangent;
		float width = height * aspectRatio;

		glm::mat4 data = setFrustrum(-width, width, -height, height, nearPlane, farPlane);

		return data;
	}


	// Processes input received from any keyboard-like input system. Accepts input parameter in the form of camera defined ENUM (to abstract it from windowing systems)
	void ProcessKeyboard(Camera_Movement direction, float deltaTime)
	{
		float velocity = MovementSpeed * deltaTime;
		if (direction == FORWARD)
			Position += Front * velocity;
		if (direction == BACKWARD)
			Position -= Front * velocity;
		if (direction == LEFT)
			Position -= Right * velocity;
		if (direction == RIGHT)
			Position += Right * velocity;
		// make sure the user stays at the ground level
		//Position.y = 0.0f; // <-- this one-liner keeps the user at the ground level (xz plane)
	}

	// Processes input received from a mouse input system. Expects the offset value in both the x and y direction.
	void ProcessMouseMovement(float xoffset, float yoffset, GLboolean constrainPitch = true)
	{
		xoffset *= MouseSensitivity;
		yoffset *= MouseSensitivity;

		Yaw += xoffset;
		Pitch += yoffset;

		// Make sure that when pitch is out of bounds, screen doesn't get flipped
		if (constrainPitch)
		{
			if (Pitch > 89.0f)
				Pitch = 89.0f;
			if (Pitch < -89.0f)
				Pitch = -89.0f;
		}

		// Update Front, Right and Up Vectors using the updated Euler angles
		updateCameraVectors();
	}

	// Processes input received from a mouse scroll-wheel event. Only requires input on the vertical wheel-axis
	void ProcessMouseScroll(float yoffset)
	{
		if (Zoom >= 1.0f && Zoom <= 45.0f)
			Zoom -= yoffset;
		if (Zoom <= 1.0f)
			Zoom = 1.0f;
		if (Zoom >= 45.0f)
			Zoom = 45.0f;
	}

public:
	// Calculates the front vector from the Camera's (updated) Euler Angles
	void updateCameraVectors()
	{
		// Calculate the new Front vector
		glm::vec3 front;
		front.x = cos(glm::radians(Yaw)) * cos(glm::radians(Pitch));
		front.y = sin(glm::radians(Pitch));
		front.z = sin(glm::radians(Yaw)) * cos(glm::radians(Pitch));
		Front = glm::normalize(front);
		// Also re-calculate the Right and Up vector
		Right = glm::normalize(glm::cross(Front, WorldUp));  // Normalize the vectors, because their length gets closer to 0 the more you look up or down which results in slower movement.
		Up = glm::normalize(glm::cross(Right, Front));
	}
};
