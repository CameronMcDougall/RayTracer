
#include <iostream>

// glm
#include <gtc/matrix_transform.hpp>

// project
#include "camera.hpp"
#include "opengl.hpp"

using namespace std;


void Camera::setPositionOrientation(const glm::vec3 &pos, float yaw, float pitch) {
	m_position = pos;
	m_yaw = yaw;
	m_pitch = pitch;

	// update rotation matrix (based on yaw and pitch)
	m_rotation = glm::rotate(glm::mat4(1), m_yaw, glm::vec3(0, 1, 0)) * glm::rotate(glm::mat4(1), m_pitch, glm::vec3(1, 0, 0));
}
float linearInterp(float negativeAng, float positiveAng, float interp){
		return (positiveAng-negativeAng)*interp;
}
Ray Camera::generateRay(const glm::vec2 &pixel) {
	//-------------------------------------------------------------
	// [Assignment 4] :
	// Generate a ray in the scene using the camera position,
	// rotation field of view on the y axis (fovy) and the image
	// size. The pixel is given in image coordinates [0, imagesize]
	// This COULD be done by first creating the ray in ViewSpace
	// then transforming it by the position and rotation to get
	// it into worldspace.
	//-------------------------------------------------------------

	// YOUR CODE GOES HERE
	// ...

	glm::vec3 origin = this->position();
	float angle_y = this->m_fovy/2;
	float angle_x = (m_image_size.x/m_image_size.y)*angle_y;
	glm::vec3 pRight = origin + glm::vec3(tan(angle_x), 0, -1);
	glm::vec3 pLeft= origin + glm::vec3(-tan(angle_x), 0, -1);
	glm::vec3 pTop = origin + glm::vec3( 0,tan(angle_y), -1);
	glm::vec3 pBottom= origin + glm::vec3( 0, -tan(angle_y),-1);
	float x = pLeft.x + (pRight.x-pLeft.x)*(pixel.x/m_image_size.x);
	float y = pBottom.y + (pTop.y-pBottom.y)*(pixel.y/m_image_size.y);
	glm::vec3 pointI = glm::vec3(x,y,-1);
	glm::vec3 D = glm::normalize(pointI - origin);
 	Ray ray;
	ray.origin = origin;
	ray.direction = glm::mat3(m_rotation) * D;
	return ray;
}
