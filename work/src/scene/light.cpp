
// glm
#include <gtc/constants.hpp>

// project
#include "light.hpp"


bool DirectionalLight::occluded(Scene *scene, const glm::vec3 &point) const {
	//-------------------------------------------------------------
	// [Assignment 4] :
	// Determine whether the given point is being occluded from
	// this directional light by an object in the scene.
	// Remember that directional lights are "infinitely" far away
	// so any object in the way would cause an occlusion.
	//-------------------------------------------------------------

	// YOUR CODE GOES HERE
	// ...
	Ray r(point, -m_direction);
	RayIntersection ri = scene->intersect(r);

	return ri.m_valid;
}


glm::vec3 DirectionalLight::incidentDirection(const glm::vec3 &) const {
	return m_direction;
}


glm::vec3 DirectionalLight::irradiance(const glm::vec3 &) const {
	return m_irradiance;
}


bool PointLight::occluded(Scene *scene, const glm::vec3 &point) const {
	//-------------------------------------------------------------
	// [Assignment 4] :
	// Determine whether the given point is being occluded from
	// this directional light by an object in the scene.
	// Remember that point lights are somewhere in the scene and
	// an occulsion has to occur somewhere between the light and
	// the given point.
	//-------------------------------------------------------------

	// YOUR CODE GOES HERE
	// ...
	Ray r(point, -incidentDirection(point));
	RayIntersection ri = scene->intersect(r);
	// if valid
	// if distance of interesetion is smaller than distance between point light pos and point the return true

	// otherwise return false
	float distInter = glm::distance(ri.m_position,point);
	float distPoint = glm::distance(point,m_position);
	return distInter < distPoint;


	//return ri.m_valid;
	//return false;
}


glm::vec3 PointLight::incidentDirection(const glm::vec3 &point) const {
	//-------------------------------------------------------------
	// [Assignment 4] :
	// Return the direction of the incoming light (light to point)
	//-------------------------------------------------------------
	return glm::normalize(point - m_position);
	// YOUR CODE GOES HERE
	// ...

//return glm::vec3(0);
}

glm::vec3 PointLight::irradiance(const glm::vec3 &point) const {
	//-------------------------------------------------------------
	// [Assignment 4] :
	// Return the total irradiance on the given point.
	// This requires you to convert the flux of the light into
	// irradiance by dividing it by the surface of the sphere
	// it illuminates. Remember that the surface area increases
	// as the sphere gets bigger, ie. the point is further away.
	//-------------------------------------------------------------

	// YOUR CODE GOES HERE
	// ...
	//radius is distance of point and pos
	float r = glm::distance(point,m_position);
	return m_flux/(4.0f * glm::pi<float>()*(r*r));
	//return glm::vec3(0);
}
