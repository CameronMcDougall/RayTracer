
// glm
#include <gtc/constants.hpp>
#include <gtc/random.hpp>

// std
#include <random>

// project
#include "scene.hpp"
#include "shape.hpp"
#include "light.hpp"
#include "material.hpp"
#include "path_tracer.hpp"



glm::vec3 SimplePathTracer::sampleRay(const Ray &ray, int) {
	// intersect ray with the scene
	RayIntersection intersect = m_scene->intersect(ray);

	// if ray hit something
	if (intersect.m_valid) {

		// simple grey shape shading
		float f = glm::abs(glm::dot(-ray.direction, intersect.m_normal));
		glm::vec3 grey(0.5, 0.5, 0.5);
		return glm::mix(grey / 2.0f, grey, f);
	}

	// no intersection - return background color
	return { 0.3f, 0.3f, 0.4f };
}


glm::vec3 calculateRVec(glm::vec3 L, glm::vec3 N) {
	glm::vec3 nPrime = N*(glm::dot(L, N));
	return nPrime * 2.0f - L;
}



glm::vec3 CorePathTracer::sampleRay(const Ray &ray, int) {
	//-------------------------------------------------------------
	// [Assignment 4] :
	// Implement a PathTracer that calculates the ambient, diffuse
	// and specular, for the given ray in the scene, using the
	// Phong lighting model. Give special consideration to objects
	// that occluded from direct lighting (shadow rays). You do
	// not need to use the depth argument for this implementation.
	//-------------------------------------------------------------

	// YOUR CODE GOES HERE
	// ...
	// intersect ray with the scene
	RayIntersection intersect = m_scene->intersect(ray);

	// if ray hit something
	if (intersect.m_valid) {
		glm::vec3 intersectionColor(0);

		//calculates the starting lighting seen at intersection
		for (int i = 0; i < m_scene->lights().size(); i++) {

			// direct
			if (!m_scene->lights().at(i)->occluded(m_scene, intersect.m_position + intersect.m_normal * 1e-4f)) {
				glm::vec3 irr = m_scene->lights().at(i)->irradiance(intersect.m_position);

				// diffuse
				float costheta = glm::dot(glm::normalize(-m_scene->lights().at(i)->incidentDirection(intersect.m_position)), intersect.m_normal);

				intersectionColor += irr*intersect.m_material->diffuse()*costheta; //diffuse

				// specular
				glm::vec3 R = calculateRVec(glm::normalize(-m_scene->lights().at(i)->incidentDirection(intersect.m_position)), intersect.m_normal);
				glm::vec3 V = -ray.direction / glm::length(ray.direction);
				float dot = pow(glm::dot(R, V), intersect.m_material->shininess());
				if(glm::dot(R,V)>=0)
					intersectionColor += irr * intersect.m_material->specular() * dot;

			}
			glm::vec3 ambient = m_scene->lights().at(i)->ambience() * intersect.m_material->diffuse();
			intersectionColor += ambient;
		}


		return intersectionColor;
	}
	// no intersection - return background color
	return { 0.3f, 0.3f, 0.4f };
}



glm::vec3 CompletionPathTracer::sampleRay(const Ray &ray, int depth) {
	//-------------------------------------------------------------
	// [Assignment 4] :
	// Using the same requirements for the CorePathTracer add in
	// a recursive element to calculate perfect specular reflection.
	// That is compute the reflection ray off your intersection and
	// sample a ray in that direction, using the result to additionally
	// light your object. To make this more realistic you may weight
	// the incoming light by the (1 - (1/shininess)).
	//-------------------------------------------------------------

	// YOUR CODE GOES HERE
	// ...
	RayIntersection intersect = m_scene->intersect(ray);

	// if ray hit something
	if (intersect.m_valid) {
		glm::vec3 intersectionColor(0);
		//glm::vec3 ill(0);
		//calculates the starting lighting seen at intersection
		for (int i = 0; i < m_scene->lights().size(); i++) {
			glm::vec3 ambient = m_scene->lights().at(i)->ambience() * intersect.m_material->diffuse();
			// direct
			intersectionColor += ambient;
			if (!m_scene->lights().at(i)->occluded(m_scene, intersect.m_position + intersect.m_normal * 1e-4f)) {
				glm::vec3 irr = m_scene->lights().at(i)->irradiance(intersect.m_position);
				// diffuse
				float costheta = glm::dot(glm::normalize(-m_scene->lights().at(i)->incidentDirection(intersect.m_position)), intersect.m_normal);
				intersectionColor += irr*intersect.m_material->diffuse()*costheta; //diffuse																   // specular
				glm::vec3 R = calculateRVec(glm::normalize(-m_scene->lights().at(i)->incidentDirection(intersect.m_position)), intersect.m_normal);
				glm::vec3 V = -ray.direction / glm::length(ray.direction);
				float dot = pow(glm::dot(R, V), intersect.m_material->shininess());
				if (glm::dot(R, V) >= 0)
					intersectionColor += irr * intersect.m_material->specular() * dot;
			}



		}
		if (depth >= 0 && intersect.m_material->shininess()&&intersect.m_distance>1e-4) {
			glm::vec3 normal = intersect.m_normal;
			float perp = 2.0 * glm::dot(ray.direction, normal);
			glm::vec3 reflectDir = ray.direction - (normal*perp);
			Ray ray2(intersect.m_position + intersect.m_normal * 1e-4f, reflectDir);
			float m1 = 1 - (1 / intersect.m_material->shininess());
			glm::vec3 mR = sampleRay(ray2, depth - 1) * intersect.m_material->specular();
			intersectionColor += m1* mR + (1 - m1)*intersectionColor;
		}
		return intersectionColor;
	}
	//return { 0,0,0 };
	// no intersection - return background color
	return { 0.3f, 0.3f, 0.4f };
}



glm::vec3 ChallengePathTracer::sampleRay(const Ray &ray, int depth) {
	//-------------------------------------------------------------
	// [Assignment 4] :
	// Implement a PathTracer that calculates the diffuse and
	// specular, for the given ray in the scene, using the
	// Phong lighting model. Give special consideration to objects
	// that occluded from direct lighting (shadow rays).
	// Implement support for textured materials (using a texture
	// for the diffuse portion of the material).
	//
	// EXTRA FOR EXPERTS :
	// Additionally implement indirect diffuse and specular instead
	// of using the ambient lighting term.
	// The diffuse is sampled from the surface hemisphere and the
	// specular is sampled from a cone of the phong lobe (which
	// gives a glossy look). For best results you need to normalize
	// the lighting (see http://www.thetenthplanet.de/archives/255)
	//-------------------------------------------------------------

	// YOUR CODE GOES HERE
	// ...
	RayIntersection intersect = m_scene->intersect(ray);

	// if ray hit something
	if (intersect.m_valid) {
		glm::vec3 intersectionColor(0);

		//glm::vec3 ill(0);
		float x = intersect.m_position.x;
		float y = intersect.m_position.y;
		//calculates the starting lighting seen at intersection
		for (int i = 0; i < m_scene->lights().size(); i++) {
			glm::vec3 ambient = m_scene->lights().at(i)->ambience() * tex.sample(x, y);
			// direct
			intersectionColor += ambient;
			if (!m_scene->lights().at(i)->occluded(m_scene, intersect.m_position + intersect.m_normal * 1e-4f)) {
				glm::vec3 irr = m_scene->lights().at(i)->irradiance(intersect.m_position);
				// diffuse

				float costheta = glm::dot(glm::normalize(-m_scene->lights().at(i)->incidentDirection(intersect.m_position)), intersect.m_normal);
				intersectionColor += irr*tex.sample(x,y) *costheta; //diffuse																   // specular
				glm::vec3 R = calculateRVec(glm::normalize(-m_scene->lights().at(i)->incidentDirection(intersect.m_position)), intersect.m_normal);
				glm::vec3 V = -ray.direction / glm::length(ray.direction);
				float dot = pow(glm::dot(R, V), intersect.m_material->shininess());
				if (glm::dot(R, V) >= 0)
					intersectionColor += irr * intersect.m_material->specular() * dot;
			}



		}
		if (depth >= 0 && intersect.m_material->shininess() && intersect.m_distance>1e-4) {
			glm::vec3 normal = intersect.m_normal;
			float perp = 2.0 * glm::dot(ray.direction, normal);
			glm::vec3 reflectDir = ray.direction - (normal*perp);
			Ray ray2(intersect.m_position + intersect.m_normal * 1e-4f, reflectDir);
			float m1 = 1 - (1 / intersect.m_material->shininess());
			glm::vec3 mR = sampleRay(ray2, depth - 1) * intersect.m_material->specular();
			intersectionColor += m1* mR + (1 - m1)*intersectionColor;
		}
		return intersectionColor;
	}
	//return { 0,0,0 };
	// no intersection - return background color
	return { 0.3f, 0.3f, 0.4f };
}
