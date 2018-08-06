
// std
#include <algorithm>
#include <utility>
#include <iostream>
// glm
#include <glm.hpp>
#include <gtc/constants.hpp>

// project
#include "shape.hpp"


RayIntersection AABB::intersect(const Ray &ray) {
	RayIntersection intersect;
	glm::vec3 rel_origin = ray.origin - m_center;

	// start magic
	// x
	float rdx_inv = 1 / ray.direction.x;
	float tx1 = (-m_halfsize.x - rel_origin.x) * rdx_inv;
	float tx2 = (m_halfsize.x - rel_origin.x) * rdx_inv;

	float tmin = std::min(tx1, tx2);
	float tmax = std::max(tx1, tx2);

	// y
	float rdy_inv = 1 / ray.direction.y;
	float ty1 = (-m_halfsize.y - rel_origin.y) * rdy_inv;
	float ty2 = (m_halfsize.y - rel_origin.y) * rdy_inv;

	tmin = std::max(tmin, std::min(ty1, ty2));
	tmax = std::min(tmax, std::max(ty1, ty2));

	// z
	float rdz_inv = 1 / ray.direction.z;
	float tz1 = (-m_halfsize.z - rel_origin.z) * rdz_inv;
	float tz2 = (m_halfsize.z - rel_origin.z) * rdz_inv;

	tmin = std::max(tmin, std::min(tz1, tz2));
	tmax = std::min(tmax, std::max(tz1, tz2));

	if (tmax < tmin) return intersect;
	// end magic

	intersect.m_distance = tmin < 0 ? tmax : tmin;
	intersect.m_position = ray.origin + intersect.m_distance * ray.direction;
	intersect.m_valid = tmax >= 0;
	glm::vec3 work_out_a_name_for_it_later = glm::abs((intersect.m_position - m_center) / m_halfsize);
	float max_v = std::max(work_out_a_name_for_it_later[0], std::max(work_out_a_name_for_it_later[1], work_out_a_name_for_it_later[2]));
	intersect.m_normal = glm::normalize(glm::mix(intersect.m_position - m_center, glm::vec3(0), glm::lessThan(work_out_a_name_for_it_later, glm::vec3(max_v))));
	intersect.m_uv_coord = (glm::abs(intersect.m_normal.x) > 0) ?
		glm::vec2(intersect.m_position.y, intersect.m_position.z) :
		glm::vec2(intersect.m_position.x, intersect.m_position.y + intersect.m_position.z);
	intersect.m_shape = this;

	return intersect;
}
float findTValueSphere(float a, float b, float c) {
	float discrim = pow(b, 2) - (4 * a * c);
	float start = -b;

	float first = (start - sqrt(discrim)) / (2 * a);
	float second = (start + sqrt(discrim)) / (2 * a);
	if (first <= 0 && second <= 0) return -1;
	if (first > 0 && second <= 0)return first;
	else if (second > 0 && first <= 0) return second;
	if (first < second) {
		return first;
	}
	else{
		return second;
	}

	//return 0;
}

bool isIntersectingSphere(float a, float b, float c) {
	//Returns if the discriminant is not negative
	return pow(b, 2) - (4 * a * c) >= 0;
}

RayIntersection Sphere::intersect(const Ray &ray) {
	RayIntersection intersect;
	//-------------------------------------------------------------
	// [Assignment 4] :
	// Implement the intersection method for Sphere that returns
	// a RayIntersection object with valid == false if the
	// ray doesn't intersect and true otherwise. If true, then
	// remember to fill out each feild in the object correctly:
	// - m_valid : true if object is itersected
	// - m_distance : distance along the ray of the intersection
	// - m_position : position on the surface of the intersection
	// - m_normal : normal on the surface of the intersection
	// - m_shape : set to "this"
	// - m_uv_coord : texture coordinates (challenge only)
	//-------------------------------------------------------------

	// YOUR CODE GOES HERE
	// ...
	glm::vec3 direction = ray.direction;
	glm::vec3 origin = ray.origin;
	glm::vec3 center = m_center;
	float radius = m_radius;
	float a = glm::dot(direction, direction);
	float b = glm::dot((origin - center) * 2.0f, direction);
	float c = glm::dot((origin - center), (origin - center)) - radius*radius;
	//sets if the intersection is valid or not
	intersect.m_valid = isIntersectingSphere(a, b, c) && findTValueSphere(a, b, c)>0;
	if (intersect.m_valid) {
		//if theres a valid intersection
		float t = findTValueSphere(a, b, c);
		intersect.m_shape = this;

		//gets the position of intersection from the + of the quadratic equation
		intersect.m_position = origin + t * direction;
		intersect.m_normal = glm::normalize(intersect.m_position - center);
		intersect.m_distance = glm::distance(origin, intersect.m_position);
		//intersect.
		glm::vec3 n = intersect.m_normal;
		float u = 0.5+ (atan2(n.z, n.x)) / (2 * glm::pi<float>()) ;
		float v = 0.5-asin(n.y)/glm::pi<float>();
		intersect.m_uv_coord= glm::vec2(u, v);
	}
	return intersect;
}
float findTPlane(glm::vec3 Q,glm::vec3 O,glm::vec3 N,glm::vec3 D) {
	return glm::dot((Q - O), N) / glm::dot(D,N);
}

bool isIntersectingPlane(glm::vec3 Q, glm::vec3 O, glm::vec3 N, glm::vec3 D) {
	if (glm::dot(D, N) != 0) {
		if (findTPlane(Q, O, N, D) > 0) {
			return true;
		}
}
	return false;
}
RayIntersection Plane::intersect(const Ray &ray) {
	RayIntersection intersect;
	glm::vec3 direction = ray.direction;
	glm::vec3 origin = ray.origin;
	float t = -1;
	float denom = glm::dot(ray.direction, glm::normalize(m_normal));
	if (glm::abs(denom) > 1e-4) {
		glm::vec3 dist = m_center - origin;
		t = glm::dot(dist, glm::normalize(m_normal)) / denom;
	}
	glm::vec3 point = origin + t * direction;
	intersect.m_valid = isIntersectingPlane(m_center, origin, m_normal, direction) && t>=0;
	if (intersect.m_valid) {
		//float t = findTPlane(m_center, origin, m_normal, direction);;


		glm::vec3 inter = origin + t * direction;
		intersect.m_position = inter;
		intersect.m_shape = this;
		intersect.m_distance = glm::distance(origin, inter);
		intersect.m_normal= glm::normalize(m_normal);
		if (denom >= 0)
			intersect.m_normal *= (-1.0f);
		float u =  0.5+(point.x- m_center.x)*m_center.x/m_center.x;
		float v = 0.5+(point.z - m_center.z)*m_center.z / m_center.z;

		intersect.m_uv_coord = glm::vec2(u, v);
	}
	return intersect;
}

RayIntersection Disk::intersect(const Ray &ray) {
	RayIntersection intersect;
	glm::vec3 direction = ray.direction;
	glm::vec3 origin = ray.origin;
	float t = -1;
	float denom = glm::dot(ray.direction, glm::normalize(m_normal));
	if (glm::abs(denom) > 1e-4) {
		glm::vec3 dist = m_center - origin;
		t = glm::dot(dist, glm::normalize(m_normal)) / denom;
	}
	glm::vec3 point = origin + t * direction;
	intersect.m_valid = isIntersectingPlane(m_center, origin, m_normal, direction)
		&& t>=0 && glm::distance(point,m_center)<= m_radius;
	if (intersect.m_valid) {
		//float t = findTPlane(m_center, origin, m_normal, direction);

		if (glm::length(point - m_center) < m_radius) {
			intersect.m_position = point;
			intersect.m_shape = this;
			intersect.m_distance = glm::distance(origin, point);
			intersect.m_normal = glm::normalize(m_normal);
			if (denom >= 0)
				intersect.m_normal *= (-1.0f);

			glm::vec3 n = intersect.m_normal;
			float u = atan(n.x) / (2 * glm::pi<float>()) + 0.5;
			float v = 0.5-asin(n.y)/glm::pi<float>();
			intersect.m_uv_coord = glm::vec2(u, v);
		}
		else {
			intersect.m_valid = false;
		}

	}
	return intersect;
}
bool intersectTriangle(glm::vec3 a, glm::vec3 b, glm::vec3 c, glm::vec3 norm) {
	float first = glm::dot(norm,a);
	float second = glm::dot(norm, b);
	float third = glm::dot(norm, c);
	return first > 0 && second > 0 && third > 0;
}
RayIntersection Triangle::intersect(const Ray &ray) {
	RayIntersection intersect;
	glm::vec3 direction = ray.direction;
	glm::vec3 origin = ray.origin;
	float t = -1;
	float denom = glm::dot(ray.direction, glm::normalize(m_normal));
	if (glm::abs(denom) > 1e-4) {
		glm::vec3 dist = m_center - origin;
		t = glm::dot(dist, glm::normalize(m_normal)) / denom;
	}
	glm::vec3 point = origin + t * direction;
	glm::vec3 a = glm::cross(points.at(1) - points.at(0), point - points.at(0));
	glm::vec3 b = glm::cross(points.at(2) - points.at(1), point - points.at(1));
	glm::vec3 c = glm::cross(points.at(0) - points.at(2), point - points.at(2));

	intersect.m_valid = isIntersectingPlane(m_center, origin, m_normal, direction) && intersectTriangle(a,b,c,m_normal) && t >= 0;
	if (intersect.m_valid) {
			intersect.m_position = point;
			intersect.m_shape = this;
			intersect.m_distance = glm::distance(origin, point);
			intersect.m_normal = glm::normalize(m_normal);
			if (denom >= 0)
				intersect.m_normal *= (-1.0f);
			glm::vec3 P = point;
			glm::vec3 A = this->points.at(0), B = this->points.at(1), C = this->points.at(2);
			float areaABC =  glm::length(glm::cross(B - A, C - A));
			float areaACP =  glm::length(glm::cross(B - A, P - A));
			float areaABP =  glm::length(glm::cross(C - A, P - A));
			float u = areaACP / areaABC;
			float v = areaABP / areaABC;
			float w = 1 - u - v;
			float xP = w * A.x + u*B.x + v * c.x;
			float yP = w * A.y + u*B.y + v * c.y;
			intersect.m_uv_coord = glm::vec2(xP, yP);

	}
	return intersect;
}
