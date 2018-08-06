
#pragma once

// glm
#include <glm.hpp>

// project
#include "ray.hpp"
#include "scene.hpp"


class Shape {
public:
	virtual RayIntersection intersect(const Ray &ray) = 0;
};


class AABB : public Shape {
private:
	glm::vec3 m_center;
	glm::vec3 m_halfsize;

public:
	AABB(const glm::vec3 &c, float hs) : m_center(c), m_halfsize(hs) { }
	AABB(const glm::vec3 &c, const glm::vec3 &hs) : m_center(c), m_halfsize(hs) { }
	virtual RayIntersection intersect(const Ray &ray) override;
};


class Sphere : public Shape {
private:
	glm::vec3 m_center;
	float m_radius;

public:
	Sphere(const glm::vec3 &c, float radius) : m_center(c), m_radius(radius) { }
	virtual RayIntersection intersect(const Ray &ray) override;
};

//-------------------------------------------------------------
// [Assignment 4] :
// Implement the following additional Shapes :
// - Plane
// - Disk
// - Triangle
// Follow the pattern shown by AABB and Sphere for implementing
// a class that subclasses Shape making sure that you implement
// the intersect method for each new Shape.
//-------------------------------------------------------------

// YOUR CODE GOES HERE
// ...

class Plane : public Shape {
private:
	glm::vec3 m_center;
	glm::vec3 m_normal;
public:
	Plane(const glm::vec3 &c) : m_center(c) { m_normal = glm::normalize(c); }
	virtual RayIntersection intersect(const Ray &ray) override;

};

class Disk : public Shape {
private:
	glm::vec3 m_center;
	float m_radius;
	glm::vec3 m_normal;
public:
	Disk(const glm::vec3 &c, float radius) : m_center(c), m_radius(radius) { m_normal = glm::normalize(m_center); }
	virtual RayIntersection intersect(const Ray &ray) override;

};

class Triangle : public Shape {
private:
	std::vector<glm::vec3> points;
	glm::vec3 m_normal;
	glm::vec3 m_center;
public:
	Triangle(const glm::vec3 &v0, const glm::vec3 &v1, const glm::vec3 &v2)
	{
		points.push_back(v0);
		points.push_back(v1);
		points.push_back(v2);
		m_normal = glm::cross(points.at(1) - points.at(0), points.at(2) - points.at(0));
		m_center = points.at(0);
	}
	virtual RayIntersection intersect(const Ray &ray) override;
};