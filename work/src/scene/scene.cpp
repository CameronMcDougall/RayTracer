
// std
#include <limits>

// glm
#include <gtc/matrix_transform.hpp>

// project
#include "scene.hpp"
#include "scene_object.hpp"
#include "light.hpp"


RayIntersection Scene::intersect(const Ray &ray) {
	RayIntersection closest_intersect;

	// go through the vector of objects and return the closest intersection
	for (std::shared_ptr<SceneObject> &object : m_objects) {
		RayIntersection intersect = object->intersect(ray);
		if (intersect.m_valid && intersect.m_distance < closest_intersect.m_distance) {
			closest_intersect = intersect;
		}
	}

	return closest_intersect;
}



Scene Scene::simpleScene() {
	std::vector<std::shared_ptr<SceneObject>> objects;
	std::vector<std::shared_ptr<Light>> lights;

	// declare materials
	std::shared_ptr<Material> shiny_red = std::make_shared<Material>(glm::vec3(1, 0, 0), 10, 0.5f, 0);
	std::shared_ptr<Material> green = std::make_shared<Material>(glm::vec3(0, 0.8f, 0), 1.05f, 0.1f, 0);

	// create a box on a sphere
	objects.push_back(std::make_shared<SceneObject>(
		std::make_shared<Sphere>(glm::vec3(0, -2, -10), 1), shiny_red
	));

	objects.push_back(std::make_shared<SceneObject>(
		std::make_shared<AABB>(glm::vec3(0, -3.5f, -10), glm::vec3(3, 0.5, 3)), green
	));

	// one directional light
	lights.push_back(std::make_shared<DirectionalLight>(glm::vec3(-1, -1, -1), glm::vec3(0.5f), glm::vec3(0.05f)));

	return Scene(objects, lights);
}



Scene Scene::lightScene() {
	std::vector<std::shared_ptr<SceneObject>> objects;
	std::vector<std::shared_ptr<Light>> lights;

	// declare materials
	std::shared_ptr<Material> shiny_red = std::make_shared<Material>(glm::vec3(1, 0, 0), 10, 0.5f, 0);
	std::shared_ptr<Material> green = std::make_shared<Material>(glm::vec3(0, 0.8f, 0), 1.05f, 0.1f, 0);

	// create a box on a sphere
	objects.push_back(std::make_shared<SceneObject>(
		std::make_shared<Sphere>(glm::vec3(0, -2, -10), 1), shiny_red
	));

	objects.push_back(std::make_shared<SceneObject>(
		std::make_shared<AABB>(glm::vec3(0, -3.5f, -10), glm::vec3(3, 0.5f, 3)), green
	));

	// wall blocking one of the point lights
	objects.push_back(std::make_shared<SceneObject>(
		std::make_shared<AABB>(glm::vec3(3.5f, 0, -10), glm::vec3(0.5f, 3, 3)), green
	));

	// one directional light
	lights.push_back(std::make_shared<DirectionalLight>(glm::vec3(-1, -1, -1), glm::vec3(0.5f), glm::vec3(0.05f)));

	// two point lights
	lights.push_back(std::make_shared<PointLight>(glm::vec3(-5, 0, -10), glm::vec3(50), glm::vec3(0.05f)));
	lights.push_back(std::make_shared<PointLight>(glm::vec3(5, 0, -10), glm::vec3(50), glm::vec3(0.05f)));

	return Scene(objects, lights);
}



Scene Scene::materialScene() {
	std::vector<std::shared_ptr<SceneObject>> objects;
	std::vector<std::shared_ptr<Light>> lights;

	// declare materials
	std::shared_ptr<Material> green = std::make_shared<Material>(glm::vec3(0, 0.8f, 0), 1.05f, 0.1f, 0);

	// create a grid of materials with varying shininess and specular ratios
	for (int shin = 0; shin <= 10; shin++) {
		for (int spec = 0; spec <= 10; spec++) {
			float shininess = 1 * exp(float(shin));
			float specular_ratio = spec / 10.f;

			std::shared_ptr<Material> m = std::make_shared<Material>(glm::vec3(1, 0, 0), shininess, specular_ratio, 0);

			objects.push_back(std::make_shared<SceneObject>(
				std::make_shared<Sphere>(glm::vec3(5.5f - shin, -2, -5.5f - spec), 0.4), m
			));
		}
	}

	objects.push_back(std::make_shared<SceneObject>(
		std::make_shared<AABB>(glm::vec3(0, -3, -10), glm::vec3(6, 0.5f, 6)), green
	));

	lights.push_back(std::make_shared<DirectionalLight>(glm::vec3(-1, -1, -1), glm::vec3(0.5f), glm::vec3(0.05f)));

	return Scene(objects, lights);
}



Scene Scene::shapeScene() {
	std::vector<std::shared_ptr<SceneObject>> objects;
	std::vector<std::shared_ptr<Light>> lights;

	//-------------------------------------------------------------
	// [Assignment 4] :
	// To show off intersection with shapes create a scene with
	// the different shapes you have implemented.
	// Core:
	//  - Sphere
	// Completion:
	//  - Plane
	//  - Disk
	//  - Triangle
	//-------------------------------------------------------------

	std::shared_ptr<Material> white1 = std::make_shared<Material>(glm::vec3(1), 10, 0.5f, 0);
	std::shared_ptr<Material> white =std::make_shared<Material>(glm::vec3(1, 0, 0), 10, 0.5f, 0);
	std::shared_ptr<Material> white2 = std::make_shared<Material>(glm::vec3(0,0,1), 10, 0.5f, 0);
	objects.push_back(std::make_shared<SceneObject>(
		std::make_shared<AABB>(glm::vec3(-2, -0.5, -5), glm::vec3(0.5)), white
	));


	objects.push_back(std::make_shared<SceneObject>(
	std::make_shared<Sphere>(glm::vec3(2, -0.5, -4), 0.5), white2
	));

	// YOUR CODE GOES HERE
	// ...

	objects.push_back(std::make_shared<SceneObject>(
	std::make_shared<Plane>(glm::vec3(0, -1, 0)), white1
		));
objects.push_back(std::make_shared<SceneObject>(
		std::make_shared<Disk>(glm::vec3(1, -0.3, -5), 0.7), white
	));
	objects.push_back(std::make_shared<SceneObject>(
		std::make_shared<Triangle>(glm::vec3(0, 1, -7), glm::vec3(-0.5, -1, -7), glm::vec3(0.5, -1, -7)), white
		));
		lights.push_back(std::make_shared<PointLight>(glm::vec3(0, 1.5f, -5), glm::vec3(50), glm::vec3(0.05f)));


		lights.push_back(std::make_shared<PointLight>(glm::vec3(-3, 1.5f, -10), glm::vec3(50), glm::vec3(0.05f)));
	lights.push_back(std::make_shared<DirectionalLight>(glm::vec3(-1, -1, -1), glm::vec3(0.5f), glm::vec3(0.05f)));

	return Scene(objects, lights);
}



Scene Scene::cornellBoxScene() {

	std::vector<std::shared_ptr<SceneObject>> objects;
	std::vector<std::shared_ptr<Light>> lights;

	auto white = std::make_shared<Material>(glm::vec3(1), 1.05f, 0.1f, 0);
	auto green = std::make_shared<Material>(glm::vec3(0, 1, 0), 1.05f, 0.1f, 0);
	auto red = std::make_shared<Material>(glm::vec3(1, 0, 0), 1.05f, 0.1f, 0);

	auto gold = std::make_shared<Material>(glm::vec3(1, 1, 0), 50, 0.8f, 1);
	auto silver = std::make_shared<Material>(glm::vec3(1, 1, 1), 1000, 0.8f, 1);
	auto blue = std::make_shared<Material>(glm::vec3(0.5f, 0.5f, 1), 1.1f, 0.1f, 0);


	// box
	//

	// bottom
	objects.push_back(std::make_shared<SceneObject>(
		std::make_shared<AABB>(glm::vec3(0, -3.2f, 0), glm::vec3(3, .2f, 13)), white
	));

	// top
	objects.push_back(std::make_shared<SceneObject>(
		std::make_shared<AABB>(glm::vec3(0, 3.2f, 0), glm::vec3(3, .2f, 13)), white
	));

	// back
	objects.push_back(std::make_shared<SceneObject>(
		std::make_shared<AABB>(glm::vec3(0, 0, -13.2f), glm::vec3(3, 3, .2f)), white
	));

	// front
	objects.push_back(std::make_shared<SceneObject>(
		std::make_shared<AABB>(glm::vec3(0, 0, 13.2f), glm::vec3(3, 3, .2f)), white
	));

	// right
	objects.push_back(std::make_shared<SceneObject>(
		std::make_shared<AABB>(glm::vec3(3.2f, 0, 0), glm::vec3(.2f, 3, 13)), green
	));

	// left
	objects.push_back(std::make_shared<SceneObject>(
		std::make_shared<AABB>(glm::vec3(-3.2f, 0, 0), glm::vec3(.2f, 3, 13)), red
	));


	// spheres
	//
	objects.push_back(std::make_shared<SceneObject>(
		std::make_shared<Sphere>(glm::vec3(1, -2, -7), 1), gold
	));

	objects.push_back(std::make_shared<SceneObject>(
		std::make_shared<Sphere>(glm::vec3(-1.25, -2.25f, -7), .75f), silver
	));

	objects.push_back(std::make_shared<SceneObject>(
		std::make_shared<Sphere>(glm::vec3(0, -1.5, -10), 1.5), blue
	));


	// lights
	lights.push_back(std::make_shared<PointLight>(glm::vec3(0, 2.5f, -10), glm::vec3(50), glm::vec3(0.05f)));

	lights.push_back(std::make_shared<PointLight>(glm::vec3(0, 2.5f, 0), glm::vec3(50), glm::vec3(0.05f)));

	lights.push_back(std::make_shared<PointLight>(glm::vec3(0, 2.5f, 10), glm::vec3(50), glm::vec3(0.05f)));

	return Scene(objects, lights);
}
