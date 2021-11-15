////////////////////////////////////////////////////////////////////////////////
// C++ include
#include <fstream>
#include <iostream>
#include <limits>
#include <memory>
#include <string>
#include <vector>
#include <stack>
#include <math.h>
#include <time.h>

// Eigen for matrix operations
#include <Eigen/Dense>
#include <Eigen/Geometry>

// Image writing library
#define STB_IMAGE_WRITE_IMPLEMENTATION // Do not include this line twice in your project!
#define PI 3.14159265
#include "stb_image_write.h"
#include "utils.h"

// JSON parser library (https://github.com/nlohmann/json)
#include "json.hpp"
using json = nlohmann::json;

// Shortcut to avoid Eigen:: everywhere, DO NOT USE IN .h
using namespace Eigen;

// define the effects enabled
struct CONFIG {
	bool shadow = false; // the shadow part
	bool reflection = false; // reflection part in ray_color
	bool refraction = false; // refraction part in ray_color
	bool dof = false; // depth of filed
} config;

////////////////////////////////////////////////////////////////////////////////
// Define types & classes
////////////////////////////////////////////////////////////////////////////////

struct Ray {
	Vector3d origin;
	Vector3d direction;
	Ray() { }
	Ray(Vector3d o, Vector3d d) : origin(o), direction(d) { }
};

struct Light {
	Vector3d position;
	Vector3d intensity;
};

struct Intersection {
	Vector3d position;
	Vector3d normal;
	double ray_param;
};

struct Camera {
	bool is_perspective;
	Vector3d position;
	double field_of_view; // between 0 and PI
	double focal_length;
	double lens_radius; // for depth of field
};

struct Material {
	Vector3d ambient_color;
	Vector3d diffuse_color;
	Vector3d specular_color;
	double specular_exponent; // Also called "shininess"

	Vector3d reflection_color;
	Vector3d refraction_color;
	double refraction_index;
};

struct Object {
	Material material;
	virtual ~Object() = default; // Classes with virtual methods should have a virtual destructor!
	virtual bool intersect(const Ray &ray, Intersection &hit) = 0;
};

// We use smart pointers to hold objects as this is a virtual class
typedef std::shared_ptr<Object> ObjectPtr;

struct Sphere : public Object {
	Vector3d position;
	double radius;

	virtual ~Sphere() = default;
	virtual bool intersect(const Ray &ray, Intersection &hit) override;
};

struct Parallelogram : public Object {
	Vector3d origin;
	Vector3d u;
	Vector3d v;

	virtual ~Parallelogram() = default;
	virtual bool intersect(const Ray &ray, Intersection &hit) override;
};

struct AABBTree {
	struct Node {
		AlignedBox3d bbox;
		int parent; // Index of the parent node (-1 for root)
		int left; // Index of the left child (-1 for a leaf)
		int right; // Index of the right child (-1 for a leaf)
		int triangle; // Index of the node triangle (-1 for internal nodes)
	};

	std::vector<Node> nodes;
	int root;

	AABBTree() = default; // Default empty constructor
	AABBTree(const MatrixXd &V, const MatrixXi &F); // Build a BVH from an existing mesh
};

struct Mesh : public Object {
	MatrixXd vertices; // n x 3 matrix (n points)
	MatrixXi facets; // m x 3 matrix (m triangles)

	AABBTree bvh;

	Mesh() = default; // Default empty constructor
	Mesh(const std::string &filename);
	virtual ~Mesh() = default;
	virtual bool intersect(const Ray &ray, Intersection &hit) override;
};

struct Scene {
	Vector3d background_color;
	Vector3d ambient_light;

	Camera camera;
	std::vector<Material> materials;
	std::vector<Light> lights;
	std::vector<ObjectPtr> objects;
};

////////////////////////////////////////////////////////////////////////////////

// Read a triangle mesh from an off file
void load_off(const std::string &filename, MatrixXd &V, MatrixXi &F) {
	std::ifstream in(filename);
	std::string token;
	in >> token;
	int nv, nf, ne;
	in >> nv >> nf >> ne;
	V.resize(nv, 3);
	F.resize(nf, 3);
	for (int i = 0; i < nv; ++i) {
		in >> V(i, 0) >> V(i, 1) >> V(i, 2);
	}
	for (int i = 0; i < nf; ++i) {
		int s;
		in >> s >> F(i, 0) >> F(i, 1) >> F(i, 2);
		assert(s == 3);
	}
}

Mesh::Mesh(const std::string &filename) {
	// Load a mesh from a file (assuming this is a .off file), and create a bvh
	load_off(filename, vertices, facets);
	bvh = AABBTree(vertices, facets);
}

////////////////////////////////////////////////////////////////////////////////
// BVH Implementation
////////////////////////////////////////////////////////////////////////////////

// Bounding box of a triangle
AlignedBox3d bbox_triangle(const Vector3d &a, const Vector3d &b, const Vector3d &c) {
	AlignedBox3d box;
	box.extend(a);
	box.extend(b);
	box.extend(c);
	return box;
}

// Recursive function for Top-Down constructor
int split(const MatrixXd &V, const MatrixXi &F, const MatrixXd &Centroids, const std::vector<int> &triangles, AABBTree *tree) {
	// base condition for Top-Down method
	int size = triangles.size();
	// std::cout << "size: " << size << " ";
	if (size == 1) {
		AlignedBox3d root_box;
		for (int i = 0; i < F.cols(); i++) {
			root_box.extend((Vector3d)V.row(F(triangles[0], i)));
		}
		AABBTree::Node root_node;
		root_node.bbox = root_box;
		root_node.triangle = triangles[0];
		tree->nodes.push_back(root_node);
		return tree->nodes.size() - 1;
	}

	// get the box contains all triangles passed in
	AlignedBox3d root_box;
	// Get the biggest box contains all vertices
	for (int i = 0; i < triangles.size(); i++) {
		for (int j = 0; j < F.cols(); j++) {
			root_box.extend((Vector3d)V.row(F(triangles[i], j)));
		}
	}
	// get the longest axis of the root box
	Vector3d sizes = root_box.sizes();
	int axis = 0;
	for (int i = 1; i < V.cols(); i++) {
		if (sizes(i) > sizes(axis)) axis = i;
	}
	// Sort the centroids
	std::vector<std::pair<double, int>> centros;
	for (int i = 0; i < triangles.size(); i++) {
		centros.push_back(std::make_pair(Centroids.row(triangles[i])(axis), triangles[i]));
	}
	std::sort(centros.begin(), centros.end(), [](const std::pair<double, int>& c1, const std::pair<double, int>& c2) {
		return c1.first < c2.first;
	});

	// new sorted triangles
	std::vector<int> new_triangles; // vector stores the index of sorted triangles
	for (std::pair<double, int> p : centros) {
		new_triangles.push_back(p.second);
	}

	std::size_t const half_size = size / 2;
	std::vector<int> left(new_triangles.begin(), new_triangles.begin() + half_size);
	std::vector<int> right(new_triangles.begin() + half_size, new_triangles.end());
	int left_node = split(V, F, Centroids, left, tree);
	int right_node = split(V, F, Centroids, right, tree);
	AABBTree::Node root_node;
	root_node.bbox = root_box;
	root_node.left = left_node;
	root_node.right = right_node;
	tree->nodes[left_node].parent = tree->nodes.size();
	tree->nodes[right_node].parent = tree->nodes[left_node].parent;
	root_node.triangle = -1; // for internal node
	tree->nodes.push_back(root_node);
	return tree->nodes.size() - 1;
}


AABBTree::AABBTree(const MatrixXd &V, const MatrixXi &F) {
	// Compute the centroids of all the triangles in the input mesh
	MatrixXd centroids(F.rows(), V.cols());
	centroids.setZero();
	for (int i = 0; i < F.rows(); ++i) {
		for (int k = 0; k < F.cols(); ++k) {
			centroids.row(i) += V.row(F(i, k));
		}
		centroids.row(i) /= F.cols();
	}

	// TODO (Assignment 3)

	// Method (1): Top-down approach.
	// Split each set of primitives into 2 sets of roughly equal size,
	// based on sorting the centroids along one direction or another.
	
	// get the index of all triangles
	std::vector<int> triangles; // vector stores the index of sorted triangles
	for (int i = 0; i < F.rows(); i++) {
		triangles.push_back(i);
	}

	// Call the recursive function
	clock_t start = clock();
	this->root = split(V, F, centroids, triangles, this);
	clock_t end = clock();
	std::cout << "time for constructing bvh: " << (double)(end - start)/CLOCKS_PER_SEC << std::endl;
	
	// Method (2): Bottom-up approach.
	// Merge nodes 2 by 2, starting from the leaves of the forest, until only 1 tree is left.
}

////////////////////////////////////////////////////////////////////////////////

bool Sphere::intersect(const Ray &ray, Intersection &hit) {
	// TODO (Assignment 2)
	double A = ray.direction.transpose() * ray.direction;
	double B = 2 * ray.direction.transpose() * (ray.origin - this->position);
	double C = (ray.origin - this->position).transpose() * (ray.origin - this->position) - this->radius * this->radius;
	double delta = B * B - 4 * A * C;
	if (delta < 0) return false;
	// Compute the intersection between the ray and the sphere
	// If the ray hits the sphere, set the result of the intersection in the
	// struct 'hit'
	double t1, t2;
	t1 = (-B - sqrt(delta)) / 2 / A;
	t2 = (-B + sqrt(delta)) / 2 / A;
	if (t2 < 0) return false;
	if (t1 < 0) {
		hit.ray_param = t2;
		hit.position = ray.origin + hit.ray_param * ray.direction;
		hit.normal = (this->position - hit.position).normalized();
	}
	else {
		hit.ray_param = t1;
		hit.position = ray.origin + hit.ray_param * ray.direction;
		hit.normal = (hit.position - this->position).normalized();
	}
	return true;
}

bool Parallelogram::intersect(const Ray &ray, Intersection &hit) {
	// TODO (Assignment 2)
	Matrix3d M;
	M << -this->u(0), -this->v(0), ray.direction(0),
		 -this->u(1), -this->v(1), ray.direction(1),
		 -this->u(2), -this->v(2), ray.direction(2);
	Vector3d res = this->origin - ray.origin;
	double t = M(2,1)*(M(0,0)*res(1)-res(0)*M(1,0))+M(1,1)*(res(0)*M(2,0)-M(0,0)*res(2))+M(0,1)*(M(1,0)*res(2)-res(1)*M(2,0));
	t = -t/M.determinant();
	if (t < 0) return false;
	double gama = M(2,2)*(M(0,0)*res(1)-res(0)*M(1,0))+M(1,2)*(res(0)*M(2,0)-M(0,0)*res(2))+M(0,2)*(M(1,0)*res(2)-res(1)*M(2,0));
	gama /= M.determinant();
	if (gama < 0 || gama > 1) return false;
	double beta = res(0)*(M(1,1)*M(2,2)-M(1,2)*M(2,1))+res(1)*(M(0,2)*M(2,1)-M(0,1)*M(2,2))+res(2)*(M(0,1)*M(1,2)-M(1,1)*M(0,2));
	beta /= M.determinant();
	if (beta < 0 || beta > 1) return false;
	hit.ray_param = t;
	hit.position = ray.origin + t * ray.direction;
	hit.normal = this->u.cross(this->v).normalized();
	if (hit.normal.dot(ray.direction) > 0) hit.normal = -hit.normal;
	return true;
}

// -----------------------------------------------------------------------------

bool intersect_triangle(const Ray &ray, const Vector3d &a, const Vector3d &b, const Vector3d &c, Intersection &hit) {
	// TODO (Assignment 3)
	Vector3d u = b - a;
	Vector3d v = c - a;
	Matrix3d M;
	M << -u(0), -v(0), ray.direction(0),
		 -u(1), -v(1), ray.direction(1),
		 -u(2), -v(2), ray.direction(2);
	Vector3d res = a - ray.origin;
	double t = M(2,1)*(M(0,0)*res(1)-res(0)*M(1,0))+M(1,1)*(res(0)*M(2,0)-M(0,0)*res(2))+M(0,1)*(M(1,0)*res(2)-res(1)*M(2,0));
	t = -t/M.determinant();
	if (t < 0) return false;
	double gama = M(2,2)*(M(0,0)*res(1)-res(0)*M(1,0))+M(1,2)*(res(0)*M(2,0)-M(0,0)*res(2))+M(0,2)*(M(1,0)*res(2)-res(1)*M(2,0));
	gama /= M.determinant();
	if (gama < 0 || gama > 1) return false;
	double beta = res(0)*(M(1,1)*M(2,2)-M(1,2)*M(2,1))+res(1)*(M(0,2)*M(2,1)-M(0,1)*M(2,2))+res(2)*(M(0,1)*M(1,2)-M(1,1)*M(0,2));
	beta /= M.determinant();
	if (beta < 0 || beta + gama > 1) return false;
	hit.ray_param = t;
	hit.position = ray.origin + t * ray.direction;
	hit.normal = u.cross(v).normalized();
	if (hit.normal.dot(ray.direction) > 0) hit.normal = -hit.normal; // N should be on the same side with ray
	return true;	
	// Compute whether the ray intersects the given triangle.
	// If you have done the parallelogram case, this should be very similar to it.
}

bool intersect_box(const Ray &ray, const AlignedBox3d &box) {
	// TODO (Assignment 3)

	// brute force
	// int corner_size = 8;
	// Intersection hit;
	// for (int i = 0; i < corner_size - 2; i++) {
	// 	for (int j = i + 1; j < corner_size - 1; j++) {
	// 		for (int k = j + 1; k < corner_size; k++) {
	// 			if (intersect_triangle(ray,
	// 								   box.corner(AlignedBox3d::CornerType(i)),
	// 								   box.corner(AlignedBox3d::CornerType(j)),
	// 								   box.corner(AlignedBox3d::CornerType(k)),
	// 								   hit)) return true;
	// 		}
	// 	}
	// }
	// return false;

	// a better way to do this
	// find the minimal coner and maximal coner
	Vector3d min = box.min();
	Vector3d max = box.max();
	double t_min = std::numeric_limits<double>::min();
	double t_max = std::numeric_limits<double>::max();
	for (int a = 0; a < 3; a++) {
		auto t0 = fmin((min(a) - ray.origin(a)) / ray.direction(a),
					   (max(a) - ray.origin(a)) / ray.direction(a));
		auto t1 = fmax((min(a) - ray.origin(a)) / ray.direction(a),
					   (max(a) - ray.origin(a)) / ray.direction(a));
		t_min = fmax(t0, t_min);
		t_max = fmin(t1, t_max);
		if (t_max <= t_min)
			return false;
	}
	return true;
	// Compute whether the ray intersects the given box.
	// There is no need to set the resulting normal and ray parameter, since
	// we are not testing with the real surface here anyway.
}

bool Mesh::intersect(const Ray &ray, Intersection &closest_hit) {
	// TODO (Assignment 3)

	// Method (1): Traverse every triangle and return the closest hit.
	// closest_hit.ray_param = std::numeric_limits<double>::max();
	// Intersection hit;
	// bool has_intersection = false;
	// for (int i = 0; i < facets.rows(); i++) {
	// 	if (intersect_triangle(ray, vertices.row(facets(i, 0)), vertices.row(facets(i, 1)), vertices.row(facets(i, 2)), hit)) {
	// 		if (hit.ray_param < closest_hit.ray_param) {
	// 			closest_hit.ray_param = hit.ray_param;
	// 			closest_hit.position = hit.position;
	// 			closest_hit.normal = hit.normal;
	// 			has_intersection = true;
	// 		}
	// 	}
	// }

	// Method (2): Traverse the BVH tree and test the intersection with a
	// triangles at the leaf nodes that intersects the input ray.
	AABBTree::Node node = bvh.nodes[bvh.root];
	
	// use two vectors to store all intersected nodes
	std::vector<AABBTree::Node> internal_nodes;
	internal_nodes.push_back(node);
	std::vector<AABBTree::Node> leaf_nodes;
	while (!internal_nodes.empty()) { // still have internal node to verify
		node = internal_nodes.back();
		internal_nodes.pop_back();
		if (node.triangle != -1) { // this is a leaf node
			leaf_nodes.push_back(node);
			continue;
		}
		if (intersect_box(ray, bvh.nodes[node.left].bbox)) {
			// intersect with the left subTree
			internal_nodes.push_back(bvh.nodes[node.left]);
		}
		if (intersect_box(ray, bvh.nodes[node.right].bbox)) {
			// intersect with the right subTree
			internal_nodes.push_back(bvh.nodes[node.right]);
		}
	}

	// search the leaf node
	closest_hit.ray_param = std::numeric_limits<double>::max();
	Intersection hit;
	int triangle;
	bool has_intersection = false;
	for (int i = 0; i < leaf_nodes.size(); i++) {
		triangle = leaf_nodes[i].triangle;
		if (intersect_triangle(ray, 
							   vertices.row(facets(triangle, 0)),
							   vertices.row(facets(triangle, 1)), 
							   vertices.row(facets(triangle, 2)), 
							   hit)) {
			if (hit.ray_param < closest_hit.ray_param) {
				closest_hit.ray_param = hit.ray_param;
				closest_hit.position = hit.position;
				closest_hit.normal = hit.normal;
				has_intersection = true;
			}
		}
	}

	return has_intersection;
}

////////////////////////////////////////////////////////////////////////////////
// Define ray-tracing functions
////////////////////////////////////////////////////////////////////////////////

// Function declaration here (could be put in a header file)
Vector3d ray_color(const Scene &scene, const Ray &ray, const Object &object, const Intersection &hit, int max_bounce);
Object * find_nearest_object(const Scene &scene, const Ray &ray, Intersection &closest_hit);
bool is_light_visible(const Scene &scene, const Ray &ray, const Light &light);
Vector3d shoot_ray(const Scene &scene, const Ray &ray, int max_bounce);

// -----------------------------------------------------------------------------

// function to compute the refract ray
bool refract(const Vector3d &d, const Vector3d &N, double n0, double n1, Vector3d &t) {
	double cos_fi_square = 1 - (n0 * n0 * (1 - d.dot(N) * d.dot(N)))/(n1 * n1);
	if (cos_fi_square < 0) return false;
	t = (n0 * (d - N * d.dot(N)) / n1 - N * sqrt(cos_fi_square)).normalized();
	return true;
}

Vector3d ray_color(const Scene &scene, const Ray &ray, const Object &obj, const Intersection &hit, int max_bounce) {
	// offset for shadow and reflection
	double epsilon = 0.001;

	// Material for hit object
	const Material &mat = obj.material;

	// Ambient light contribution
	Vector3d ambient_color = obj.material.ambient_color.array() * scene.ambient_light.array();

	// Punctual lights contribution (direct lighting)
	Vector3d lights_color(0, 0, 0);
	for (const Light &light : scene.lights) {
		Vector3d Li = (light.position - hit.position).normalized();
		Vector3d N = hit.normal;

		// TODO (Assignment 2, shadow rays)
		if (config.shadow) {
			Ray shadow_ray(hit.position + epsilon * Li, Li); // set some offset so won't intersect with itself
			if (!is_light_visible(scene, shadow_ray, light)) continue;
		}

		// Diffuse contribution
		Vector3d diffuse = mat.diffuse_color * std::max(Li.dot(N), 0.0);

		// TODO (Assignment 2, specular contribution)
		Vector3d h = (Li - ray.direction).normalized();
		Vector3d specular = mat.specular_color * pow(std::max(h.dot(N), 0.0), mat.specular_exponent);

		// Attenuate lights according to the squared distance to the lights
		Vector3d D = light.position - hit.position;
		lights_color += (diffuse + specular).cwiseProduct(light.intensity) /  D.squaredNorm();
	}

	// TODO (Assignment 2, reflected ray)
	Vector3d reflection_color(0, 0, 0);
	//if km is not zero(black)
	if (config.reflection && mat.reflection_color.norm() > 0 && max_bounce > 0) {
		Vector3d r = ray.direction - 2 * ray.direction.dot(hit.normal) * hit.normal;
		Ray reflection_ray(hit.position + epsilon * r, r);
		reflection_color += mat.reflection_color.cwiseProduct(shoot_ray(scene, reflection_ray, max_bounce - 1));
	}

	// TODO (Assignment 2, refracted ray)
	Vector3d refraction_color(0, 0, 0);
	Vector3d t;
	if (config.refraction && mat.refraction_color.norm() > 0 && max_bounce > 0) {
		if (ray.direction.dot(hit.normal) < 0) { // air to object
			if (refract(ray.direction, hit.normal, 1.0, mat.refraction_index, t)) {
				Ray refraction_ray(hit.position + epsilon * t, t);
				refraction_color += mat.refraction_color.cwiseProduct(shoot_ray(scene, refraction_ray, max_bounce - 1));
			}
		}
		else { // object to air
			if (refract(ray.direction, hit.normal, mat.refraction_index, 1.0, t)) {
				Ray refraction_ray(hit.position + epsilon * t, t);
				refraction_color += mat.refraction_color.cwiseProduct(shoot_ray(scene, refraction_ray, max_bounce - 1));
			}
		}
	}

	// Rendering equation
	Vector3d C = ambient_color + lights_color + reflection_color + refraction_color;

	return C;
}

// -----------------------------------------------------------------------------

Object * find_nearest_object(const Scene &scene, const Ray &ray, Intersection &closest_hit) {
	int closest_index = -1;
	// TODO (Assignment 2, find nearest hit)
	closest_hit.ray_param = std::numeric_limits<double>::max();
	Intersection hit;
	for (int i = 0; i < scene.objects.size(); i++) {
		if (scene.objects[i]->intersect(ray, hit)) {
			if (hit.ray_param < closest_hit.ray_param) {
				closest_hit.position = hit.position;
				closest_hit.normal = hit.normal;
				closest_hit.ray_param = hit.ray_param;
				closest_index = i;
			}
		}
	}
	if (closest_index < 0) {
		// Return a NULL pointer
		return nullptr;
	} else {
		// Return a pointer to the hit object. Don't forget to set 'closest_hit' accordingly!
		return scene.objects[closest_index].get();
	}
}

bool is_light_visible(const Scene &scene, const Ray &ray, const Light &light) {
	// TODO (Assignment 2, shadow ray)
	Intersection shadow_ray_hit;
	for (int i = 0; i < scene.objects.size(); i++) {
		if (scene.objects[i].get()->intersect(ray, shadow_ray_hit)) {
			// if the hit point is between the light and ray origin
			if (shadow_ray_hit.ray_param < (light.position - ray.origin).norm()) return false;
		}
	}
	return true;
}

Vector3d shoot_ray(const Scene &scene, const Ray &ray, int max_bounce) {
	Intersection hit;
	if (Object * obj = find_nearest_object(scene, ray, hit)) {
		// 'obj' is not null and points to the object of the scene hit by the ray
		return ray_color(scene, ray, *obj, hit, max_bounce);
	} else {
		// 'obj' is null, we must return the background color
		return scene.background_color;
	}
}

////////////////////////////////////////////////////////////////////////////////

void render_scene(const Scene &scene) {
	std::cout << "Simple ray tracer." << std::endl;

	int w = 640;
	int h = 480;
	MatrixXd R = MatrixXd::Zero(w, h);
	MatrixXd G = MatrixXd::Zero(w, h);
	MatrixXd B = MatrixXd::Zero(w, h);
	MatrixXd A = MatrixXd::Zero(w, h); // Store the alpha mask

	// The camera always points in the direction -z
	// The sensor grid is at a distance 'focal_length' from the camera center,
	// and covers an viewing angle given by 'field_of_view'.
	double aspect_ratio = double(w) / double(h);
	double scale_y = scene.camera.focal_length * tan(scene.camera.field_of_view / 2); // TODO: Stretch the pixel grid by the proper amount here
	double scale_x = scale_y * aspect_ratio; //

	// The pixel grid through which we shoot rays is at a distance 'focal_length'
	// from the sensor, and is scaled from the canonical [-1,1] in order
	// to produce the target field of view.
	Vector3d grid_origin(-scale_x, scale_y, -scene.camera.focal_length);
	Vector3d x_displacement(2.0/w*scale_x, 0, 0);
	Vector3d y_displacement(0, -2.0/h*scale_y, 0);

	// sample camera offset for Depth of field
	int sample_times = 20;
	std::vector<Vector3d> offsets;
	if (config.dof) {
		for (int i = 0; i < sample_times; i++) {
			double radius = (rand() % 100) / 100. * scene.camera.lens_radius;
			double angle = (rand() % 100) / 100. * 2 * PI;
			offsets.push_back(Vector3d(radius * sin(angle), radius * cos(angle), 0));
		}
	}
	
	for (unsigned i = 0; i < w; ++i) {
		std::cout << std::fixed << std::setprecision(2);
		std::cout << "Ray tracing: " << (100.0 * i) / w << "%\r" << std::flush;
		for (unsigned j = 0; j < h; ++j) {
			// TODO (Assignment 2, depth of field)
			Vector3d shift = grid_origin + (i+0.5)*x_displacement + (j+0.5)*y_displacement;

			// Prepare the ray
			Ray ray;
			int max_bounce = 5;
			Vector3d C(0,0,0);

			if (scene.camera.is_perspective) {
				// Perspective camera
				// TODO (Assignment 2, perspective camera)
				if (config.dof) {
					// Depth of field version
					for (int i = 0; i < sample_times; i++) {	
						Vector3d pixel = scene.camera.position + shift;
						ray.origin = scene.camera.position + offsets[i];
						ray.direction = (pixel - ray.origin).normalized();
						C += shoot_ray(scene, ray, max_bounce);
					}
					C /= sample_times;
				}
				else {
					// without depth_of_fields
					Vector3d pixel = scene.camera.position + shift;
					ray.origin = scene.camera.position;
					ray.direction = (pixel - ray.origin).normalized();
					C += shoot_ray(scene, ray, max_bounce);
				}
			} else {
				// Orthographic camera
				ray.origin = scene.camera.position + Vector3d(shift[0], shift[1], 0);
				ray.direction = Vector3d(0, 0, -1);
				C += shoot_ray(scene, ray, max_bounce);
			}

			R(i, j) = C(0);
			G(i, j) = C(1);
			B(i, j) = C(2);
			A(i, j) = 1;
		}
	}

	std::cout << "Ray tracing: 100%  " << std::endl;

	// Save to png
	const std::string filename("raytrace.png");
	write_matrix_to_png(R, G, B, A, filename);
}

////////////////////////////////////////////////////////////////////////////////

Scene load_scene(const std::string &filename) {
	Scene scene;

	// Load json data from scene file
	json data;
	std::ifstream in(filename);
	in >> data;

	// Helper function to read a Vector3d from a json array
	auto read_vec3 = [] (const json &x) {
		return Vector3d(x[0], x[1], x[2]);
	};

	// Read scene info
	scene.background_color = read_vec3(data["Scene"]["Background"]);
	scene.ambient_light = read_vec3(data["Scene"]["Ambient"]);

	// Read camera info
	scene.camera.is_perspective = data["Camera"]["IsPerspective"];
	scene.camera.position = read_vec3(data["Camera"]["Position"]);
	scene.camera.field_of_view = data["Camera"]["FieldOfView"];
	scene.camera.focal_length = data["Camera"]["FocalLength"];
	scene.camera.lens_radius = data["Camera"]["LensRadius"];

	// Read materials
	for (const auto &entry : data["Materials"]) {
		Material mat;
		mat.ambient_color = read_vec3(entry["Ambient"]);
		mat.diffuse_color = read_vec3(entry["Diffuse"]);
		mat.specular_color = read_vec3(entry["Specular"]);
		mat.reflection_color = read_vec3(entry["Mirror"]);
		mat.refraction_color = read_vec3(entry["Refraction"]);
		mat.refraction_index = entry["RefractionIndex"];
		mat.specular_exponent = entry["Shininess"];
		scene.materials.push_back(mat);
	}

	// Read lights
	for (const auto &entry : data["Lights"]) {
		Light light;
		light.position = read_vec3(entry["Position"]);
		light.intensity = read_vec3(entry["Color"]);
		scene.lights.push_back(light);
	}

	// Read objects
	for (const auto &entry : data["Objects"]) {
		ObjectPtr object;
		if (entry["Type"] == "Sphere") {
			auto sphere = std::make_shared<Sphere>();
			sphere->position = read_vec3(entry["Position"]);
			sphere->radius = entry["Radius"];
			object = sphere;
		} else if (entry["Type"] == "Parallelogram") {
			// TODO
		} else if (entry["Type"] == "Mesh") {
			// Load mesh from a file
			std::string filename = std::string(DATA_DIR) + entry["Path"].get<std::string>();
			object = std::make_shared<Mesh>(filename);
		}
		object->material = scene.materials[entry["Material"]];
		scene.objects.push_back(object);
	}

	return scene;
}

////////////////////////////////////////////////////////////////////////////////

int main(int argc, char *argv[]) {
	if (argc < 2) {
		std::cerr << "Usage: " << argv[0] << " scene.json" << std::endl;
		return 1;
	}
	Scene scene = load_scene(argv[1]);
	render_scene(scene);
	return 0;
}
