////////////////////////////////////////////////////////////////////////////////
// C++ include
#include <fstream>
#include <iostream>
#include <limits>
#include <memory>
#include <string>
#include <vector>
#include <math.h>
#include <gif.h>

// Eigen for matrix operations
#include <Eigen/Dense>

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

struct Scene {
	Vector3d background_color;
	Vector3d ambient_light;

	Camera camera;
	std::vector<Material> materials;
	std::vector<Light> lights;
	std::vector<ObjectPtr> objects;
};

////////////////////////////////////////////////////////////////////////////////

bool Sphere::intersect(const Ray &ray, Intersection &hit) {
	// TODO:
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
	// TODO
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
	return true;
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

		// TODO: Shoot a shadow ray to determine if the light should affect the intersection point
		Ray shadow_ray(hit.position + epsilon * Li, Li); // set some offset so won't intersect with itself
		if (!is_light_visible(scene, shadow_ray, light)) continue;
		
		// Diffuse contribution
		Vector3d diffuse = mat.diffuse_color * std::max(Li.dot(N), 0.0);

		// TODO: Specular contribution
		Vector3d h = (Li - ray.direction).normalized();
		Vector3d specular = mat.specular_color * pow(std::max(h.dot(N), 0.0), mat.specular_exponent);

		// Attenuate lights according to the squared distance to the lights
		Vector3d D = light.position - hit.position;
		lights_color += (diffuse + specular).cwiseProduct(light.intensity) / D.squaredNorm();
	}

	// TODO: Compute the color of the reflected ray and add its contribution to the current point color.
	Vector3d reflection_color(0, 0, 0);
	//if km is not zero(black)
	if (mat.reflection_color.norm() > 0 && max_bounce > 0) {
		Vector3d r = ray.direction - 2 * ray.direction.dot(hit.normal) * hit.normal;
		Ray reflection_ray(hit.position + epsilon * r, r);
		reflection_color += mat.reflection_color.cwiseProduct(shoot_ray(scene, reflection_ray, max_bounce - 1));
	}
	
	// TODO: Compute the color of the refracted ray and add its contribution to the current point color.
	//       Make sure to check for total internal reflection before shooting a new ray.
	Vector3d refraction_color(0, 0, 0);
	Vector3d t;
	if (mat.refraction_color.norm() > 0 && max_bounce > 0) {
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
	// TODO:
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
	// Find the object in the scene that intersects the ray first
	// The function must return 'nullptr' if no object is hit, otherwise it must
	// return a pointer to the hit object, and set the parameters of the argument
	// 'hit' to their expected values.
	if (closest_index < 0) {
		// Return a NULL pointer
		return nullptr;
	} else {
		// Return a pointer to the hit object. Don't forget to set 'closest_hit' accordingly!
		return scene.objects[closest_index].get();
	}
}

bool is_light_visible(const Scene &scene, const Ray &ray, const Light &light) {
	// TODO: Determine if the light is visible here
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

	// sample camera offset
	int sample_times = 20;
	std::vector<Vector3d> offsets;
	for (int i = 0; i < sample_times; i++) {
		double radius = (rand() % 100) / 100. * scene.camera.lens_radius;
		double angle = (rand() % 100) / 100. * 2 * PI;
		offsets.push_back(Vector3d(radius * sin(angle), radius * cos(angle), 0));
	}

	for (unsigned i = 0; i < w; ++i) {
		for (unsigned j = 0; j < h; ++j) {
			// TODO: Implement depth of field
			Vector3d shift = grid_origin + (i+0.5)*x_displacement + (j+0.5)*y_displacement;

			// Prepare the ray
			Ray ray;
			int max_bounce = 5;
			Vector3d C(0,0,0);

			if (scene.camera.is_perspective) {
				// Perspective camera
				// TODO
				for (int i = 0; i < sample_times; i++) {	
					Vector3d pixel = scene.camera.position + shift;
					ray.origin = scene.camera.position + offsets[i];
					ray.direction = (pixel - ray.origin).normalized();
					C += shoot_ray(scene, ray, max_bounce);
				}
				C /= sample_times;
				
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

	// Save to png
	const std::string filename("raytrace.png");
	write_matrix_to_png(R, G, B, A, filename);
}

// animation
void animation(Scene &scene) {
	std::cout << "Amimation ray tracing" << std::endl;

	const char* fileName = "out.gif";
	std::vector<uint8_t> image;

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

	// sample camera offset
	int sample_times = 20;
	std::vector<Vector3d> offsets;
	for (int i = 0; i < sample_times; i++) {
		double radius = (rand() % 100) / 100. * scene.camera.lens_radius;
		double angle = (rand() % 100) / 100. * 2 * PI;
		offsets.push_back(Vector3d(radius * sin(angle), radius * cos(angle), 0));
	}

	int delay = 15; // Milliseconds to wait between frames
	GifWriter g;
	GifBegin(&g, fileName, R.rows(), R.cols(), delay);

	for (unsigned k = 0; k < 20; k++) {
		for (unsigned i = 0; i < w; ++i) {
			for (unsigned j = 0; j < h; ++j) {
				// TODO: Implement depth of field
				Vector3d shift = grid_origin + (i+0.5)*x_displacement + (j+0.5)*y_displacement;

				// Prepare the ray
				Ray ray;
				int max_bounce = 5;
				Vector3d C(0,0,0);

				if (scene.camera.is_perspective) {
					// Perspective camera
					// TODO
					for (int i = 0; i < sample_times; i++) {	
						ray.origin = scene.camera.position + offsets[i];
						ray.direction = shift.normalized();
						C += shoot_ray(scene, ray, max_bounce);
					}
					C /= sample_times;
					
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
		write_matrix_to_uint8(R, G, B, A, image);    
    	GifWriteFrame(&g, image.data(), R.rows(), R.cols(), delay);
		// Matrix3d rotate;
		// double alpha = PI / 60;
		// rotate << cos(alpha), 0, -sin(alpha),
		// 		  0, 1, 0,
		// 		  sin(alpha), 0, cos(alpha);
		//scene.camera.position = rotate * scene.camera.position;
		scene.camera.position(2) += 0.2;
	}
	GifEnd(&g);
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
	render_scene(scene); // simple ray tracing one picture
	// animation(scene); // sample an animation gif
	return 0;
}
