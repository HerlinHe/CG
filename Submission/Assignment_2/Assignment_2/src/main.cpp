// C++ include
#include <iostream>
#include <string>
#include <vector>

// Utilities for the Assignment
#include "utils.h"

// Image writing library
#define STB_IMAGE_WRITE_IMPLEMENTATION // Do not include this line twice in your project!
#include "stb_image_write.h"


// Shortcut to avoid Eigen:: everywhere, DO NOT USE IN .h
using namespace Eigen;

double t, gama, beta;

void raytrace_sphere() {
	std::cout << "Simple ray tracer, one sphere with orthographic projection" << std::endl;

	const std::string filename("sphere_orthographic.png");
	MatrixXd C = MatrixXd::Zero(800,800); // Store the color
	MatrixXd A = MatrixXd::Zero(800,800); // Store the alpha mask

	// The camera is orthographic, pointing in the direction -z and covering the unit square (-1,1) in x and y
	Vector3d origin(-1,1,1);
	Vector3d x_displacement(2.0/C.cols(),0,0);
	Vector3d y_displacement(0,-2.0/C.rows(),0);

	// Single light source
	const Vector3d light_position(-1,1,1);

	for (unsigned i=0; i < C.cols(); ++i) {
		for (unsigned j=0; j < C.rows(); ++j) {
			// Prepare the ray
			Vector3d ray_origin = origin + double(i)*x_displacement + double(j)*y_displacement;
			Vector3d ray_direction = RowVector3d(0,0,-1);

			// Intersect with the sphere
			// NOTE: this is a special case of a sphere centered in the origin and for orthographic rays aligned with the z axis
			Vector2d ray_on_xy(ray_origin(0),ray_origin(1));
			const double sphere_radius = 0.9;

			if (ray_on_xy.norm() < sphere_radius) {
				// The ray hit the sphere, compute the exact intersection point
				Vector3d ray_intersection(ray_on_xy(0),ray_on_xy(1),sqrt(sphere_radius*sphere_radius - ray_on_xy.squaredNorm()));

				// Compute normal at the intersection point
				Vector3d ray_normal = ray_intersection.normalized();

				// Simple diffuse model
				C(i,j) = (light_position-ray_intersection).normalized().transpose() * ray_normal;

				// Clamp to zero
				C(i,j) = std::max(C(i,j),0.);

				// Disable the alpha mask for this pixel
				A(i,j) = 1;
			}
		}
	}

	// Save to png
	write_matrix_to_png(C,C,C,A,filename);

}

bool ray_parallelogram_intersection(Vector3d &e, Vector3d &d, Vector3d &a, Vector3d &u, Vector3d &v) {
	Matrix3d M;
	M << -u(0), -v(0), d(0), 
		 -u(1), -v(1), d(1),
	     -u(2), -v(2), d(2);
	Vector3d res = a - e;
	t = M(2,1)*(M(0,0)*res(1)-res(0)*M(1,0))+M(1,1)*(res(0)*M(2,0)-M(0,0)*res(2))+M(0,1)*(M(1,0)*res(2)-res(1)*M(2,0));
	t = -t/M.determinant();
	if (t < 0) return false;
	gama = M(2,2)*(M(0,0)*res(1)-res(0)*M(1,0))+M(1,2)*(res(0)*M(2,0)-M(0,0)*res(2))+M(0,2)*(M(1,0)*res(2)-res(1)*M(2,0));
	gama /= M.determinant();
	if (gama < 0 || gama > 1) return false;
	beta = res(0)*(M(1,1)*M(2,2)-M(1,2)*M(2,1))+res(1)*(M(0,2)*M(2,1)-M(0,1)*M(2,2))+res(2)*(M(0,1)*M(1,2)-M(1,1)*M(0,2));
	beta /= M.determinant();
	if (beta < 0 || beta > 1) return false;
	return true;
}

bool ray_sphere_intersection(Vector3d &e, Vector3d &d, Vector3d &c, double r) {
	double A = d.transpose() * d;
	double B = 2 * d.transpose() * (e - c);
	double C = (e - c).transpose() * (e - c) - r * r;
	double delta = B * B - 4 * A * C;
	if (delta < 0) return false;
	t = -B - sqrt(delta);
	return true;
}

void raytrace_parallelogram() {
	std::cout << "Simple ray tracer, one parallelogram with orthographic projection" << std::endl;

	const std::string filename("plane_orthographic.png");
	MatrixXd C = MatrixXd::Zero(800,800); // Store the color
	MatrixXd A = MatrixXd::Zero(800,800); // Store the alpha mask

	// The camera is orthographic, pointing in the direction -z and covering the unit square (-1,1) in x and y
	Vector3d origin(-1,1,1);
	Vector3d x_displacement(2.0/C.cols(),0,0);
	Vector3d y_displacement(0,-2.0/C.rows(),0);

	// TODO: Parameters of the parallelogram (position of the lower-left corner + two sides)
	Vector3d pgram_origin(-0.9,-0.5,0);
	Vector3d pgram_u(1.4,0,0);
	Vector3d pgram_v(0.4,1,0);

	// Single light source
	const Vector3d light_position(-1,1,1);

	for (unsigned i=0; i < C.cols(); ++i) {
		for (unsigned j=0; j < C.rows(); ++j) {
			// Prepare the ray
			Vector3d ray_origin = origin + double(i)*x_displacement + double(j)*y_displacement;
			Vector3d ray_direction = RowVector3d(0,0,-1);

			// TODO: Check if the ray intersects with the parallelogram
			if (ray_parallelogram_intersection(ray_origin, ray_direction, pgram_origin, pgram_u, pgram_v)) {
				// TODO: The ray hit the parallelogram, compute the exact intersection point
				Vector3d ray_intersection = ray_origin + t * ray_direction;

				// TODO: Compute normal at the intersection point
				Vector3d ray_normal = pgram_u.cross(pgram_v).normalized();

				// Simple diffuse model
				C(i,j) = (light_position-ray_intersection).normalized().transpose() * ray_normal;

				// Clamp to zero
				C(i,j) = std::max(C(i,j),0.);

				// Disable the alpha mask for this pixel
				A(i,j) = 1;
			}
		}
	}

	// Save to png
	write_matrix_to_png(C,C,C,A,filename);
}

void raytrace_perspective() {
	std::cout << "Simple ray tracer, one parallelogram with perspective projection" << std::endl;

	const std::string filename("plane_perspective.png");
	MatrixXd C = MatrixXd::Zero(800,800); // Store the color
	MatrixXd A = MatrixXd::Zero(800,800); // Store the alpha mask

	// The camera is perspective, pointing in the direction -z and covering the unit square (-1,1) in x and y
	Vector3d origin(-1,1,1);
	Vector3d x_displacement(2.0/C.cols(),0,0);
	Vector3d y_displacement(0,-2.0/C.rows(),0);

	// TODO: Parameters of the parallelogram (position of the lower-left corner + two sides)
	Vector3d pgram_origin(-1.8,-0.5,0);
	Vector3d pgram_u(1.4,0,0);
	Vector3d pgram_v(0.4,1,0);

	// define the center of the sphere
	Vector3d sphere_center(0.9,0,0);
	const double sphere_radius = 0.7;

	// Single light source
	const Vector3d light_position(-1,1,1);

	for (unsigned i=0; i < C.cols(); ++i) {
		for (unsigned j=0; j < C.rows(); ++j) {
			// TODO: Prepare the ray (origin point and direction)
			Vector3d camera_origin(0, 0, 2);
			Vector3d ray_origin = origin + double(i)*x_displacement + double(j)*y_displacement;
			Vector3d ray_direction = (ray_origin - camera_origin).normalized();

			// TODO: Check if the ray intersects with the parallelogram
			if (ray_parallelogram_intersection(ray_origin, ray_direction, pgram_origin, pgram_u, pgram_v)) {
				// TODO: The ray hit the parallelogram, compute the exact intersection point
				Vector3d ray_intersection = ray_origin + t * ray_direction;

				// TODO: Compute normal at the intersection point
				Vector3d ray_normal = pgram_u.cross(pgram_v).normalized();

				// Simple diffuse model
				C(i,j) = (light_position-ray_intersection).normalized().transpose() * ray_normal;

				// Clamp to zero
				C(i,j) = std::max(C(i,j),0.);

				// Disable the alpha mask for this pixel
				A(i,j) = 1;
			}

			// check if the ray intersects with the sphere
			if (ray_sphere_intersection(ray_origin, ray_direction, sphere_center, sphere_radius)) {
				Vector3d ray_intersection = ray_origin + t * ray_direction;
				Vector3d ray_normal = (ray_intersection - sphere_center).normalized();

				// Simple diffuse model
				C(i,j) = (light_position-ray_intersection).normalized().transpose() * ray_normal;

				// Clamp to zero
				C(i,j) = std::max(C(i,j),0.);

				// Disable the alpha mask for this pixel
				A(i,j) = 1;
			}
		}
	}

	// Save to png
	write_matrix_to_png(C,C,C,A,filename);
}

void raytrace_shading(){
	std::cout << "Simple ray tracer, one sphere with different shading" << std::endl;

	const std::string filename("shading.png");
	MatrixXd R = MatrixXd::Zero(800,800); // Store the color
	MatrixXd G = MatrixXd::Zero(800,800);
	MatrixXd B = MatrixXd::Zero(800,800);
	MatrixXd A = MatrixXd::Zero(800,800); // Store the alpha mask
	
	// define the RGB color of sphere
	double cr = 0.;
	double cg = 255.;
	double cb = 0.;
	const Vector3d color(cr/255.,cg/255.,cb/255.);

	// Define the sphere
	Vector3d sphere_center(0,0,0);
	const double sphere_radius = 0.9;

	// Define the parallelogram
	Vector3d pgram_origin(-1.8,-1.,0);
	Vector3d pgram_u(2.8,0,0);
	Vector3d pgram_v(0.8,2,0);

	// The camera is perspective, pointing in the direction -z and covering the unit square (-1,1) in x and y
	Vector3d origin(-1,1,1);
	Vector3d x_displacement(2.0/R.cols(),0,0);
	Vector3d y_displacement(0,-2.0/R.rows(),0);

	// Single light source
	const Vector3d light_position(-1,1,1);
	double ambient = 0.1;
	double diffuse_coeff = 0.5;
	double specular_coeff = 0.5;;
	MatrixXd diffuse = MatrixXd::Zero(800, 800);
	MatrixXd specular = MatrixXd::Zero(800, 800);


	for (unsigned i=0; i < R.cols(); ++i) {
		for (unsigned j=0; j < R.rows(); ++j) {
			// Prepare the ray
			Vector3d camera_origin(0,0,2);
			Vector3d ray_origin = origin + double(i)*x_displacement + double(j)*y_displacement;
			Vector3d ray_direction = (ray_origin - camera_origin).normalized();

			// Intersect with the sphere
			if (ray_sphere_intersection(ray_origin, ray_direction, sphere_center, sphere_radius)) {
				// The ray hit the sphere, compute the exact intersection point
				Vector3d ray_intersection = ray_origin + t * ray_direction;

				// Compute normal at the intersection point
				Vector3d ray_normal = ray_intersection.normalized();

				// TODO: Add shading parameter here
				double p = 100;
				diffuse(i,j) = (light_position - ray_intersection).normalized().transpose() * ray_normal;
				Vector3d h = ((light_position - ray_intersection).normalized() + (ray_origin - ray_intersection).normalized()).normalized();
				specular(i,j) = h.transpose() * ray_normal;
				diffuse(i,j) = std::max(diffuse(i,j), 0.);
				specular(i,j) = std::max(specular(i,j),0.);
				specular(i,j) = pow(specular(i,j), p);

				// Final Shading Equation
				double l = ambient + diffuse_coeff * diffuse(i,j) + specular_coeff * specular(i,j);

				// Clamp to zero
				l = std::max(l, 0.);

				// Set the RGB
				R(i,j) = l * color(0);
				G(i,j) = l * color(1);
				B(i,j) = l * color(2);

				// Disable the alpha mask for this pixel
				A(i,j) = 1;
			}

			// Intersect with the parallelogram
			// if (ray_parallelogram_intersection(ray_origin, ray_direction, pgram_origin, pgram_u, pgram_v)) {
			// 	// The ray hit the sphere, compute the exact intersection point
			// 	Vector3d ray_intersection = ray_origin + t * ray_direction;

			// 	// Compute normal at the intersection point
			// 	Vector3d ray_normal = pgram_u.cross(pgram_v).normalized();

			// 	// TODO: Add shading parameter here
			// 	double p = 10;
			// 	diffuse(i,j) = (light_position - ray_intersection).normalized().transpose() * ray_normal;
			// 	Vector3d h = ((light_position - ray_intersection).normalized() + (ray_origin - ray_intersection).normalized()).normalized();
			// 	specular(i,j) = h.transpose() * ray_normal;
			// 	diffuse(i,j) = std::max(diffuse(i,j), 0.);
			// 	specular(i,j) = std::max(specular(i,j),0.);
			// 	specular(i,j) = pow(specular(i,j), p);

			// 	// Final Shading Equation
			// 	double l = ambient + diffuse_coeff * diffuse(i,j) + specular_coeff * specular(i,j);

			// 	// Clamp to zero
			// 	l = std::max(l, 0.);

			// 	// Set the RGB
			// 	R(i,j) = l * color(0);
			// 	G(i,j) = l * color(1);
			// 	B(i,j) = l * color(2);

			// 	// Disable the alpha mask for this pixel
			// 	A(i,j) = 1;
			// }
		}
	}

	// Save to png
	write_matrix_to_png(R,G,B,A,filename);
}


int main() {
	raytrace_sphere();
	raytrace_parallelogram();
	raytrace_perspective();
	raytrace_shading();

	return 0;
}
