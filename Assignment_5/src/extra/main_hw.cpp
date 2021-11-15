// C++ include
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <math.h>
#include <gif.h>

// Utilities for the Assignment
#include "raster.h"

// Eigen for matrix operations
#include <Eigen/Dense>
#include <Eigen/Geometry>

// Image writing library
#define STB_IMAGE_WRITE_IMPLEMENTATION // Do not include this line twice in your project!
#include "stb_image_write.h"

#define PI 3.14159265

using namespace std;
using namespace Eigen;

struct Mesh {
	MatrixXf vertices; // n x 3 matrix (n points)
	MatrixXi facets; // m x 3 matrix (m triangles)
    MatrixXf vertices_norm;

	Mesh() = default; // Default empty constructor
	Mesh(const std::string &filename);
	virtual ~Mesh() = default;
};

// Read a triangle mesh from an off file
void load_off(const std::string &filename, MatrixXf &V, MatrixXi &F, MatrixXf &N) {
	std::ifstream in(filename);
	string token;
	in >> token;
	int nv, nf, ne;
	in >> nv >> nf >> ne;
	V.resize(nv, 3);
	F.resize(nf, 3);
    N.resize(nv, 3);
	for (int i = 0; i < nv; ++i) {
		in >> V(i, 0) >> V(i, 1) >> V(i, 2);
        N.row(i) << 0, 0, 0;
	}
	for (int i = 0; i < nf; ++i) {
		int s;
		in >> s >> F(i, 0) >> F(i, 1) >> F(i, 2);
		assert(s == 3);
	}
}

Mesh::Mesh(const std::string &filename) {
	// Load a mesh from a file (assuming this is a .off file), and create a bvh
	load_off(filename, vertices, facets, vertices_norm);
}

Mesh load_scene() {
	Mesh scene;

    // Load mesh from a file
    string filename = string(DATA_DIR) + "bunny.off";
    scene = Mesh(filename);

	return scene;
}

void flat_render(
    Program &program,
    Mesh &scene,
    UniformAttributes &uniform,
    vector<VertexAttributes> &vertices,
    Matrix<FrameBufferAttributes,Eigen::Dynamic,Eigen::Dynamic> &frameBuffer) {
    // Render the flat
    uniform.color << 0.8, 0.8, 0.8, 1;
    for (int i = 0; i < scene.facets.rows(); i++) {
        vertices.clear();
        Vector3f vertex1 = scene.vertices.row(scene.facets(i,0));
        Vector3f vertex2 = scene.vertices.row(scene.facets(i,1));
        Vector3f vertex3 = scene.vertices.row(scene.facets(i,2));
        VertexAttributes va1 = VertexAttributes(vertex1(0),vertex1(1),vertex1(2));
        VertexAttributes va2 = VertexAttributes(vertex2(0),vertex2(1),vertex2(2));
        VertexAttributes va3 = VertexAttributes(vertex3(0),vertex3(1),vertex3(2));
        Vector3f facet_norm = (vertex2 - vertex1).cross(vertex3 - vertex1).normalized();
        va1.norm = facet_norm;
        va2.norm = facet_norm;
        va3.norm = facet_norm;
        vertices.push_back(va1);
	    vertices.push_back(va2);
	    vertices.push_back(va3);
        rasterize_triangles(program,uniform,vertices,frameBuffer);
    }
}

void wireframe_render(
    Program &program,
    Mesh &scene,
    UniformAttributes &uniform,
    vector<VertexAttributes> &vertices,
    Matrix<FrameBufferAttributes,Eigen::Dynamic,Eigen::Dynamic> &frameBuffer) {
    // Render the wireframe
    uniform.color << 0, 0, 0, 1;
    double offset = 0.0003; // prevent "z-fighting"
    for (int i = 0; i < scene.facets.rows(); i++) {
        vertices.clear();
        Vector3f vertex1 = scene.vertices.row(scene.facets(i,0));
        Vector3f vertex2 = scene.vertices.row(scene.facets(i,1));
        Vector3f vertex3 = scene.vertices.row(scene.facets(i,2));
        vertices.push_back(VertexAttributes(vertex1(0),vertex1(1),vertex1(2)+offset));
	    vertices.push_back(VertexAttributes(vertex2(0),vertex2(1),vertex2(2)+offset));
        vertices.push_back(VertexAttributes(vertex2(0),vertex2(1),vertex2(2)+offset)); 
	    vertices.push_back(VertexAttributes(vertex3(0),vertex3(1),vertex3(2)+offset));
        vertices.push_back(VertexAttributes(vertex3(0),vertex3(1),vertex3(2)+offset));
        vertices.push_back(VertexAttributes(vertex1(0),vertex1(1),vertex1(2)+offset));
        rasterize_lines(program,uniform,vertices,0.5,frameBuffer);
    }
}

void per_vertex_render(
    Program &program,
    Mesh &scene,
    UniformAttributes &uniform,
    vector<VertexAttributes> &vertices,
    Matrix<FrameBufferAttributes,Eigen::Dynamic,Eigen::Dynamic> &frameBuffer) {
    // Per-vetex shading
    // Compute the norm of each vertex
    for (int i = 0; i < scene.facets.rows(); i++) {
        Vector3f vertex1 = scene.vertices.row(scene.facets(i,0));
        Vector3f vertex2 = scene.vertices.row(scene.facets(i,1));
        Vector3f vertex3 = scene.vertices.row(scene.facets(i,2));
        Vector3f facet_norm = (vertex2 - vertex1).cross(vertex3 - vertex1).normalized();
        scene.vertices_norm.row(scene.facets(i,0)) += facet_norm;
        scene.vertices_norm.row(scene.facets(i,1)) += facet_norm;
        scene.vertices_norm.row(scene.facets(i,2)) += facet_norm;
    }
    // Normalize the norm of vertex
    for (int i = 0; i < scene.vertices_norm.rows(); i++) {
        Vector3f norm = scene.vertices_norm.row(i).normalized();
        scene.vertices_norm(i, 0) = norm(0);
        scene.vertices_norm(i, 1) = norm(1);
        scene.vertices_norm(i, 2) = norm(2);
    }

    // Render the per-vertex scene
    uniform.color << 0.8, 0.8, 0.8, 1;
    for (int i = 0; i < scene.facets.rows(); i++) {
        vertices.clear();
        Vector3f vertex1 = scene.vertices.row(scene.facets(i,0));
        Vector3f vertex2 = scene.vertices.row(scene.facets(i,1));
        Vector3f vertex3 = scene.vertices.row(scene.facets(i,2));
        VertexAttributes va1 = VertexAttributes(vertex1(0),vertex1(1),vertex1(2));
        VertexAttributes va2 = VertexAttributes(vertex2(0),vertex2(1),vertex2(2));
        VertexAttributes va3 = VertexAttributes(vertex3(0),vertex3(1),vertex3(2));
        va1.norm = scene.vertices_norm.row(scene.facets(i,0));
        va2.norm = scene.vertices_norm.row(scene.facets(i,1));
        va3.norm = scene.vertices_norm.row(scene.facets(i,2));
        vertices.push_back(va1);
	    vertices.push_back(va2);
	    vertices.push_back(va3);
        rasterize_triangles(program,uniform,vertices,frameBuffer);
    }
}

void animation(
    Program &program,
    Mesh &scene,
    UniformAttributes &uniform,
    vector<VertexAttributes> &vertices,
    string shading_type,
    int image_size) {
    std::cout << "Amimation" << std::endl;
    const char* fileName = "out.gif";
    std::vector<uint8_t> image;
    
    int delay = 10; // Milliseconds to wait between frames
    int gif_size = 40;
	GifWriter g;
	GifBegin(&g, fileName, image_size, image_size, delay);
    
    // Move the bunny to the center
    Matrix4f translation;
    translation <<
    1, 0, 0, 0.018,
    0, 1, 0, -0.105,
    0, 0, 1, 0,
    0, 0, 0, 1;

    for (int i = 0; i < gif_size; i++) {
        // scale
        Matrix4f scale;
        float scale_coeff = 5+7.0*i/gif_size;
        scale <<
        scale_coeff, 0, 0, 0,
        0, scale_coeff, 0, 0,
        0, 0, scale_coeff, 0,
        0, 0, 0, 1;

        // rotate
        float angle = 360.0*(i+1)/gif_size; // define the rotate angle
        Matrix4f rotation;
        rotation <<
        cos(angle*PI/180), -sin(angle*PI/180), 0, 0,
        sin(angle*PI/180), cos(angle*PI/180), 0, 0,
        0, 0, 1, 0,
        0, 0, 0, 1;

        uniform.view = rotation * scale * translation;
        if (shading_type == "per_vertex") {
            Matrix<FrameBufferAttributes,Eigen::Dynamic,Eigen::Dynamic> frameBuffer(image_size, image_size);
            per_vertex_render(program, scene, uniform, vertices, frameBuffer);
            framebuffer_to_uint8(frameBuffer,image);
            GifWriteFrame(&g, image.data(), frameBuffer.rows(), frameBuffer.cols(), delay);
        }
        if (shading_type == "flat") {
            Matrix<FrameBufferAttributes,Eigen::Dynamic,Eigen::Dynamic> frameBuffer(image_size, image_size);
            flat_render(program, scene, uniform, vertices, frameBuffer);
            wireframe_render(program, scene, uniform, vertices, frameBuffer);
            framebuffer_to_uint8(frameBuffer,image);
            GifWriteFrame(&g, image.data(), frameBuffer.rows(), frameBuffer.cols(), delay);
        }
        if (shading_type == "wireframe") {
            Matrix<FrameBufferAttributes,Eigen::Dynamic,Eigen::Dynamic> frameBuffer(image_size, image_size);
            // Set background color to white
            for (int i = 0; i < frameBuffer.rows(); i++) {
                for (int j = 0; j < frameBuffer.cols(); j++) {
                    frameBuffer(i,j).color << 255, 255, 255, 255;
                }
            }
            wireframe_render(program, scene, uniform, vertices, frameBuffer);
            framebuffer_to_uint8(frameBuffer,image);
            GifWriteFrame(&g, image.data(), frameBuffer.rows(), frameBuffer.cols(), delay);
        }
    }
    GifEnd(&g);
}

int main() 
{
    Mesh scene = load_scene();

	// The Framebuffer storing the image rendered by the rasterizer
	Matrix<FrameBufferAttributes,Eigen::Dynamic,Eigen::Dynamic> frameBuffer(500,500);

    // // Set background color to white
    // for (int i = 0; i < frameBuffer.rows(); i++) {
    //     for (int j = 0; j < frameBuffer.cols(); j++) {
    //         frameBuffer(i,j).color << 255, 255, 255, 255;
    //     }
    // }

	// Global Constants (empty in this example)
	UniformAttributes uniform;

	// Basic rasterization program
	Program program;

	// The vertex shader is the identity
	program.VertexShader = [](const VertexAttributes& va, const UniformAttributes& uniform)
	{
		VertexAttributes out;
		out.position = uniform.view*va.position;
        Vector3f hit;
        hit << va.position(0), va.position(1), va.position(2);
        Vector3f Li = (uniform.light - hit).normalized();
        float diffuse = std::max(Li.dot(va.norm), (float)0.0);
        Vector3f h = (Li - uniform.direction).normalized();
        float specular = pow(std::max(h.dot(va.norm), (float)0.0), uniform.specular_exponent);
        float coeff = diffuse + specular;
        coeff = std::min(coeff, (float)1.0);
        Vector4f color = uniform.color * coeff;
        out.color << color(0), color(1), color(2);
		return out;
	};

	// The fragment shader uses a fixed color
	program.FragmentShader = [](const VertexAttributes& va, const UniformAttributes& uniform)
	{
		FragmentAttributes out(va.color(0),va.color(1),va.color(2),1);
		out.position = va.position;
		return out;
	};

	// The blending shader converts colors between 0 and 1 to uint8
	program.BlendingShader = [](const FragmentAttributes& fa, const FrameBufferAttributes& previous)
	{
		if (fa.position[2] > previous.depth)
		{
			FrameBufferAttributes out(fa.color[0]*255, fa.color[1]*255, fa.color[2]*255, fa.color[3]*255);
			out.depth = fa.position[2];
			return out;
		}
		else
			return previous;
	};

    // Define the transformation

    Matrix4f scale;
    scale <<
    12, 0, 0, 0,
    0, 12, 0, 0,
    0, 0, 12, 0,
    0, 0, 0, 1;

    Matrix4f translation;
    translation <<
    1, 0, 0, 0.016,
    0, 1, 0, -0.11,
    0, 0, 1, 0,
    0, 0, 0, 1;

    uniform.view = scale * translation;
    uniform.light << 0, 0, 10;
    uniform.light_intensity << 16, 16, 16, 16;
    uniform.direction << 0, 0, -1;
    uniform.specular_exponent = 1000;

	// One triangle in the center of the screen
	vector<VertexAttributes> vertices;
    // flat_render(program, scene, uniform, vertices, frameBuffer);
    // wireframe_render(program, scene, uniform, vertices, frameBuffer);
    per_vertex_render(program, scene, uniform, vertices, frameBuffer);

	vector<uint8_t> image;

    // Animation produce gif
    // animation(program, scene, uniform, vertices, "wireframe", frameBuffer.rows());
    // animation(program, scene, uniform, vertices, "flat", frameBuffer.rows());
    animation(program, scene, uniform, vertices, "per_vertex", frameBuffer.rows());

    framebuffer_to_uint8(frameBuffer,image);
	stbi_write_png("triangle.png", frameBuffer.rows(), frameBuffer.cols(), 4, image.data(), frameBuffer.rows()*4);
	return 0;
}