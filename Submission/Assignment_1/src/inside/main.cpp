////////////////////////////////////////////////////////////////////////////////
#include <algorithm>
#include <complex>
#include <fstream>
#include <iostream>
#include <numeric>
#include <vector>
////////////////////////////////////////////////////////////////////////////////

typedef std::complex<double> Point;
typedef std::vector<Point> Polygon;

double inline det(const Point &u, const Point &v)
{
	// TODO
	return std::real(u) * std::imag(v) - std::real(v) * std::imag(u);
}

// Return true iff [a,b] intersects [c,d], and store the intersection in ans
bool intersect_segment(const Point &a, const Point &b, const Point &c, const Point &d, Point &ans)
{
	// TODO
	// a + miu * (b - a) is the intersection
	// c + phi * (d - c) is the intersection
	Point l1 = b - a;
	Point l2 = d - c;
	double miu;
	double phi;
	miu = det(c - a, l2) / det(l1, l2);
	phi = det(l1, a - c) / det(l1, l2);
	if (miu >= 0 && miu <= 1 && phi >= 0 && phi <= 1)
	{
		ans = Point(a.real() + miu * l1.real(), a.imag() + miu * l1.imag());
		return true;
	}
	return false;
}

////////////////////////////////////////////////////////////////////////////////

bool is_inside(const Polygon &poly, const Point &query)
{
	// 1. Compute bounding box and set coordinate of a point outside the polygon
	// TODO
	double minX = std::numeric_limits<double>::max();
	double minY = std::numeric_limits<double>::max();
	for (const auto p : poly)
	{
		minX = minX > p.real() ? p.real() : minX;
		minY = minY > p.imag() ? p.imag() : minY;
	}
	Point outside(minX, minY);
	// 2. Cast a ray from the query point to the 'outside' point, count number of intersections
	// TODO
	int count = 0;
	for (int i = 0; i < poly.size(); i++)
	{
		Point intersection;
		if (intersect_segment(poly[i], poly[(i + 1) % poly.size()], outside, query, intersection))
		{
			count++;
		}
	}
	if (count % 2 == 0)
		return false;
	return true;
}

////////////////////////////////////////////////////////////////////////////////

std::vector<Point> load_xyz(const std::string &filename)
{
	std::vector<Point> points;
	std::ifstream in(filename);
	// TODO
	if (!in.is_open())
	{
		throw std::runtime_error("failed to open file " + filename);
	}
	std::string line;
	int size;
	if (getline(in, line))
	{
		size = stoi(line);
	}
	for (int i = 0; i < size; i++)
	{
		if (getline(in, line))
		{
			std::istringstream str(line);
			std::string out;
			str >> out;
			double x = stod(out);
			str >> out;
			double y = stod(out);
			points.push_back(Point(x, y));
		}
	}
	return points;
}

Polygon load_obj(const std::string &filename)
{
	std::ifstream in(filename);
	// TODO
	Polygon poly;
	if (!in.is_open())
	{
		throw std::runtime_error("failed to open file " + filename);
	}
	std::string line;
	while (getline(in, line))
	{
		if (line[0] == 'v')
		{
			std::istringstream str(line);
			std::string out;
			str >> out;
			str >> out;
			double x = stod(out);
			str >> out;
			double y = stod(out);
			poly.push_back(Point(x, y));
		}
	}
	return poly;
}

void save_xyz(const std::string &filename, const std::vector<Point> &points)
{
	// TODO
	std::ofstream out(filename);
	if (!out.is_open())
	{
		throw std::runtime_error("failed to open file " + filename);
	}
	out << std::fixed;
	out << points.size() << "\n";
	for (const auto &v : points)
	{
		out << v.real() << ' ' << v.imag() << " 0\n";
	}
	out << std::endl;
}

////////////////////////////////////////////////////////////////////////////////

int main(int argc, char *argv[])
{
	if (argc <= 3)
	{
		std::cerr << "Usage: " << argv[0] << " points.xyz poly.obj result.xyz" << std::endl;
	}
	std::vector<Point> points = load_xyz(argv[1]);
	Polygon poly = load_obj(argv[2]);
	std::vector<Point> result;
	for (size_t i = 0; i < points.size(); ++i)
	{
		if (is_inside(poly, points[i]))
		{
			result.push_back(points[i]);
		}
	}
	save_xyz(argv[3], result);
	return 0;
}
