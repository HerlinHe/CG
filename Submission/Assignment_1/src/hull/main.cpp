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
	return u.real() * v.imag() - v.real() * u.imag();
}

struct Compare
{
	Point p0; // Leftmost point of the poly
	bool operator()(const Point &p1, const Point &p2)
	{
		// TODO
		Point l1 = p1 - p0;
		Point l2 = p2 - p0;
		return std::arg(l1) < std::arg(l2);
	}
};

bool inline salientAngle(Point &a, Point &b, Point &c)
{
	// TODO
	Point l1 = b - a;
	Point l2 = c - b;

	if (det(l1, l2) <= 0)
	{
		return false;
	}
	else
	{
		return true;
	}
}

/////////////////////// /////////////////////////////////////////////////////////

Polygon convex_hull(std::vector<Point> &points)
{
	Compare order;
	// TODO
	Point lower_left = points[0];
	for (std::vector<Point>::iterator it = points.begin(); it != points.end(); it++)
	{
		if ((*it).imag() < lower_left.imag())
		{
			lower_left = *it;
		}
		else if ((*it).imag() == lower_left.imag() && (*it).real() < lower_left.real())
		{
			lower_left = *it;
		}
	}
	order.p0 = lower_left;
	std::sort(points.begin(), points.end(), order);
	Polygon hull;
	// TODO
	// use salientAngle(a, b, c) here
	hull.push_back(points[0]);
	hull.push_back(points[1]);
	for (int i = 2; i < points.size(); i++)
	{
		int size = hull.size();
		while (size > 1)
		{
			Point p1 = hull[size - 2];
			Point p2 = hull[size - 1];
			if (salientAngle(p1, p2, points[i]))
			{
				break;
			}
			else
			{
				hull.pop_back();
				size--;
			}
		}
		hull.push_back(points[i]);
	}
	return hull;
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

void save_obj(const std::string &filename, Polygon &poly)
{
	std::ofstream out(filename);
	if (!out.is_open())
	{
		throw std::runtime_error("failed to open file " + filename);
	}
	out << std::fixed;
	for (const auto &v : poly)
	{
		out << "v " << v.real() << ' ' << v.imag() << " 0\n";
	}
	for (size_t i = 0; i < poly.size(); ++i)
	{
		out << "l " << i + 1 << ' ' << 1 + (i + 1) % poly.size() << "\n";
	}
	out << std::endl;
}

////////////////////////////////////////////////////////////////////////////////

int main(int argc, char *argv[])
{
	if (argc <= 2)
	{
		std::cerr << "Usage: " << argv[0] << " points.xyz output.obj" << std::endl;
	}
	std::vector<Point> points = load_xyz(argv[1]);
	Polygon hull = convex_hull(points);
	save_obj(argv[2], hull);
	return 0;
}
