// based on NAMD.
#ifndef ___CLASS_LATTICE
#define ___CLASS_LATTICE

#include <Eigen/Core>
#include <Eigen/Geometry>
#include "common.hpp"
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <sstream>

class Lattice
{
	private:
		std::string filename;
		Eigen::Vector3d o, a1, a2, a3, g1, g2, g3, b1, b2, b3;

	public:
		Lattice(std::string filename) 
		{
			this -> filename = filename;
		}

		bool set()
		{
			std::ifstream fi(filename.c_str());
			if (!fi)
			{
				std::cerr << "xsc file not exists.\n";
				return false;
			}
			std::string s;
			while (std::getline(fi, s)) {
			if (s[0] == '#') continue;
			std::istringstream is(s);
			int itmp;
			is >> itmp 
				>> a1.x() >> a1.y() >> a1.z()
				>> a2.x() >> a2.y() >> a2.z()
				>> a3.x() >> a3.y() >> a3.z()	
				>> o.x()  >> o.y()  >> o.z();
			}
			recalculate();
			return true;
		}

	Lattice () {}
	Lattice(double x, double y, double z)
	{
		a1.x() = x  ; a1.y() = 0.0; a1.z() = 0.0;
		a2.x() = 0.0; a2.y() = y  ; a2.z() = 0.0;
		a3.x() = 0.0; a3.y() = 0.0; a3.z() =   z;
		o = V3ZERO;
		recalculate();
	}

	Eigen::Vector3d delta(const Eigen::Vector3d& pos1, const Eigen::Vector3d& pos2) const
	{
		Eigen::Vector3d diff = pos1 - pos2;
		diff -=  a1 * floor(0.5 + b1.dot(diff))
			+a2 * floor(0.5 + b2.dot(diff))
			+a3 * floor(0.5 + b3.dot(diff));
		return diff;
	}

	Eigen::Vector3d wrap_delta(const Eigen::Vector3d& pos1) const
	{
		Eigen::Vector3d diff = pos1 - o;
		Eigen::Vector3d result(0., 0., 0.);
		result -= a1 * floor(0.5 + b1.dot(diff))
			+ a2 * floor(0.5 + b2.dot(diff))
			+ a3 * floor(0.5 + b3.dot(diff));
		return result;
	}
	double volume() { return a1.dot(a2.cross(a3)); }

	std::string _filename() { return filename; }

	Eigen::Vector3d _a1() { return a1; }
	Eigen::Vector3d _a2() { return a2; }
	Eigen::Vector3d _a3() { return a3; }
	Eigen::Vector3d _o() { return o; }
	Eigen::Vector3d _g1() { return g1; }
	Eigen::Vector3d _g2() { return g2; }
	Eigen::Vector3d _g3() { return g3; }
	
	double _x() { return a1.x(); }
	double _y() { return a2.y(); }
	double _z() { return a3.z(); }

	void _x(const double x) { a1.x() = x; }
	void _y(const double y) { a2.y() = y; }
	void _z(const double z) { a3.z() = z; }
	
	private:
	void recalculate()
	{
		double rvol2pi = 2. * N_PI / volume();
		g1 = a2.cross(a3) * rvol2pi;
		g2 = a3.cross(a1) * rvol2pi;
		g3 = a1.cross(a2) * rvol2pi;

		b1.x() = 1. / a1.x();
		b2.y() = 1. / a2.y();
		b3.z() = 1. / a3.z();
		b1.y()=b1.z()=b2.x()=b2.z()=b3.x()=b3.y()=0.;
	}
};
#endif
