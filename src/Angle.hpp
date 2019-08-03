#ifndef ___CLASS_ANGLE
#define ___CLASS_ANGLE

#include <cmath>
#include <Eigen/Core>
#include "common.hpp"
#include "Atom.hpp"

class Angle
{
	// V(angle) = Ktheta(theta - theta0)**2
	// V(Urey-Bradley) = Kub(S - S0)**2
	// Ktheta: kcal/mole/rad**2
	// theta0: degrees
	// Kub: kcal/mol/A**2 (Urey-Bradley)
	// S0: A
		
	public:

	Atom* ptr_atom1;
	Atom* ptr_atom2;
	Atom* ptr_atom3;
	
	std::string at1;
	std::string at2;
	std::string at3;

	double Kt;
	double t0;
	double t;
	double Kub;
	double s0;
	double s;

	Eigen::Vector3d r12;
	Eigen::Vector3d r32;
	Eigen::Vector3d r13;

	double invabsr12;
	double invabsr32;

	double r12r32;
	
	double energy;
	double energy_ub;

	double diff;
	
	Eigen::Vector3d force1;
	Eigen::Vector3d force2;
	Eigen::Vector3d force3;
	Eigen::Vector3d force_ub;


	//constructor
	Angle(std::string at1, std::string at2, std::string at3, double Kt, double t0, double Kub, double s0)
	{
		this -> at1 = at1;
		this -> at2 = at2;
		this -> at3 = at3;
		this -> Kt = Kt;
		this -> t0 = t0;
		this -> Kub = Kub;
		this -> s0 = s0;
	}
	Angle(Atom* ptr_atom1, Atom* ptr_atom2, Atom* ptr_atom3, double Kt, double t0, double Kub, double s0)
	{
		this -> ptr_atom1 = ptr_atom1;
		this -> ptr_atom2 = ptr_atom2;
		this -> ptr_atom3 = ptr_atom3;
		this -> Kt = Kt;
		this -> t0 = t0 * DEG2RAD;
		this -> Kub = Kub;
		this -> s0 = s0;

		// initialization
		t = 0.0;
		s = 0.0;
		energy = 0.0;
		energy_ub = 0.0;
		diff = 0.0;
		r12.x() = 0.0;
		r12.y() = 0.0;
		r12.z() = 0.0;
		r32.x() = 0.0;
		r32.y() = 0.0;
		r32.z() = 0.0;
		r13.x() = 0.0;
		r13.y() = 0.0;
		r13.z() = 0.0;
		force1.x() = 0.0;
		force1.y() = 0.0;
		force1.z() = 0.0;
		force2.x() = 0.0;
		force2.y() = 0.0;
		force2.z() = 0.0;
		force3.x() = 0.0;
		force3.y() = 0.0;
		force3.z() = 0.0;
		force_ub.x() = 0.0;
		force_ub.y() = 0.0;
		force_ub.z() = 0.0;
	}

	// member
	inline void calc_force()
	{
		// based on NAMD code !

		r12 = ptr_atom1 -> position - ptr_atom2 -> position;
		r32 = ptr_atom3 -> position - ptr_atom2 -> position;
	
		invabsr12 = 1.0 / r12.norm();
		invabsr32 = 1.0 / r32.norm();

		double cos_theta = r12.dot(r32) * invabsr12 * invabsr32;

		// 僅かな誤差を（あれば）除く
		if (cos_theta > 1.0) cos_theta = 1.0;
		else if (cos_theta < -1.0) cos_theta = -1.0;

		r12r32 = r12.dot(r32);

		t = acos(cos_theta);
		
		diff = t - t0;

		energy = Kt * diff * diff;

		// constant factor 2K(theta-theta0)/sin(theta)の計算
		double sin_theta = sqrt(1.0 - cos_theta * cos_theta);
		if (sin_theta < 1.e-6)
		{
			// NAMD によると、bonds がほとんどparallelになっている
			// 場合を考慮する
			if (diff < 0.)
				diff = 2.0 * Kt;
			else
				diff = -2.0 * Kt;
		}
		else
		{
			diff *= -2.0 * Kt / sin_theta;
		}

		double c1 = diff * invabsr12;
	        double c2 = diff * invabsr32;

		force1 = c1 * (invabsr12 * cos_theta * r12 - invabsr32 * r32);
		force3 = c2 * (invabsr32 * cos_theta * r32 - invabsr12 * r12);
		force2 = -force1 - force3;
	}
	inline void calc_force_ub()
	{
		r13 = ptr_atom1 -> position - ptr_atom3 -> position;
	
		s = r13.norm();
		energy_ub = Kub * (s - s0) * (s - s0);
		force_ub  = -2.0 * Kub * (s - s0) / s * r13;
	}
};
#endif
