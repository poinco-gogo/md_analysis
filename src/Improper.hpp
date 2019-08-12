// based on implementation of NAMD
#ifndef ___CLASS_IMPROPER
#define ___CLASS_IMPROPER

#include <Eigen/Core>
#include <Eigen/Geometry>
#include "common.hpp"
#include "Atom.hpp"
#include <cmath>

class Improper
{
	// V(improper) = Kpsi(psi-psi0)**2
	// Kpsi: kcal/mole/rad**2
	// psi0: degrees
	private:
	double chi;
	double inA, inB, inC;
	double cos_psi, sin_psi;
	Eigen::Vector3d hA, hB, hC;

	public:
	bool is_wildcard;
	std::string at1, at2, at3, at4;
	
	Atom* ptr_atom1;
	Atom* ptr_atom2;
	Atom* ptr_atom3;
	Atom* ptr_atom4;

	double energy;
	double k;
	double Kpsi;
	double psi0;
	double psi;

	Eigen::Vector3d r12, r23, r34;

	Eigen::Vector3d force1, force2, force3, force4;

	// constructor
	Improper(std::string at1, std::string at2, std::string at3, std::string at4, double Kpsi, double psi0)
	{
		this -> at1 = at1;
		this -> at2 = at2;
		this -> at3 = at3;
		this -> at4 = at4;
		this -> Kpsi = Kpsi;
		this -> psi0 = psi0;
	}
	Improper(Atom* ptr_atom1, Atom* ptr_atom2, Atom* ptr_atom3, Atom* ptr_atom4, double Kpsi, double psi0)
	{
		this -> ptr_atom1 = ptr_atom1;
		this -> ptr_atom2 = ptr_atom2;
		this -> ptr_atom3 = ptr_atom3;
		this -> ptr_atom4 = ptr_atom4;
		this -> Kpsi = Kpsi;
		this -> psi0 = psi0 * DEG2RAD;
	}

	inline void calc_force()
	{
		calc_chi();

		double diff = chi - psi0;

		// minimum image convention
		if (diff < -N_PI)     diff += N_TWOPI;
		else if (diff > N_PI) diff -= N_TWOPI;

		energy = Kpsi * diff * diff;

		double k = 2. * Kpsi * diff;

		calc_derivatives();
	}

	inline void calc_chi()
	{
		r12 = ptr_atom1->position - ptr_atom2->position;
		r23 = ptr_atom2->position - ptr_atom3->position;
		r34 = ptr_atom3->position - ptr_atom4->position;

		Eigen::Vector3d A = r12.cross(r23);
		Eigen::Vector3d B = r23.cross(r34);
		Eigen::Vector3d C = r23.cross(  A);

		double inA = 1. / A.norm();
		double inB = 1. / B.norm();
		double inC = 1. / C.norm();

		hA = A * inA;
		hB = B * inB;
		hC = C * inC;

		cos_psi = hA.dot(hB);
		sin_psi = hC.dot(hB);

		chi = -atan2(sin_psi, cos_psi);
	}

	inline void calc_derivatives()
	{
		if (std::abs(sin_psi) > 0.1)
		{
			k /= -sin_psi;

			Eigen::Vector3d r13 = r12 + r23;

			Eigen::Vector3d v1 = inA * (cos_psi * hA - hB);
			Eigen::Vector3d v2 = inB * (cos_psi * hB - hA);

			force1 = k * r23.cross(v1);
			force2 = k * (r34.cross(v2) - r13.cross(v1));
			force4 = k * r23.cross(v2);
			force3 = -force1 - force2 - force4;
		}
		else
		{
			k /= cos_psi;

			Eigen::Vector3d r24 = r23 + r34;

			Eigen::Vector3d v3 = inB * (sin_psi * hB - hC);
			Eigen::Vector3d v4 = inC * (sin_psi * hC - hB);

			force1.x() 
			= k * ((r23.y()*r23.y()+r23.z()*r23.z())*v4.x()
				       - r23.x()*r23.y()*v4.y()
				       - r23.x()*r23.z()*v4.z());
			force1.y() = k * (-r23.y()*r23.x()*v4.x()
			+(r23.z()*r23.z() + r23.x()*r23.x())*v4.y()
				       - r23.y()*r23.z()*v4.z());
			force1.z() = k * (-r23.z()*r23.x()*v4.x()
				       - r23.z()*r23.y()*v4.y()
			+(r23.x()*r23.x() + r23.y()*r23.y())*v4.z());
			
			force3.x()
			= k * ((r12.y()*r23.y()+r12.z()*r23.z())*v4.x()
			+(-2.*r23.x()*r12.y()+r12.x()*r23.y())*v4.y()
			+(-2.*r12.z()*r23.x()+r12.x()*r23.z())*v4.z()
			+ r24.z()*v3.y() - r24.y()*v3.z());
			force3.y()
			= k * ((-2.*r12.x()*r23.y()+r12.y()*r23.x())*v4.x()
			+(r12.z()*r23.z()+r12.x()*r23.x())*v4.y()
			+(r12.y()*r23.z()-2.*r12.z()*r23.y())*v4.z()
			- r24.z()*v3.x()+r24.x()*v3.z());
			force3.z()
			= k * ((-2.*r12.x()*r23.z()+r12.z()*r23.x())*v4.x()
			-(2.*r12.y()*r23.z()-r12.z()*r23.y())*v4.y()
			+(r12.x()*r23.x()+r12.y()*r23.y())*v4.z()
			+ r24.y()*v3.x()-r24.x()*v3.y());
			
			force4 = k * r23.cross(v3);

			force2 = -force1 - force3 - force4;
		}
	}
};
#endif
