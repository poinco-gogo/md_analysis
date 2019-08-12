#ifndef ___CLASS_DIHEDRAL
#define ___CLASS_DIHEDRAL

#include <cmath>
#include <iomanip>
#include <Eigen/Core>
#include <Eigen/Geometry>
#include "Atom.hpp"
#include "common.hpp"

class Dihedral
{
	// V(dihedral) = Kchi(1 + cos(n(chi) - delta))
	// Kchi: kcal/mole
	// n: periodicity
	// delta: degrees
	private:
	double inA, inB, inC;
	double sin_phi, cos_phi;
	Eigen::Vector3d hA, hB, hC;
		
	public:
	bool is_wildcard;
	Atom* ptr_atom1;
	Atom* ptr_atom2;
	Atom* ptr_atom3;
	Atom* ptr_atom4;
	std::string at1, at2, at3, at4;
	Eigen::Vector3d r12, r23, r34;
	double k;
	double Kchi;
	double n;
	double delta;
	double chi;
	double energy;
	Eigen::Vector3d force1, force2, force3, force4;

	//constructor
	Dihedral(){}; // do not remove this...cmap uses it
	Dihedral(std::string at1, std::string at2, std::string at3, std::string at4, double Kchi, double n, double delta)
	{
		this -> at1 = at1;
		this -> at2 = at2;
		this -> at3 = at3;
		this -> at4 = at4;
		this -> Kchi = Kchi;
		this -> n = n;
		this -> delta = delta;
	}
	Dihedral(Atom* ptr_atom1, Atom* ptr_atom2, Atom* ptr_atom3, Atom* ptr_atom4, double Kchi, double n, double delta)
	{
		this -> ptr_atom1 = ptr_atom1;
		this -> ptr_atom2 = ptr_atom2;
		this -> ptr_atom3 = ptr_atom3;
		this -> ptr_atom4 = ptr_atom4;
		this -> Kchi = Kchi;
		this -> n = n;
		this -> delta = delta * DEG2RAD;
		energy = 0.;
	}

	// ref. NAMD ComputeDihedrals.C
	inline void calc_force()
	{
		calc_chi();

		energy = Kchi * (1. + cos(n * chi - delta));
		
		k = n * Kchi * sin(n * chi - delta);

		calc_derivatives();

	}
	inline void calc_derivatives_for_cmap()
	{
		// do not remove this subroutine ! referenced from Cmap.hpp
		calc_chi();

		k = 1.0;

		calc_derivatives();
	}

	inline void calc_chi()
	{
		r12 = ptr_atom1 -> position - ptr_atom2 -> position;
		r23 = ptr_atom2 -> position - ptr_atom3 -> position;
		r34 = ptr_atom3 -> position - ptr_atom4 -> position;
		
		Eigen::Vector3d A = r12.cross(r23);
		Eigen::Vector3d B = r23.cross(r34);
		Eigen::Vector3d C = r23.cross(  A);

		inA = 1. / A.norm();
		inB = 1. / B.norm();
		inC = 1. / C.norm();

		hA = A * inA;
		hB = B * inB;
		hC = C * inC;

		cos_phi = hA.dot(hB);
		sin_phi = hC.dot(hB);
		chi = -atan2(sin_phi, cos_phi);
	}

	inline void calc_derivatives()
	{
		if (std::abs(sin_phi) > 0.1)
		{
			k /= -sin_phi;

			Eigen::Vector3d r13 = r12 + r23;

			Eigen::Vector3d v1 = inA * (cos_phi * hA - hB);
			Eigen::Vector3d v2 = inB * (cos_phi * hB - hA);

			force1 = k * r23.cross(v1);
			force2 = k * (r34.cross(v2) - r13.cross(v1));
			force4 = k * r23.cross(v2);
			force3 = -force1 - force2 - force4;
		}
		else
		{
			k /= cos_phi;

			Eigen::Vector3d r24 = r23 + r34;

			Eigen::Vector3d v3 = inB * (sin_phi * hB - hC);
			Eigen::Vector3d v4 = inC * (sin_phi * hC - hB);

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

	inline void show_profile()
	{
		std::cout 
		<< std::setw(5) << ptr_atom1 -> PDBAtomName
		<< std::setw(5) << ptr_atom2 -> PDBAtomName
		<< std::setw(5) << ptr_atom3 -> PDBAtomName
		<< std::setw(5) << ptr_atom4 -> PDBAtomName
		<< std::setw(7) << ptr_atom1 -> PSFIndex
		<< std::setw(5) << ptr_atom2 -> PSFIndex
		<< std::setw(5) << ptr_atom3 -> PSFIndex
		<< std::setw(5) << ptr_atom4 -> PSFIndex
		<< '\n';
	}
	bool is_phi_of(int tgt_resID)
	{
		// if phi    C    N   CA    C
		// if ALAD   CLP  NL  CA    CRP
		if (ptr_atom3->PSFResID != tgt_resID)
			return false;
		return ((
		   ptr_atom1->PDBAtomName=="C" 
		&& ptr_atom2->PDBAtomName=="N"
		&& ptr_atom3->PDBAtomName=="CA" 
		&& ptr_atom4->PDBAtomName=="C") ||
		(  
		   ptr_atom1->PDBAtomName=="CLP"
		&& ptr_atom2->PDBAtomName=="NL"
		&& ptr_atom3->PDBAtomName=="CA"
		&& ptr_atom4->PDBAtomName=="CRP"
		))	? true : false;

	}
	bool is_psi_of(int tgt_resID)
	{
		// if psi    N   CA    C    N
		// if ALAD   NL  CA    CRP  NR
		if (ptr_atom2->PSFResID != tgt_resID)
			return false;
		return ((
		   ptr_atom1->PDBAtomName=="N" 
		&& ptr_atom2->PDBAtomName=="CA"
		&& ptr_atom3->PDBAtomName=="C" 
		&& ptr_atom4->PDBAtomName=="N") ||
		(		
		   ptr_atom1->PDBAtomName=="NL" 
		&& ptr_atom2->PDBAtomName=="CA"
		&& ptr_atom3->PDBAtomName=="CRP" 
		&& ptr_atom4->PDBAtomName=="NR"
		)) ? true : false;
	}

	Eigen::Vector3d metric(int PSFIndex)
	{
		Eigen::Vector3d vtmp(0., 0., 0.);

		if (ptr_atom1->PSFIndex == PSFIndex)
			return force1;
		else if (ptr_atom2->PSFIndex == PSFIndex)
			return force2;
		else if (ptr_atom3->PSFIndex == PSFIndex)
			return force3;
		else if (ptr_atom4->PSFIndex == PSFIndex)
			return force4;
		else
			return vtmp;
	}
};
#endif
