#ifndef ___CLASS_BOND
#define ___CLASS_BOND

#include <Eigen/Core>
#include <Eigen/Geometry>
#include "Atom.hpp"

class Bond
{
	// V(bond) = Kb(b - b0)**2
	// Kb: kcal/mole/A**2
	// b0: A
	
	public:

	//constructor
	Bond(std::string at1, std::string at2, double Kb, double b0)
	{
		this -> at1 = at1;
		this -> at2 = at2;
		this -> Kb = Kb;
		this -> b0 = b0;
	}
	Bond(Atom* ptr_atom1, Atom* ptr_atom2, double Kb, double b0)
	{
		this -> ptr_atom1 = ptr_atom1;
		this -> ptr_atom2 = ptr_atom2;
		this -> Kb    = Kb;
		this -> b0    = b0;
		
		// initialization
		b       = 0.0;
		energy  = 0.0;
		r12.x()   = 0.0;
		r12.y()   = 0.0;
		r12.z()   = 0.0;
		force.x() = 0.0;
		force.y() = 0.0;
		force.z() = 0.0;
		is_XH = false;
	}

	
	// member
	
	bool is_XH;
	
	std::string at1;
	std::string at2;
	
	Atom*   ptr_atom1;
	Atom*   ptr_atom2;
	
	double  Kb;
	double  b0;
	double  b;
	double  energy;

	Eigen::Vector3d r12;
	Eigen::Vector3d force;
/*
	inline void calc_energy()
	{
		r12 = ptr_atom1 -> position - ptr_atom2 -> position;
		b = vabs(r12);
		vbond = Kb * (b - b0) * (b - b0);
	}
*/	
	inline void calc_force()
	{
		r12 = ptr_atom1 -> position - ptr_atom2 -> position;
		b = r12.norm();
		energy = Kb * (b - b0) * (b - b0);
		// これは粒子 1 に作用する force である。
		// 粒子 2 に作用する force = -vfrc
		force = -2.0 * Kb * (b - b0) / b * r12;
	}
	inline void set_flag_XH()
	{
		int mass1 = ptr_atom1 -> mass + 0.5;
		int mass2 = ptr_atom2 -> mass + 0.5;
		if (mass1 == 1 || mass2 == 1)
			is_XH = true;
		else
			is_XH = false;
	}
	inline bool is_same_as(Bond& bond)
	{
		return (at1 == bond.at1 && at2 == bond.at2 && Kb == bond.Kb
				&& b0 == bond.b0);
	}
	inline bool contain(int index)
	{
		return (ptr_atom1->PSFIndex == index || ptr_atom2->PSFIndex == index);
	}
};
#endif
