#ifndef ___CLASS_COMPUTE_DIHEDRAL
#define ___CLASS_COMPUTE_DIHEDRAL

#include <vector>
#include "Dihedral.hpp"

class ComputeDihedral
{
	private:

	double sum_vdihedral;
	std::vector<Atom>* ptr_atomVector;
	std::vector<Dihedral>* ptr_dihedralVector;
	
	public:

	ComputeDihedral(){};
	ComputeDihedral(std::vector<Dihedral>& dihedralVector);

	void reset();

	double compute_force();

	void show_dihedral(const Dihedral& dihed);
	void show_all_dihedrals();
};
#endif
