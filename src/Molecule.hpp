#ifndef ___CLASS_MOLECULE
#define ___CLASS_MOLECULE

#include <vector>
#include "Atom.hpp"

class Molecule
{
	public:

	int mol_index;

	Eigen::Vector3d   vcom;

	std::vector<Atom*> ptr_atoms;

	void show_mol_info();

	Molecule(int mol_index);
};
#endif
