#ifndef ___CLASS_MOLECULE
#define ___CLASS_MOLECULE

#include <vector>
#include <string>
#include "Atom.hpp"
#include "common.hpp"

class Molecule
{
	private:

	bool is_water;
	bool is_protein;

	public:
	
	std::string mol_name;

	int mol_index;

	Eigen::Vector3d   vcom;

	std::vector<Atom*> ptr_atoms;

	Molecule(int mol_index, std::string mol_name);

	void show_mol_info();
	
	void calc_com();

	bool _is_water() { return is_water; }
};
#endif
