#include <iostream>
#include <iomanip>
#include "Molecule.hpp"
using namespace std;

Molecule::Molecule(int mol_index)
{
	this -> mol_index = mol_index;
}

void  Molecule::show_mol_info()
{
	for (auto& ptr_atom: ptr_atoms)
	{
		Atom& atom = *ptr_atom;

		cout
		<< atom.PSFIndex << ' ' << atom.PDBAtomName
		<< " @" << atom.PSFResName << atom.PSFResID
		<< "  " << atom.charge << 'e'
		<< "  " << atom.mass   << "amu\n";
	}
}
