#include <iostream>
#include <iomanip>
#include "Molecule.hpp"
using namespace std;

Molecule::Molecule(int mol_index, string mol_name)
{
	this -> mol_index = mol_index;

	this -> mol_name  = mol_name;
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

void Molecule::calc_com()
{
	this -> vcom = V3ZERO;

	double tot_mass = 0;
	for (auto& ptr_atom: ptr_atoms)
	{
		Atom& atom = *ptr_atom;

		this->vcom += atom.mass * atom.position;

		tot_mass += atom.mass;
	}

	this->vcom /= tot_mass;
}
