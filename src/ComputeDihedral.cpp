#include <iostream>
#include <iomanip>
#include "Dihedral.hpp"
#include "ComputeDihedral.hpp"
using namespace std;

ComputeDihedral::ComputeDihedral(vector<Dihedral>& dihedralVector)
{
	this -> ptr_dihedralVector = &dihedralVector;

	reset();
}

void ComputeDihedral::reset()
{
	sum_vdihedral = 0.0;
}

double ComputeDihedral::compute_force()
{
	reset();

	for (auto& dihed: *ptr_dihedralVector)
	{
		dihed.calc_force();

		dihed.ptr_atom1->force += dihed.force1;
		dihed.ptr_atom2->force += dihed.force2;
		dihed.ptr_atom3->force += dihed.force3;
		dihed.ptr_atom4->force += dihed.force4;

		sum_vdihedral += dihed.energy;
	}

	return sum_vdihedral;
}

void ComputeDihedral::show_dihedral(const Dihedral& dihed)
{
	cout
	<< setw(10) << fixed << setprecision(3)
	<< dihed.chi * RAD2DEG
	<< setw(5)
	<< dihed.ptr_atom1 -> PDBAtomName
	<< setw(5)
	<< dihed.ptr_atom2 -> PDBAtomName
	<< setw(5)
	<< dihed.ptr_atom3 -> PDBAtomName
	<< setw(5)
	<< dihed.ptr_atom4 -> PDBAtomName
	<< setw(7)
	<< dihed.ptr_atom1 -> PSFIndex
	<< setw(7)
	<< dihed.ptr_atom2 -> PSFIndex
	<< setw(7)
	<< dihed.ptr_atom3 -> PSFIndex
	<< setw(7)
	<< dihed.ptr_atom4 -> PSFIndex
	<< '\n';
}

void ComputeDihedral::show_all_dihedrals()
{
	int icnt = 0;
	for (auto& dihed: *ptr_dihedralVector)
	{
		cout << setw(6) << ++icnt;
		show_dihedral(dihed);
	}
}
