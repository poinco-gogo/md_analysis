#include <iostream>
#include <iomanip>
#include <cstdlib>
#include "common.hpp"
#include "PDB.hpp"
#include "PSF.hpp"
#include "Dihedral.hpp"
#include "ComputeDihedral.hpp"
#include "LoadParm.hpp"
using namespace std;

int main(int argc, char** argv)
{
	if (argc < 4)
	{
		cout <<
		"\nD_ORACLE6\n"
		"\nOUTPUT ALL DIHEDRALS IN SYSTEM:\n"
		"\nusage: ./a.out psf pdb parm_list [ResID]\n\n";
		return 1;
	}

	cout << "REMARK ";
	for (int i = 0; i < argc; i++)
		cout << argv[i] << ' ';
	cout << "\nREMARK D_ORACLE6\n";

	vector<Atom> atomVector;
	PSF PSFFile(argv[1], &atomVector);
	PDB PDBFile(argv[2]);
	if (!PDBFile.LoadCoords(atomVector))
		return 1;

	LoadParm ALL22;
	if (!ALL22.open_fi(argv[3])) return 1;
	if (!PSFFile.set_dihedral_parm(ALL22.dihedralParmVector)) return 1;

	cout << "REMARK =============================\n";
	
	ComputeDihedral JOB(PSFFile.dihedralVector);
	
	vector<int> tgt_dihedNo(2, -1);
	if (argv[4])
	{
		int ResID = atoi(argv[4]);
		for (int i = 0; i < PSFFile.dihedralVector.size(); i++)
                {
			Dihedral& di = PSFFile.dihedralVector[i];
			if (di.is_phi_of(ResID))
				tgt_dihedNo[0] = i + 1;
			if (di.is_psi_of(ResID))
				tgt_dihedNo[1] = i + 1;
		}
		if (tgt_dihedNo[0] < 0 || tgt_dihedNo[1] < 0)
		{
			cerr << "error: something is wrong.\n";
			return 0;
		}

		Dihedral& di1 = PSFFile.dihedralVector[tgt_dihedNo[0] - 1];
		Dihedral& di2 = PSFFile.dihedralVector[tgt_dihedNo[1] - 1];
		
		Atom* ptr_atom3 = di1.ptr_atom3;

		cout 
		<< "REMARK OUTPUT DIHEDRAL ANGLE OF " 
		<< ptr_atom3 -> PSFResName << ptr_atom3 -> PSFResID << '\n';
		
		cout << "REMARK DIHED No. " << tgt_dihedNo[0] << ':';
		di1.calc_chi();
		JOB.show_dihedral(di1);

		cout << "REMARK DIHED No. " << tgt_dihedNo[1] << ':';
		di2.calc_chi();
		JOB.show_dihedral(di2);

		return 1;
	}

	// actually not needed, but...
	JOB.reset();

	JOB.show_all_dihedrals();
}
