#include <iostream>
#include <sstream>
#include "common.hpp"
#include "Dihedral.hpp"
#include "LoadParm.hpp"
#include "PSF.hpp"
#include "ReadDCD.hpp"
using namespace std;
int main(int argc, char** argv)
{
	if (argc < 5)
	{
		cout << "\nD_TRJ_DIHEDRAL\n"
		"\nDIHEDRAL ANGLE CALCULATION ALONG DCD\n"
		"\nusage: ./a.out psf dcd parm dihedralNo. [dihedralNo.]\n"
		"\nnote: you may specify residue No. by -r 43, for example."
		"\nif -r option is used, dihedralNo. is ignored.\n\n";
		return 1;
	}

	cout << "REMARK ";
	for (int i = 0; i < argc; i++)
		cout << argv[i] << ' ';
	cout << "\nREMARK D_TRJ_DIHEDRAL\n";
	
	bool is_r = false;

	int tgt_resID = -1;
	vector<int> tgt_dihedNo(2, -1);

	string s(argv[4]);
	istringstream is(s);
	is >> tgt_dihedNo[0];

	if (!is)
	{
		if (is.str() == "-r")
		{
			cout << "REMARK -r option detected.\n";
			is_r = true;
		}
		else
		{	
			cerr << "error: improper argment : " << argv[4] << '\n';
			return 0;
		}
	}
	if (argv[5])
	{
		string s(argv[5]);
		istringstream is(s);
		if (is_r)
			is >> tgt_resID;
		else
			is >> tgt_dihedNo[1];

		if (!is)
		{
			cerr << "error: improper argument:" << argv[5] << '\n';
			return 0;
		}
	}
	else if (is_r)
	{
		cerr << "error: you must specify resID.\n";
		return 0;
	}
	
	PSF PSFFile(argv[1]);
	ReadDCD DCD(&PSFFile.atomVector);
	LoadParm ALL22(argv[3]);
	if (!PSFFile.set_dihedral_parm(ALL22.dihedralParmVector)) return 1;
	
	if (is_r)
	{
		for (int i = 0; i < PSFFile.dihedralVector.size(); i++)
		{
			Dihedral& di = PSFFile.dihedralVector[i];

			if (di.is_phi_of(tgt_resID))
				tgt_dihedNo[0] = i + 1;
			if (di.is_psi_of(tgt_resID))
				tgt_dihedNo[1] = i + 1;
		}
	}

	if (tgt_dihedNo[0] < 0 || (argv[5] && tgt_dihedNo[1] < 0))
	{
		cerr << "error: something is wrong.\n";
		return 0;
	}
	
	if (is_r)
	{
		string s =
		PSFFile.dihedralVector[tgt_dihedNo[0]-1].ptr_atom3->PSFResName;

		cout << 
		"REMARK Output dihedral angle OF " << s << tgt_resID << '\n';
	}

	cout << "REMARK Output " << tgt_dihedNo[0] << " th dihedral angle\n";
	
	cout << "REMARK Output " << tgt_dihedNo[1] << " th dihedral angle\n";

	
	cout << "REMARK Dihedral profile :\n"
	<< "REMARK Dihedral No." << tgt_dihedNo[0];
	PSFFile.dihedralVector[tgt_dihedNo[0] - 1].show_profile();
	if (tgt_dihedNo[1] > 0)
	{
		cout << "REMARK Dihedral No." << tgt_dihedNo[1];
		PSFFile.dihedralVector[tgt_dihedNo[1] - 1].show_profile();
	}

	if (!DCD.open_dcd_read(argv[2]))
		return 1;

	cout << fixed << setprecision(3);
	while(DCD.read_1step())
	{
		PSFFile.dihedralVector[tgt_dihedNo[0] - 1].calc_chi();

		cout << setw(10) << DCD._nsteps()
		<< setw(10)
		<< PSFFile.dihedralVector[tgt_dihedNo[0] - 1].chi * RAD2DEG;

		if (tgt_dihedNo[1] > 0)
		{
			PSFFile.dihedralVector[tgt_dihedNo[1] - 1].calc_chi();
			cout
			<< setw(10)
			<< PSFFile.dihedralVector[tgt_dihedNo[1] - 1].chi * RAD2DEG;
		}

		cout << '\n';
	}
}
