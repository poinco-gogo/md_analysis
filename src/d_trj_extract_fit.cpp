#include <iostream>
#include <iomanip>
#include <sstream>
#include "PSF.hpp"
#include "common.hpp"
#include "PDB.hpp"
#include "ReadDCD.hpp"
#include "Index.hpp"
#include "ComputeRMSD.hpp"
using namespace std;
int main(int argc, char** argv)
{
	if (argc < 6)
	{
		cout <<
		"\n\nusage: ./a.out psf dcd cv_index fitfile fit.ind\n\n";
		return 1;
	}

	output_args(argc, argv);

	PSF PSFFile(argv[1]);
	ReadDCD DCD(&PSFFile.atomVector);

	if (!DCD.open_dcd_read(argv[2]))
		return 1;

	vector<int> cv_index;
	Index IND(argv[3]);
	if (!IND.load_indexes(&cv_index))
		return 1;

	cout << "REMARK " << cv_index.size() << " indexes found.\n";

	PDB PDBFile(argv[4]);
	vector<Atom> refAtomVector;
	refAtomVector.resize(PSFFile.atomVector.size());
	if (!PDBFile.LoadCoords(refAtomVector))
		return 1;

	vector<int> fit_index;
	Index IND2(argv[5]);
	if (!IND2.load_indexes(&fit_index))
		return 1;

	cout << "REMARK " << fit_index.size() << " indexes found.\n";

	ComputeRMSD SOLVER(&fit_index, &refAtomVector);
	SOLVER.remove_ref_com();

	cout << setprecision(8) << scientific;
	while(DCD.read_1step())
	{
		SOLVER.set_tgtAtomVector(&PSFFile.atomVector);

		SOLVER.compute_rmsd();

		SOLVER.add_ref_com(PSFFile.atomVector);

		cout << setw(16) << DCD._nsteps();

		for (auto& id: cv_index)
		{
			Atom& atom = PSFFile.atomVector[id - 1];

			cout 
				<< setw(16) << atom.position.x()
				<< setw(16) << atom.position.y()
				<< setw(16) << atom.position.z();
		}
		cout << '\n';
	}
}
