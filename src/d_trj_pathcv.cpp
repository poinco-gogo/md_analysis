#include <iostream>
#include <cstdlib>
#include <iomanip>
#include <fstream>
#include "PSF.hpp"
#include "PDB.hpp"
#include "Index.hpp"
#include "ComputePathCV.hpp"
#include "ComputeRMSD.hpp"
#include "ReadDCD.hpp"

using namespace std;

int main(int argc, char** argv)
{
	if (argc != 6)
	{
		cout <<
		"\nPathCV analysis\n"
		"\nusage: ./a.out psf dcd fitfile ind path.dat\n\n";
		return 0;
	}

	output_args(argc, argv);

	vector<Atom> atomVector;

	PSF PSFFile(argv[1], &atomVector);

	ReadDCD DCD(&atomVector);
	if (!DCD.open_dcd_read(argv[2]))
		return 1;

	vector<Atom> refAtomVector;
	refAtomVector.resize(atomVector.size());

	PDB PDBFile(argv[3]);
	if (!PDBFile.LoadCoords(refAtomVector))
		return 1;

	vector<int> alignIndex;
	Index IDX(argv[4]);
	if (!IDX.load_indexes(&alignIndex))
		return 0;

	cout << "REMARK CV atoms: ";
	IDX.show_atoms(&alignIndex, &atomVector);

	ComputePathCV JOB(argv[5], &alignIndex, &atomVector);

	ComputeRMSD SOLVER(&alignIndex, &refAtomVector);

	SOLVER.remove_ref_com();

	cout << scientific << setprecision(12);
	while (DCD.read_1step())
	{
		SOLVER.set_tgtAtomVector(&atomVector);

		SOLVER.compute_rmsd();

		SOLVER.add_ref_com(atomVector);

		JOB.calc_pathcv();

		cout 
		<< setw(12) << DCD._nsteps()
		<< setw(20) << JOB._progress()
		<< setw(20) << JOB._distance() << '\n';
	}
}
