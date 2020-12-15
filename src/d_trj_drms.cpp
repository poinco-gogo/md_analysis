#include <iostream>
#include <cstdlib>
#include <iomanip>
#include <fstream>
#include "PSF.hpp"
#include "PDB.hpp"
#include "Index.hpp"
#include "ReadDCD.hpp"

using namespace std;

int main(int argc, char** argv)
{
	if (argc < 6)
	{
		cout <<
		"\nDrms calculation.\n"
		"\nusage: ./a.out psf dcd ref ind1 ind2\n\n";
		return 0;
	}

	output_args(argc, argv);

	vector<Atom> atomVector;

	PSF PSFFile(argv[1]);
	atomVector = PSFFile.atomVector;

	ReadDCD DCD(&atomVector);
	if (!DCD.open_dcd_read(argv[2]))
		return 1;

	vector<Atom> refAtomVector;
	refAtomVector.resize(atomVector.size());

	PDB PDBFile(argv[3]);
	if (!PDBFile.LoadCoords(refAtomVector))
		return 1;

	vector<int> sel1, sel2;

	Index IDX1(argv[4]);
	Index IDX2(argv[5]);
	if (!IDX1.load_indexes(&sel1) || !IDX2.load_indexes(&sel2))
		return 0;

	cout << "REMARK Selection 1 contains following "
		<< sel1.size() << " atoms:\n";
	for (auto& i: sel1) PSFFile.showAtom(atomVector[i - 1]);
	cout << "REMARK Selection 2 contains following "
		<< sel2.size() << " atoms:\n";
	for (auto& i: sel2) PSFFile.showAtom(atomVector[i - 1]);

	const unsigned int npairs = sel1.size() * sel2.size();

	Eigen::MatrixXd dref(sel1.size(), sel2.size());

	for (int i = 0; i < sel1.size(); i++)
	{
		Atom& iatom = refAtomVector[sel1[i] - 1];

		for (int j = 0; j < sel2.size(); j++)
		{
			Atom& jatom = refAtomVector[sel2[j] - 1];

			dref(i, j) = (iatom.position - jatom.position).norm();
		}
	}

	cout << scientific << setprecision(8);
	while (DCD.read_1step())
	{
		double drms = 0.;

		for (int i = 0; i < sel1.size(); i++)
		{
			Atom& iatom = atomVector[sel1[i] - 1];

			for (int j = 0; j < sel2.size(); j++)
			{
				Atom& jatom = atomVector[sel2[j] - 1];

				double dij
				= (iatom.position - jatom.position).norm();

				double del = dref(i, j) - dij;
				drms += del * del;
			}
		}

		drms = sqrt(drms/npairs);

		cout 
		<< setw(10) << DCD._nsteps()
		<< setw(16) << drms << '\n';
	}
}
