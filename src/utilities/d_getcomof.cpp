#include <Eigen/Core>
#include <iostream>
#include <cstdlib>
#include <iomanip>
#include "PSF.hpp"
#include "PDB.hpp"
#include "Index.hpp"
#include "common.hpp"
using namespace std;
int main(int argc, char** argv)
{
	if (argc < 4)
	{
		cout << 
		"\nGET COM COORDINATE.\n"
		"\nusage: ./a.out psf pdb ind\n\n";
		return 1;
	}

	output_args(argc, argv);

	vector<Atom> atomVector;
	PSF PSFFile(argv[1], &atomVector);
	PDB PDBFile(argv[2]);
	if (!PDBFile.LoadCoords(atomVector))
		return 1;

	vector<int> ids;
	Index IDX(argv[3]);
	if (!IDX.load_indexes(&ids))
		return 1;

	cout << "REMARK============\n";

	Eigen::Vector3d com = PSFFile.getcomof(ids);

	cout << setprecision(6) << fixed;

	cout << com.transpose() << '\n';
}
