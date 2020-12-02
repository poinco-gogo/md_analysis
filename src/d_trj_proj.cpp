#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <ctime>
#include <cstdlib>
#include <Eigen/Eigenvalues>
#include "MatrixUtil.hpp"
#include "PSF.hpp"
#include "PDB.hpp"
#include "NAMDBin.hpp"
#include "ReadDCD.hpp"
#include "Index.hpp"
#include "common.hpp"
using namespace std;
int main(int argc, char** argv)
{
	if (argc != 7)
	{
		cout <<
		"\nPCA projection\n"
		"\nusage: ./a.out psf[natom] dcd ave[pdb/coor] ind evecs pcNo\n\n";
		return 1;
	}

	output_args(argc, argv);

	vector<Atom> atomVector;

	int natom = 0;
	string spsf(argv[1]);
        string spsfex = spsf.substr(spsf.find_last_of(".") + 1);
        if (spsfex == "psf")
        {
                PSF PSFFile(argv[1], &atomVector);
                natom = atomVector.size();
        }
        else
        {
                natom = atoi(argv[1]);
        }

	vector<int> alignIndex;
	Index IDX(argv[4]);
	if (!IDX.load_indexes(&alignIndex))
		return 1;
	cout << "REMARK " << alignIndex.size() << " indexes loaded.\n";

	ReadDCD TRJFile(&atomVector);
	if (!TRJFile.open_dcd_read(argv[2]))
		return 0;

	vector<Atom> refAtomVector(natom);
	string sref(argv[3]);
	string srefex = sref.substr(sref.find_last_of(".") + 1);
	if (srefex == "pdb")
	{
		PDB PDBFile(argv[3]);
		if (!PDBFile.LoadCoords(refAtomVector))
			return 0;
	}
	else
	{
		NAMDBin JOB(argv[3], "coor");
		JOB.read_fi(refAtomVector);
	}

	const int dim = alignIndex.size() * 3;
	Eigen::MatrixXd m(dim, dim);
	ifstream fi(argv[5]);
	if (!fi)
	{
		cerr << "\nerror: Could not open file " << argv[4] << '\n';
		return 1;
	}
	for (int i = 0; i < dim; i++)
	{
		for (int j = 0; j < dim; j++)
		{
			double dtmp;
			fi >> dtmp;
			m(i, j) = dtmp;
		}
	}

	const int pcNo = atoi(argv[6]);
	cout << "REMARK project onto pc No. " << pcNo << '\n';
	cout << "REMARK ==============\n";

	while (TRJFile.read_1step())
	{
		double sum = 0;

		for (int i = 0; i < alignIndex.size(); i++)
		{
			int id = alignIndex[i] - 1;

			Eigen::Vector3d del
			= atomVector[id].position - refAtomVector[id].position;

			sum +=
			  m(3 * i    , dim - pcNo) * del.x()
			+ m(3 * i + 1, dim - pcNo) * del.y()
			+ m(3 * i + 2, dim - pcNo) * del.z();
		}

		cout
			<< setw(12) << TRJFile._nsteps()
			<< setw(16) << sum
			<< '\n';
	}
}
