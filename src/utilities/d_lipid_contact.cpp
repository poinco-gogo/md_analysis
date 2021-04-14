#include <iostream>
#include <fstream>
#include <sstream>
#include "PSF.hpp"
#include "ReadDCD.hpp"
#include "common.hpp"
#include "Lattice.hpp"
using namespace std;

int main(int argc, char** argv)
{
	if (argc != 4)
	{
		cout << 
		"\nD_LIPID_CONTACT\n"
		"\nCalculate lipid contact number\n"
		"\nusage: ./a.out psf dcd cutoff\n\n";
		return 1;
	}

	output_args(argc, argv);

	PSF PSFFile(argv[1]);
	int natom = PSFFile.atomVector.size();

	vector<int> iNH3L, iE219CD, iprthv;

	for (int i = 0; i < natom; i++)
	{
		Atom& atom = PSFFile.atomVector[i];

		string seg = atom.PSFSegmentName;

		if (atom.PSFAtomName == "NH3L")
			iNH3L.push_back(i);
		if (
				atom.PSFResID == 219 &&
				atom.PDBAtomName == "CD" &&
				(seg == "PRB" || seg == "PRD")
				)
			iE219CD.push_back(i);
		if (seg == "PRT"
				|| seg == "PRA"
				|| seg == "PRB"
				|| seg == "PRC"
				|| seg == "PRD"
				|| seg == "PRI"
				|| seg == "PRJ"
				|| seg == "PRK"
				|| seg == "PRL")
		{
			if (atom.mass > 2)
				iprthv.push_back(i);
		}

	}

	cout << "REMARK E219CD atom(s): " << iE219CD.size() << '\n';
	cout << "REMARK   NH3L atom(s): " << iNH3L.size() << '\n';
	cout << "REMARK  prthv atom(s): " << iprthv.size() << '\n';

	double cutoff = atof(argv[3]);

	cout << "REMARK  cutoff value: " << cutoff << '\n';
	cout << "REMARK ==========================================\n";

	ReadDCD DCD(&PSFFile.atomVector);
	if (!DCD.open_dcd_read(argv[2]))
		return 0;

	cout << setprecision(8) << scientific;
	while(DCD.read_1step())
	{
		const double lx = DCD._boxx();
		const double ly = DCD._boxy();
		const double lz = DCD._boxz();

		Lattice lattice(lx, ly, lz);

		Eigen::Vector3d com = PSFFile.getcomof_zerobased(iprthv);

		double sum1 = 0;
		double sum2 = 0;
		for (auto& i: iNH3L)
		{
			Atom& N = PSFFile.atomVector[i];
			Eigen::Vector3d del = N.position - com;
			N.position += lattice.wrap_delta(del);

			for (auto& j: iE219CD)
			{
				Atom& CD = PSFFile.atomVector[j];

				double d = (N.position - CD.position).norm();

				double r = d / cutoff;
				double r2 = r * r;
				double r6 = r2 * r2 * r2;
				double r12 = r6 * r6;
				double r10 = r2 * r2 * r6;
				double r20 = r10 * r10;

				double numer1 = 1.0 - r6;
				double denom1 = 1.0 - r12;

				double numer2 = 1.0 - r10;
				double denom2 = 1.0 - r20;

				sum1 += numer1 / denom1;
				sum2 += numer2 / denom2;
			}
		}

		cout
			<< setw(10) << DCD._nsteps()
			<< setw(16) << sum1
			<< setw(16) << sum2
			<< '\n';
	}
}
