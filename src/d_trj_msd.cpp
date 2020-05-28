#include <iostream>
#include <iomanip>
#include <cstdlib>
#include "PSF.hpp"
#include "ReadDCD.hpp"
using namespace std;
int main(int argc, char** argv)
{
	if (argc < 6)
	{
		cout << 
		"\nGET MEAN SQUARED DISPLACEMENT FROM DCD\n"
		"\nusage: ./a.out psf dcd resname delta(in MD steps) dimension('xy' or 'xyz')\n\n";
		return 1;
	}
	cout << "REMARK ";
	for (int i = 0; i < argc; i++)
		cout << argv[i] << " ";
	cout << "\nREMARK GET MEAN SQUARED DISPLACEMENT FROM DCD.\n";
	vector<Atom> atomVector;
	PSF PSFFile(argv[1], &atomVector);
	ReadDCD DCD(&atomVector);
	if (!DCD.open_dcd_read(argv[2]))
		return 0;

	PSFFile.define_molecules();

	string stgt(argv[3]);
	vector<int> idx;
	for (auto& mol: PSFFile.moleculeVector)
	{
		if (mol.mol_name == stgt)
			idx.push_back(mol.mol_index);
	}
	if (!idx.size())
	{
		cerr << "\nerror: no molecule named \"" << stgt
			<< "\" found.\n\n";
		return 0;
	}
	cout << "REMARK Calculate MSD for " << stgt << " molecule(s)\n";
	cout << "REMARK " << idx.size() << " target molecule(s) found.\n";

	string sdim(argv[5]);
	if (sdim == "xy")
		cout << "REMARK Calculate MSD in X-Y plane\n";
	else if (sdim == "xyz")
		cout << "REMARK Calculate MSD in 3D space\n";
	else
		die("error: unknown argument \"" + sdim + "\"");

	cout << "REMARK============\n";

	vector< vector<double> > _x(idx.size());
	vector< vector<double> > _y(idx.size());
	vector< vector<double> > _z(idx.size());

	while (DCD.read_1step())
	{
		// loop over selected molecule
		int icnt = 0;
		for (auto& i: idx)
		{
			Molecule& mol = PSFFile.moleculeVector[i - 1];

			mol.calc_com();

			_x[icnt].push_back(mol.vcom.x());
			_y[icnt].push_back(mol.vcom.y());
			_z[icnt].push_back(mol.vcom.z());

			++icnt;
		}
	}

	cout << setprecision(8) << scientific;
	int delta = atoi(argv[4]);
	for (int idel = 1; idel <= delta; idel++)
	{
		double vmsd = 0;
		for (int istep = 1; istep <= DCD._nsets() - idel; istep++)
		{
			for (int m = 0; m < idx.size(); m++)
			{
				double xdel
				= _x[m][istep + idel - 1] - _x[m][istep - 1];
				double ydel
				= _y[m][istep + idel - 1] - _y[m][istep - 1];

				vmsd += xdel * xdel + ydel * ydel;

				if (sdim == "xyz")
				{
					double zdel =
					_z[m][istep + idel - 1] - _z[m][istep - 1];
					vmsd += zdel * zdel;
				}
			}
		}
		vmsd /= DCD._nsets() - idel;
		vmsd /= idx.size(); // average over molecules

		cout 
		<< setw(10) << idel
		<< setw(16) << vmsd
		<< '\n';
	}
}
