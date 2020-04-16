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
		"\nusage: ./a.out psf dcd resid delta(in MD steps) dimension('xy' or 'xyz')\n\n";
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

	int tgtID = atoi(argv[3]);
	vector<int> idx;
	for (auto& atom: atomVector)
	{
		if (atom.PSFResID == tgtID)
			idx.push_back(atom.PSFIndex);
	}
	cout << "REMARK " << idx.size() << " target atom(s) found.\n";

	string sdim(argv[5]);
	if (sdim == "xy")
		cout << "REMARK Calculate MSD in X-Y plane\n";
	else if (sdim == "xyz")
		cout << "REMARK Calculate MSD in 3D space\n";
	else
		die("error: unknown argument \"" + sdim + "\"");

	cout << "REMARK============\n";
	
	vector<double> _x, _y, _z;

	while (DCD.read_1step())
	{
		Eigen::Vector3d com = V3ZERO;

		double totmass = 0.0;
		for (auto& i: idx)
		{
			Atom& at = atomVector[i - 1];

			com += at.position * at.mass;

			totmass += at.mass;
		}
		com /= totmass;

		_x.push_back(com.x());
		_y.push_back(com.y());
		_z.push_back(com.z());
	}

	cout << setprecision(8) << scientific;
	int delta = atoi(argv[4]);
	for (int idel = 1; idel <= delta; idel++)
	{
		double vmsd = 0;
		for (int istep = 1; istep <= _x.size() - idel; istep++)
		{
			double xdel = _x[istep + idel - 1] - _x[istep - 1];
			double ydel = _y[istep + idel - 1] - _y[istep - 1];

			vmsd += xdel * xdel + ydel * ydel;

			if (sdim == "xyz")
			{
				double zdel =
					_z[istep + idel - 1] - _z[istep - 1];
				vmsd += zdel * zdel;
			}
		}
		vmsd /= _x.size() - idel;

		cout 
		<< setw(10) << idel
		<< setw(16) << vmsd
		<< '\n';
	}
}
