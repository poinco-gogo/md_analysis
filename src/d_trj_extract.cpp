#include <iostream>
#include <iomanip>
#include <cstdlib>
#include "PSF.hpp"
#include "ReadDCD.hpp"
#include "Index.hpp"
using namespace std;
int main(int argc, char** argv)
{
	if (argc < 4)
	{
		cout << 
		"\nusage: ./a.out psf dcd ind\n\n";
		return 1;
	}
	cout << "REMARK ";
	for (int i = 0; i < argc; i++)
		cout << argv[i] << ' ';
	cout << '\n';

	vector<Atom> atomVector;

	PSF PSFFile(argv[1], &atomVector);

	ReadDCD DCD(&atomVector);
	if (!DCD.open_dcd_read(argv[2]))
		return 0;

	Index idx(argv[3]);
	vector<int> cvs;
	if (!idx.load_indexes(&cvs))
		return 0;
	cout << "REMARK " << cvs.size() << " indexes loaded.\n";

	cout << setprecision(8) << scientific;
	while (DCD.read_1step())
	{
		cout << setw(16) << DCD._nsteps();

		for (auto& i: cvs)
		{
			Atom& atom = atomVector[i - 1];

			cout
				<< setw(16) << atom.position.x()
				<< setw(16) << atom.position.y()
				<< setw(16) << atom.position.z();
		}
		cout << '\n';
	}
}
