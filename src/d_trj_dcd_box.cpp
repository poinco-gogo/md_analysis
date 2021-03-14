#include <iostream>
#include <iomanip>
#include <cstdlib>
#include "PSF.hpp"
#include "ReadDCD.hpp"
using namespace std;
int main(int argc, char** argv)
{
	if (argc < 3)
	{
		cout << 
		"\nusage: ./a.out psf[natom] dcd\n\n";
		return 1;
	}

	int natom = 0;
	vector<Atom> atomVector;

	string spsf(argv[1]);
	string spsfex = spsf.substr(spsf.find_last_of(".") + 1);
	if (spsfex == "psf")
	{
		PSF PSFFile(argv[1]);
		atomVector = PSFFile.atomVector;
	}
	else
	{
		natom = atoi(argv[1]);
		atomVector.resize(natom);
	}

	ReadDCD DCD(&atomVector);
	if (!DCD.open_dcd_read(argv[2]))
		return 0;

	cout << setprecision(8) << scientific;
	while (DCD.read_1step())
	{
		cout 
		<< setw(10) << DCD._nsteps()
		<< setw(16) << DCD._boxx()
		<< setw(16) << DCD._boxy()
		<< setw(16) << DCD._boxz()
		<< '\n';
	}
}
