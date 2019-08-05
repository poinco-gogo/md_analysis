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
		"\nusage: ./a.out psf dcd\n\n";
		return 1;
	}
	vector<Atom> atomVector;
	PSF PSFFile(argv[1], &atomVector);
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
