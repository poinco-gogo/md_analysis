#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include "PSF.hpp"
#include "PDB.hpp"
#include "ComputeSASA.hpp"

using namespace std;

int main(int argc, char** argv)
{
	if (argc != 6)
	{
		cout <<
		"\nSASA calculation\n"
		"\nusage: ./a.out psf pdb npoints resid mode\n"
		"\n   mode = 1: gaussian random numbers\n"
		  "   mode = 2: Fibonacci sphere algorithm\n\n";
		return 0;
	}

	output_args(argc, argv);

	PSF PSFFile(argv[1]);

	PDB PDBFile(argv[2]);
	if (!PDBFile.LoadCoords(PSFFile.atomVector))
		return 0;

	unsigned int npoints = atoi(argv[3]);
	unsigned int resid   = atoi(argv[4]);
	unsigned int mode    = atoi(argv[5]);

	ComputeSASA JOB(npoints, mode, &PSFFile.atomVector);

	cout
		<< setprecision(2) << fixed
		<< setw(16) << JOB.calc_sasa(resid)
		<< '\n';
}
