#include <iostream>
#include <iomanip>
#include <cstdlib>
#include "PSF.hpp"
#include "ReadDCD.hpp"
using namespace std;
int main(int argc, char** argv)
{
	if (argc < 5)
	{
		cout << 
		"\nGET DISTANCE BETWEEN 2 ATOMS FROM DCD\n"
		"\nusage: ./a.out psf dcd id1 id2\n\n";
		return 1;
	}
	cout << "REMARK ";
	for (int i = 0; i < argc; i++)
		cout << argv[i] << " ";
	cout << "\nREMARK GET DISTANCE FROM DCD.\n";
	vector<Atom> atomVector;
	PSF PSFFile(argv[1], &atomVector);
	ReadDCD DCD(&atomVector);
	if (!DCD.open_dcd_read(argv[2]))
		return 0;

	int id1 = atoi(argv[3]);
	int id2 = atoi(argv[4]);

	Atom& atom1 = atomVector[id1 - 1];
	Atom& atom2 = atomVector[id2 - 1];

	cout << "REMARK atom1:";
	PSFFile.showAtom(atom1);
	
	cout << "REMARK atom2:";
	PSFFile.showAtom(atom2);

	cout << "REMARK============\n";
	
	cout << setprecision(8) << scientific;
	while (DCD.read_1step())
	{
		cout 
		<< setw(10) << DCD._nsteps()
		<< setw(16) << (atom1.position - atom2.position).norm()
		<< '\n';
	}
}
