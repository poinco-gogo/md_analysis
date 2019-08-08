#include <Eigen/Core>
#include <iostream>
#include <cstdlib>
#include <iomanip>
#include "PSF.hpp"
#include "PDB.hpp"
using namespace std;
int main(int argc, char** argv)
{
	if (argc < 5)
	{
		cout << 
		"\nGET DISTANCE BETWEEN 2 ATOMS.\n"
		"\nusage: ./a.out psf pdb id1 id2 [need unit vector?]\n\n";
		return 1;
	}
	cout << "REMARK ";
	for (int i = 0; i < argc; i++)
		cout << argv[i] << " ";
	cout << "\nREMARK GET DISTANCE\n";
	vector<Atom> atomVector;
	PSF PSFFile(argv[1], &atomVector);
	PDB PDBFile(argv[2]);
	if (!PDBFile.LoadCoords(atomVector))
		return 1;

	int id1 = atoi(argv[3]);
	int id2 = atoi(argv[4]);

	Atom& atom1 = atomVector[id1 - 1];
	Atom& atom2 = atomVector[id2 - 1];

	cout << "REMARK atom1:";
	PSFFile.showAtom(atom1);
	
	cout << "REMARK atom2:";
	PSFFile.showAtom(atom2);

	cout << "REMARK============\n";
	Eigen::Vector3d r12 = atom1.position - atom2.position;
	double len = r12.norm();
	cout << setprecision(8) << scientific;
	cout << "REMARK r12       = ";
	cout << r12.transpose() << '\n';
	cout << "REMARK |r12|     = " << len << '\n';
	if (argv[5])
	{
		string s(argv[5]);
		if (s != "yes" && s != "true" && s != "on")
		{
			cerr << "unknown option \"" << argv[5] << "\" detected.\n";
			return 1;
		}
		cout << "REMARK r12/|r12| = ";
		cout << (r12 / len).transpose() << '\n';
		cout << "REMARK |r12/|r12|| = "
	     	<< (r12 / len).norm()
	     	<< '\n';
	}
}
