#include <iostream>
#include <fstream>
#include "PSF.hpp"
using namespace std;
int main (int argc, char** argv)
{
	if (argc < 1)
	{
		cout << "\n./a.out psf\n\n";
		return 1;
	}

	PSF PSFFile(argv[1]);

	vector<int> lipid_p;

	int icnt = 0;
	for (auto& atom: PSFFile.atomVector)
	{
		++icnt;

		if (atom.PDBAtomName == "P")
			lipid_p.push_back(icnt);
	}

	ofstream fp("p.ind");

	for (auto& id: lipid_p)
	{
		fp << id << " ";
	}

	fp.close();
}
