#include <iostream>
#include <iomanip>
#include "PSF.hpp"
#include "Index.hpp"
using namespace std;
int main (int argc, char** argv)
{
	if (argc < 2)
	{
		cout << "\nusage: ./a.out psf ind\n\n";
		return 1;
	}

	PSF PSFFile(argv[1]);

	vector<int> vind;
	Index IDX(argv[2]);
	if (!IDX.load_indexes(&vind))
		return 1;

	int icnt = 0;
	for (auto& i: vind)
	{
		Atom& atom = PSFFile.atomVector[i - 1];

		cout << setw(8) << ++icnt << setw(6) << atom.PSFSegmentName
			<< "   ";
		PSFFile.showAtom(atom);
	}
}
