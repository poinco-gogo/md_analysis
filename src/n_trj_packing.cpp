#include <iostream>
#include <iomanip>
#include <cstdlib>
#include "PSF.hpp"
#include "Index.hpp"
#include "ReadDCD.hpp"
using namespace std;
int main(int argc, char** argv)
{
	if (argc < 5)
	{
		cout << 
		"\nPACKING SCORE CALCULATION\n"
		"\nusage: ./a.out psf dcd ind1 ind2\n\n";
		return 1;
	}
	cout << "REMARK ";
	for (int i = 0; i < argc; i++)
		cout << argv[i] << " ";
	cout << "\nREMARK PACKING SCORE CALCULATION\n";
	vector<Atom> atomVector;
	PSF PSFFile(argv[1], &atomVector);
	ReadDCD DCD(&atomVector);
	if (!DCD.open_dcd_read(argv[2]))
		return 0;

	Index job1(argv[3]);
	Index job2(argv[4]);
	vector<int> ind1, ind2;
	if (!job1.load_indexes(&ind1)) return 0;
	if (!job2.load_indexes(&ind2)) return 0;

	cout << "REMARK " << ind1.size() << " indexes found.\n";
	cout << "REMARK " << ind2.size() << " indexes found.\n";

	const int num_pairs = ind1.size() * ind2.size();

	cout << "REMARK============\n";
	
	cout << setprecision(8) << scientific;
	while (DCD.read_1step())
	{
		double sum = 0;

		for (auto& i: ind1)
		{
			Atom& iat = atomVector[i];

			for (auto& j: ind2)
			{
				Atom& jat = atomVector[j];

				double d2
				= (iat.position - jat.position).squaredNorm();

				double d6 = d2 * d2 * d2;

				sum += 1./ d6;
			}
		}

		sum /= num_pairs;

		cout 
		<< setw(10) << DCD._nsteps()
		<< setw(16) << sum
		<< '\n';
	}
}
