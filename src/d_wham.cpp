#include <iostream>
#include <sstream>
#include <cstdlib>
#include "ComputeWHAM.hpp"
using namespace std;
int main(int argc, char** argv)
{
	if (argc < 2)
	{
		cout 
		<< "\nD_WHAM\n"
		<< "\nusage: ./a.out metadatafile min max nbin tol temperature [P/Ppi]\n\n";
		return 1;
	}
	cout << "REMARK ";
	for (int i = 0; i < argc; i++)
		cout << argv[i] << " ";
	cout << '\n';

	istringstream ismin(argv[2]);
	istringstream ismax(argv[3]);
	istringstream isbin(argv[4]);
	istringstream istol(argv[5]);
	istringstream istmp(argv[6]);
	double vmin, vmax, tol, temperature;
	int nbin;
	ismin >> vmin; ismax >> vmax;
	istol >> tol;  istmp >> temperature;
	isbin >> nbin;
	string speriod = "no";
	if (argv[7])
	{
		speriod = argv[7];
	}

	ComputeWHAM JOB(argv[1], vmin, vmax, nbin, tol, temperature, speriod);

	do {
		JOB.wham_iteration();
	} while (!JOB.check_convergence());
	
	JOB.output_results();
}
