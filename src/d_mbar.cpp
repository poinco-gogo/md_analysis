#include <iostream>
#include <sstream>
#include <cstdlib>
#include "common.hpp"
#include "ComputeMBAR.hpp"
using namespace std;
int main(int argc, char** argv)
{
	if (argc < 5)
	{
		cout 
		<< "\nMBAR analysis of umbrella sampling data.\n"
		<< "\nusage: ./a.out metadatafile ncvs tol temperature [P/Ppi]\n\n";
		return 1;
	}

	output_args(argc, argv);

	int    ncvs        = atoi(argv[2]);
	double tol         = atof(argv[3]);
	double temperature = atof(argv[4]);
	string speriod = "no";
	if (argv[7])
	{
		speriod = argv[7];
	}

	ComputeMBAR JOB(argv[1], ncvs, tol, temperature, speriod);

	do {
		JOB.mbar_iteration();
	} while (!JOB.check_convergence());
	
	JOB.output_results();

}
