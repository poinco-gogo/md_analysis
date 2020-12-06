#include <iostream>
#include <sstream>
#include <cstdlib>
#include "common.hpp"
#include "ComputeMBAR.hpp"
using namespace std;
int main(int argc, char** argv)
{
	if (argc < 8)
	{
		cout 
		<< "\nMBAR analysis of umbrella sampling data.\n"
		<< "\nusage: ./a.out metadatafile ndim vmin vmax nbin tol temperature [P/Ppi]\n\n";
		return 1;
	}

	output_args(argc, argv);

	int    ndim        = atoi(argv[2]);
	double vmin        = atof(argv[3]);
	double vmax        = atof(argv[4]);
	int    nbin        = atoi(argv[5]);
	double tol         = atof(argv[6]);
	double temperature = atof(argv[7]);
	string speriod = "no";
	if (argv[8])
	{
		speriod = argv[8];
	}

	ComputeMBAR JOB(argv[1], ndim, vmin, vmax, nbin, tol, temperature, speriod);

	do {
		JOB.mbar_iteration();
	} while (!JOB.check_convergence());
	
	JOB.calc_weights();

	JOB.output_results();
}