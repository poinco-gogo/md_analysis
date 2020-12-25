#include <iostream>
#include <sstream>
#include <cstdlib>
#include "common.hpp"
#include "ComputeMBAR.hpp"
using namespace std;
int main(int argc, char** argv)
{
	if (argc < 9)
	{
		cout 
		<< "\nMBAR analysis of umbrella sampling data.\n"
		<< "\nusage: ./a.out metadatafile ndim vmin vmax nbin tol temperature ofilename [P/Ppi]\n"
		<< "\n   outputfile:  ofilename.fene\n"
		<<   "                ofilename{}.weight\n"
		<<   "                ofilename.pmf\n\n";
		return 1;
	}

	output_args(argc, argv);

	unsigned int ndim  = atoi(argv[2]);
	double vmin        = atof(argv[3]);
	double vmax        = atof(argv[4]);
	unsigned int nbin  = atoi(argv[5]);
	double tol         = atof(argv[6]);
	double temperature = atof(argv[7]);
	string speriod = "no";
	if (argv[9])
	{
		speriod = argv[9];
	}

	if (ndim > 1)
	{
		cerr << "\nndim > 1 is not supported yet:)\n\n";
		return 0;
	}

	ComputeMBAR JOB(argv[1], ndim, vmin, vmax, nbin, tol, temperature, argv[8], speriod);

	do {
		JOB.mbar_iteration();
	} while (!JOB.check_convergence());
	
	JOB.calc_unbiasing_weights();

	JOB.output_results();
}
