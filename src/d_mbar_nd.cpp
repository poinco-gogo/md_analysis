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
		<< "\nusage: ./a.out metadatafile ndim nbias tol temperature ofilename nself\n"
		<< "\n   outputfile:  ofilename.fene\n"
		<<   "                ofilename{}.weight\n"
		<<   "                ofilename.pmf\n\n";
		return 1;
	}

	output_args(argc, argv);

	unsigned int ndim  = atoi(argv[2]);
	unsigned int nbias = atoi(argv[3]);
	double tol         = atof(argv[4]);
	double temperature = atof(argv[5]);
	unsigned int nself = atoi(argv[7]);
	string speriod = "no";

	ComputeMBAR JOB(argv[1], ndim, nbias, tol, temperature, argv[6], nself, speriod);

	do {
		JOB.mbar_iteration();
	} while (!JOB.check_convergence());

	JOB.calc_unbiasing_weights();

	JOB.output_results();
}
