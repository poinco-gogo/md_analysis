#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include "common.hpp"
#include "ComputeHistogram.hpp"
using namespace std;
int main(int argc, char** argv)
{
	if (argc < 6)
	{
		cout 
		<< "\nCalculate weighted 1d histogram.\n"
		<< "\nusage: ./a.out metadatafile vmin vmax nbin temperature\n\n";
		return 1;
	}

	output_args(argc, argv);

	double vmin        = atof(argv[2]);
	double vmax        = atof(argv[3]);
	unsigned int nbin  = atoi(argv[4]);
	double temperature = atof(argv[5]);
	double kbT         = BOLTZMAN * temperature;
	bool normalize     = false;

	ComputeHistogram JOB(vmin, vmax, nbin, normalize);
	string s;
	ifstream fi(argv[1]);
	while (getline(fi, s))
	{
		istringstream is(s);
		string stmp;
		is >> stmp;
		JOB.load_data(stmp);
		is >> stmp;
		JOB.load_weight(stmp);
	}

	JOB.calc_weighted_histogram();

	JOB.output_pmf(kbT);
}
