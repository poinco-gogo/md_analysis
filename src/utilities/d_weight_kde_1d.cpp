#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <cmath>
#include "ComputeKDE.hpp"
#include "FileIO.hpp"
#include "common.hpp"
using namespace std;
int main (int argc, char** argv)
{
	if (argc < 8)
	{
		cout <<
			"\nWeighted Kernel Density Estimator for 1D valuable\n"
			"\nusage: ./a.out metafile min max bin band_width temperature cutoff\n\n";
		return 0;
	}

	output_args(argc, argv);

	istringstream ismin(argv[2]);
	istringstream ismax(argv[3]);
	istringstream isbin(argv[4]);
	istringstream isbnd(argv[5]);
	double vmin, vmax, h;
	int nbin;
	ismin >> vmin; ismax >> vmax;
	isbin >> nbin;
	isbnd >> h;
	double temperature = atof(argv[6]);
	double kbT = temperature * BOLTZMAN;
	double cutoff = atof(argv[7]);

	double w = (vmax - vmin) / nbin;

	vector<double> data, weight, distance;

	ifstream fi(argv[1]);
	string s;
	while (getline(fi, s))
	{
		istringstream is(s);
		string val;
		is >> val;
		FileIO fi(val);
		if (!fi.load_data(2, &data)) return 1;
		is >> val;
		FileIO fj(val);
		if (!fj.load_data(2, &weight)) return 1;
		is >> val;
		FileIO fk(val);
		if (!fk.load_data(2, &distance)) return 1;
	}

	ComputeKDE JOB(&data, &weight, &distance, h, cutoff);

	vector<double> bincenters, prob;

	cout << setprecision(4) << scientific;
	for (int i = 0; i < nbin; i++)
	{
		double x = vmin + w / 2. * (2 * i + 1);

		bincenters.push_back( x );

		prob.push_back( -kbT * log( JOB.estimate_gauss_weight( x ) ) );
	}

	double pivot = 9999;
	for (int i = 0; i < nbin; i++)
	{
		if ( isfinite(prob[i]) && prob[i] < pivot)
			pivot = prob[i];
	}

	for (int i = 0; i < nbin; i++)
	{
		cout
			<< setw(16) << bincenters[i]
			<< setw(16) << prob[i] - pivot
			<< '\n';
	}
}
