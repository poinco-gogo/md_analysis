#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include "ComputeKDE.hpp"
#include "FileIO.hpp"
#include "common.hpp"
using namespace std;
int main (int argc, char** argv)
{
	if (argc < 7)
	{
		cout <<
			"\nKernel Density Estimator for 1D valuable\n"
			"\nusage: ./a.out file min max bin col band_width\n\n";
		return 0;
	}
	output_args(argc, argv);

	istringstream ismin(argv[2]);
	istringstream ismax(argv[3]);
	istringstream isbin(argv[4]);
	istringstream iscol(argv[5]);
	istringstream isbnd(argv[6]);
	double vmin, vmax, h;
	int nbin, ncol;
	ismin >> vmin; ismax >> vmax;
	isbin >> nbin; iscol >> ncol;
	isbnd >> h;

	double w = (vmax - vmin) / nbin;

	vector<double> data;

	FileIO fi(argv[1]);
	if (!fi.load_data(ncol, &data)) return 1;

	ComputeKDE JOB(&data, h);

	cout << setprecision(4) << scientific;
	for (int i = 0; i < nbin; i++)
	{
		double x = vmin + w / 2. * (2 * i + 1);

		cout
			<< setw(12) << x
			<< setw(16) << JOB.estimate_gauss( x )
			<< '\n';
	}
}
