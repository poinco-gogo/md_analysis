#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include "ComputeKDE.hpp"
#include "common.hpp"
using namespace std;
int main (int argc, char** argv)
{
	if (argc < 7)
	{
		cout << "\nKernel Density Estimator for 2D valuables\n"
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

	ifstream fi(argv[1]);

	vector<double> xdata, ydata;

	string s;
	while ( getline(fi, s) )
	{
		istringstream is(s);
		double xval, yval;
		is >> xval, yval;
		xdata.push_back(xval);
		ydata.push_back(yval);
	}

	ComputeKDE JOB(&xdata, &ydata, h);

	cout << setprecision(4) << scientific;
	for (int i = 0; i < nbin; i++)
	{
		double x = vmin + w / 2. * (2 * i + 1);

		for (int j = 0; j < nbin; j++)
		{
			double y = vmin + w / 2. * (2 * j + 1);

			cout
			<< setw(12) << x
			<< setw(12) << y
			<< setw(16) << JOB.estimate_gauss( x , y )
			<< '\n';
		}
	}
}
