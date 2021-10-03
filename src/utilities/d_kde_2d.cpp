#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <cmath>
#include <limits>
#include "ComputeKDE.hpp"
#include "FileIO.hpp"
#include "common.hpp"
using namespace std;
int main (int argc, char** argv)
{
	if (argc < 3)
	{
		cout <<
			"\nKernel Density Estimator for 2D valuables\n"
			"\nusage: ./a.out metafile temperature\n\n";
		return 0;
	}

	output_args(argc, argv);

	double temperature = atof(argv[2]);
	double kbT = temperature * BOLTZMAN;

	double min1, max1, band_width1;
	double min2, max2, band_width2;
	int    nbin1, nbin2;

	ifstream fi(argv[1]);
	string s;
	getline(fi, s); // min1 max1 nbin1 band_width1
	istringstream is1(s);
	is1 >> min1 >> max1 >> nbin1 >> band_width1;
	getline(fi, s); // min2 max2 nbin2 band_width2
	istringstream is2(s);
	is2 >> min2 >> max2 >> nbin2 >> band_width2;
	double w1 = (max1 - min1) / nbin1;
	double w2 = (max2 - min2) / nbin2;

	vector<double> data1, data2;

	while (getline(fi, s))
	{
		istringstream is(s);
		string val;
		is >> val;
		FileIO fi(val);
		if (!fi.load_data(2, &data1)) return 1;
		is >> val;
		FileIO fj(val);
		if (!fj.load_data(2, &data2)) return 1;
	}

	if (data1.size() != data2.size())
		die("data1.size() != data2.size()");

	ComputeKDE JOB(&data1, &data2,
			min1, max1, nbin1, band_width1,
			min2, max2, nbin2, band_width2
		);

	vector<double> bincenters1, bincenters2;

	vector< vector<double> > pmf;

	cout << setprecision(4) << scientific;
	for (int i = 0; i < nbin1; i++)
	{
		double x = min1 + w1 / 2. * (2 * i + 1);

		bincenters1.push_back( x );

		vector<double> vtmp;

		for (int j = 0; j < nbin2; j++)
		{
			double y = min2 + w2 / 2. * (2 * j + 1);

			bincenters2.push_back( y );

			vtmp.push_back( -kbT * log( JOB.estimate_gauss( x, y ) ) );
		}

		pmf.push_back(vtmp);
	}

	double pivot = 1e10;
	for (int i = 0; i < nbin1; i++)
	{
		for (int j = 0; j < nbin2; j++)
		{
			if ( isfinite(pmf[i][j]) && pmf[i][j] < pivot)
				pivot = pmf[i][j];
		}
	}

	ofstream fo1("ax1.mat");
	ofstream fo2("ax2.mat");
	ofstream fo3("tmp.mat");

	fo1 << setprecision(8) << scientific;
	fo2 << setprecision(8) << scientific;
	fo3 << setprecision(8) << scientific;

	for (int i = 0; i < nbin1; i++)
		fo1 << setw(16) << bincenters1[i];

	for (int i = 0; i < nbin2; i++)
		fo2 << setw(16) << bincenters2[i];

	float inf = std::numeric_limits<float>::infinity();
	for (int i = 0; i < nbin1; i++)
	{
		for (int j = 0; j < nbin2; j++)
		{
			fo3 << setw(16) << pmf[i][j] - pivot;
		}
		fo3 << '\n';
	}

/*
	for (int i = 0; i < nbin1; i++)
	{
		for (int j = 0; j < nbin2; j++)
		{
		cout
			<< setw(16) << bincenters1[i]
			<< setw(16) << bincenters2[j]
			<< setw(16) << pmf[i][j] - pivot
			<< '\n';
		}
	}
*/
}
