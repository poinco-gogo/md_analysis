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
		cout << "\nusage: ./a.out file min max bin col band_width\n\n";
		return 0;
	}

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

	vector<double> data;

	string s;
	while ( getline(fi, s) )
	{
		if (s.empty() || is_comment(s))
			continue;

		istringstream is(s);
		double val;
		for (int i = 0; i < ncol; i++)
			is >> val;
		data.push_back(val);
	}

	ComputeKDE JOB(&data, h);

	cout << setprecision(4) << scientific;
	for (int i = 0; i < nbin; i++)
	{
		double x = vmin + w * i;

		cout
			<< setw(12) << x
			<< setw(16) << JOB.estimate_gauss( x )
			<< '\n';
	}
}
