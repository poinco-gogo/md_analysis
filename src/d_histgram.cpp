#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <vector>
#include "ComputeHistgram.hpp"
#include "common.hpp"
using namespace std;
int main (int argc, char** argv)
{	
	if (argc < 7)
	{
		cout << "\nD_HISTGRAM\n"
			"\nusage: ./a.out file min max bin col normalize?\n\n";
		return 1;
	}
	istringstream ismin(argv[2]);
	istringstream ismax(argv[3]);
	istringstream isbin(argv[4]);
	istringstream iscol(argv[5]);
	string        snorm(argv[6]);
	double vmin, vmax;
	int nbin, ncol;
	ismin >> vmin; ismax >> vmax;
	isbin >> nbin; iscol >> ncol;
	bool normalize = snorm == "yes" ? true : false;

	vector<double> data;

	string s;
	ifstream fi(argv[1]);
	while (getline(fi, s))
	{
		if (s.empty() || is_comment(s))
			continue;

		istringstream is(s);
		double val;
		for (int i = 0; i < ncol; i++)
			is >> val;

		data.push_back( val );
	}

	ComputeHistgram JOB(&data, vmin, vmax, nbin, normalize);

	JOB.calc_histgram();

	JOB.output();
}
