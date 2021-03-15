#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <iomanip>
#include <cstdlib>
#include "FileIO.hpp"

using namespace std;
int main (int argc, char** argv)
{
	if (argc < 3)
	{
		cout << "\n./a.out metadatafile cutoff\n\n";
		return 1;
	}

	double cutoff = atof(argv[2]);

	vector<double> boxx, boxy, weight, distance;

	ifstream fi(argv[1]);
	string s;
	while (getline(fi, s))
	{
		istringstream is(s);
		string val;
		is >> val;
		FileIO fj(val);
		if (!fj.load_data(2, &boxx)) return 1;
		FileIO fk(val);
		if (!fk.load_data(3, &boxy)) return 1;
		is >> val;
		FileIO fl(val);
		if (!fl.load_data(2, &weight)) return 1;
		is >> val;
		FileIO fm(val);
		if (!fm.load_data(2, &distance)) return 1;
	}

	cout << "REMARK Total sample size: " << boxx.size() << '\n';

	double sumx = 0;
	double sumy = 0;
	for (int i = 0; i < boxx.size(); i++)
	{
		if (distance[i] > cutoff) continue;

		sumx += boxx[i] * weight[i];
		sumy += boxy[i] * weight[i];
	}

	cout
		<< "REMARK Results =========================\n"
		<< setprecision(8) << scientific
		<< '\n'
		<< setw(16) << sumx
		<< setw(16) << sumy
		<< '\n';
}
