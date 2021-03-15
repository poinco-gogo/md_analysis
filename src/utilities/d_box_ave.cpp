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
	if (argc < 2)
	{
		cout << "\n./a.out metadatafile\n\n";
		return 1;
	}

	vector<double> boxx, boxy;

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
	}

	cout << "REMARK Total sample size: " << boxx.size() << '\n';

	double sumx = 0;
	double sumy = 0;
	for (int i = 0; i < boxx.size(); i++)
	{
		sumx += boxx[i];
		sumy += boxy[i];
	}

	double denom = static_cast<double>(boxx.size());
	sumx /= denom;
	sumy /= denom;

	cout
		<< "REMARK Results =========================\n"
		<< setprecision(8) << scientific
		<< '\n'
		<< setw(16) << sumx
		<< setw(16) << sumy
		<< '\n';
}
