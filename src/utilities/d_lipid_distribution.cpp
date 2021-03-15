#include <iostream>
#include <fstream>
#include <sstream>
#include "PSF.hpp"
#include "ReadDCD.hpp"
#include "common.hpp"
using namespace std;
int main(int argc, char** argv)
{
	if (argc != 4)
	{
		cout << 
		"\nD_LIPID_DISTRIBUTION\n"
		"\nCalculate distribution of lipid\n"
		"\nusage: ./a.out natom metadatafile nbin\n\n";
		return 1;
	}

	output_args(argc, argv);

	int natom = 0;
	vector<Atom> atomVector;

	natom = atoi(argv[1]);
	atomVector.resize(natom);

	ifstream fi(argv[2]);
	string s;
	getline(fi, s);
	istringstream is(s);
	double boxx, boxy;
	is >> boxx, boxy;
	double vmin = -0.5 * boxx;
	double vmax =  0.5 * boxx;
	int nbin = atoi(argv[3]);
	double w    = (vmax - vmin) / nbin;
	vector<double> bincenters(nbin, 0.);
	for (int i = 0; i < nbin; i++)
		bincenters[i] = vmin + w / 2. * (2. * i + 1);

	vector< vector<double> > mesh;
	for (int i = 0; i < nbin; i++)
	{
		vector<double> vtmp(nbin, 0.);
		mesh.push_back(vtmp);
	}

	double tot_step = 0;
	while (getline(fi, s))
	{
		istringstream is(s);
		string val;
		is >> val;
		ReadDCD TRJFile(&atomVector);
		if (!TRJFile.open_dcd_read(val))
			return 0;
	
		while(TRJFile.read_1step())
		{
			for (int i = 0; i < natom; i++)
			{
				bool goto_next = false;

				for (int x = 0; x < nbin; x++)
				{
					for (int y = 0; y < nbin; y++)
					{
						if ( abs( atomVector[i].position.x() - bincenters[x] ) < w * 0.5 && abs( atomVector[i].position.y() - bincenters[y] ) < w * 0.5 )
						mesh[x][y] += 1.0;
						goto_next = true;
						break;
					}

					if (goto_next) break;
				}
			}
		}

		tot_step += static_cast<double>(TRJFile._nsteps());
	}

	cout << "REMARK Total step: " << tot_step << '\n';

	//  average over step
	for (int x = 0; x < nbin; x++)
	{
		for (int y = 0; y < nbin; y++)
		{
			mesh[x][y] /= tot_step;
		}
	}

	ofstream fo1("ax1.mat");
	ofstream fo2("ax2.mat");
	ofstream fo3("tmp.mat");

	fo1 << setprecision(8) << scientific;
	fo2 << setprecision(8) << scientific;
	fo3 << setprecision(8) << scientific;

	for (int i = 0; i < nbin; i++)
		fo1 << setw(16) << bincenters[i];

	for (int i = 0; i < nbin; i++)
		fo2 << setw(16) << bincenters[i];

	for (int i = 0; i < nbin; i++)
	{
		for (int j = 0; j < nbin; j++)
		{
			fo3 << setw(16) << mesh[i][j];
		}
		fo3 << '\n';
	}
}
