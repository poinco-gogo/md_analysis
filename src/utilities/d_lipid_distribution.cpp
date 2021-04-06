#include <iostream>
#include <fstream>
#include <sstream>
#include "PSF.hpp"
#include "ReadDCD.hpp"
#include "common.hpp"
#include "Lattice.hpp"
using namespace std;

Eigen::Vector3d getcom(string filename);
double          getboxsize(string filename);

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
	Eigen::Vector3d com = getcom(s);
	cout << "REMARK Com coordinate: " << com.transpose() << '\n';
	getline(fi, s);
	double boxx = getboxsize(s);
	double boxy = boxx;
	cout << "REMARK Box size: " << boxx << '\n';
	double xmin = -0.5 * boxx + com.x();
	double xmax =  0.5 * boxx + com.x();
	double ymin = -0.5 * boxy + com.y();
	double ymax =  0.5 * boxy + com.y();
	int nbin = atoi(argv[3]);
	double w    = boxx / nbin;
	vector<double> xbincenters(nbin, 0.);
	vector<double> ybincenters(nbin, 0.);
	for (int i = 0; i < nbin; i++)
	{
		xbincenters[i] = xmin + w / 2. * (2. * i + 1);
		ybincenters[i] = ymin + w / 2. * (2. * i + 1);
	}

	vector< vector<double> > mesh1, mesh2;
	for (int i = 0; i < nbin; i++)
	{
		vector<double> vtmp(nbin, 0.);
		mesh1.push_back(vtmp);
		mesh2.push_back(vtmp);
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
			const double lx = TRJFile._boxx();
			const double ly = TRJFile._boxy();
			const double lz = TRJFile._boxz();

			Lattice lattice(lx, ly, lz);

			for (int i = 0; i < natom; i++)
			{
				Eigen::Vector3d& lp = atomVector[i].position;
				Eigen::Vector3d del = lp - com;
				Eigen::Vector3d del2 =  lattice.wrap_delta(del);
				lp.x() += del2.x();
				lp.y() += del2.y();

				bool goto_next = false;

				for (int x = 0; x < nbin; x++)
				{
					for (int y = 0; y < nbin; y++)
					{
						if ( abs( lp.x() - xbincenters[x] ) < w * 0.5 && abs( lp.y() - ybincenters[y] ) < w * 0.5 )
						{
							if ( lp.z() > 0)
								mesh1[x][y] += 1.0;
							else
								mesh2[x][y] += 1.0;
							goto_next = true;
							break;
						}
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
			mesh1[x][y] /= (tot_step*w*w);
			mesh2[x][y] /= (tot_step*w*w);
		}
	}

	ofstream fo1("ax1.mat");
	ofstream fo2("ax2.mat");
	ofstream fo3("ms1.mat");
	ofstream fo4("ms2.mat");

	fo1 << setprecision(8) << scientific;
	fo2 << setprecision(8) << scientific;
	fo3 << setprecision(8) << scientific;

	for (int i = 0; i < nbin; i++)
		fo1 << setw(16) << xbincenters[i];

	for (int i = 0; i < nbin; i++)
		fo2 << setw(16) << ybincenters[i];

	for (int i = 0; i < nbin; i++)
	{
		for (int j = 0; j < nbin; j++)
		{
			fo3 << setw(16) << mesh1[i][j];
			fo4 << setw(16) << mesh2[i][j];
		}
		fo3 << '\n';
		fo4 << '\n';
	}
}
Eigen::Vector3d getcom(string filename)
{
	ifstream fi(filename.c_str());

	string s;
	for (int i = 0; i < 11; i++)
		getline(fi, s);

	double x, y, z;
	istringstream is(s);
	is >> x >> y >> z;

	Eigen::Vector3d vtmp(x, y, z);

	return vtmp;
}

double getboxsize(string filename)
{
	ifstream fi(filename.c_str());

	string s;
	for (int i = 0; i < 4; i++)
		getline(fi, s);

	double boxx;
	istringstream is(s);
	is >> boxx;

	return boxx;
}
