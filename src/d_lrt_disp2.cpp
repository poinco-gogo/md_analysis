#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <vector>
#include <Eigen/Core>
#include "MatrixUtil.hpp"
using namespace std;
int main(int argc, char** argv)
{
	if (argc != 4)
	{
		cout <<
		"\n\nusage: ./a.out natom vcv force\n\n";
		return 1;
	}

	double kBT = 0.6;
	
	int natom = atoi(argv[1]);

	// Load matrix file
	double** MATRIX;
	string filename(argv[2]);
	MatrixUtil JOB(filename, &MATRIX);
	JOB.LoadMatrix();
	
	// Load site force
	vector<int>     vind;
	vector<int>     vres;
	vector<string>  vatn;
	vector<Eigen::Vector3d> vfrc;
	ifstream ff(argv[3]);
	if (!ff)
	{
		cerr << "error: file " << argv[3] << " not exists.\n";
		return 0;
	}
	else
	{
		string s;
		while (getline(ff, s))
		{
			if (s.empty() || s.find("REMARK", 0) != string::npos)
				continue;

			int itmp, itmp2;
			string stmp;
			Eigen::Vector3d vtmp(0., 0., 0.);
			istringstream is(s);
			is >> itmp >> itmp2 >> stmp 
			>> vtmp.x() >> vtmp.y() >> vtmp.z();
			if (!is)
			{
				cerr << "error: something bad happen!\n";
				return 0;
			}
			vind.push_back(itmp);
			vres.push_back(itmp2);
			vatn.push_back(stmp);
			vfrc.push_back(vtmp);
		}
		ff.close();
	}

	cout << "REMARK NUMBER OF FORCE = " << vfrc.size() << '\n';

	if (vfrc.size() != natom)
	{
		cerr << "\n\nerror: # of forces != natom\n";
		return 0;
	}

	if (vfrc.size() != JOB._dimension() / 3)
	{
		cerr << "error: vfrc.size != JOB.dimension / 3\n";
		return 0;
	}

	cout << "REMARK OUTPUT LRT DISPLACEMENT VECTOR DX DY DZ vabs(DR)\n";

	for (int i = 0; i < natom; i++)
	{
		Eigen::Vector3d dr(0., 0., 0.);
		
		for (int j = 0; j < natom; j++)
		{
			dr.x() += MATRIX[i * 3    ][j * 3    ] * vfrc[j].x()
				+ MATRIX[i * 3    ][j * 3 + 1] * vfrc[j].y()
				+ MATRIX[i * 3    ][j * 3 + 2] * vfrc[j].z();
			dr.y() += MATRIX[i * 3 + 1][j * 3    ] * vfrc[j].x()
				+ MATRIX[i * 3 + 1][j * 3 + 1] * vfrc[j].y()
				+ MATRIX[i * 3 + 1][j * 3 + 2] * vfrc[j].z();
			dr.z() += MATRIX[i * 3 + 2][j * 3    ] * vfrc[j].x()
				+ MATRIX[i * 3 + 2][j * 3 + 1] * vfrc[j].y()
				+ MATRIX[i * 3 + 2][j * 3 + 2] * vfrc[j].z();	
		}
		// output displacement vector in the force style;

		dr = dr / kBT;

		cout << setw(5) << vind[i]
		<< setw(4) << vres[i]
		<< setw(5) << vatn[i]
		<< setprecision(8) << scientific
		<< setw(16) << dr.x()
		<< setw(16) << dr.y()
		<< setw(16) << dr.z()
		<< setw(16) << dr.norm()
		<< '\n';
	}
}
