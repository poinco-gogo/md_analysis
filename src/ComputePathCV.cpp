#include <iostream>
#include <fstream>
#include <cmath>
#include "common.hpp"
#include "ComputePathCV.hpp"

using namespace std;

ComputePathCV::ComputePathCV(string filename, vector<int>* ptr_indexVector, vector<Atom>* ptr_atomVector)
{
	this->filename = filename;

	this->ptr_indexVector = ptr_indexVector;

	this->ptr_atomVector = ptr_atomVector;

	this->progress = 0.;
	this->distance = 0.;

	open_images();

	calc_lambda();
}

void ComputePathCV::open_images()
{
	ifstream fi(filename);
	if (!fi) 
	{
		cerr << "\nerror: Could not open file " << filename << '\n';
		die("terminating...");
	}

	string s;
	while (getline(fi, s))
	{
		int itmp;
		double dtmp;
		istringstream is(s);
		is >> itmp;
		vector<double> vtmp;
		while (is >> dtmp)
			vtmp.push_back(dtmp);

		images.push_back(vtmp);
	}

	this->nrep = images.size();
	cout << "REMARK Number of replicas: " << nrep << '\n';

	this->ndim = images[0].size();
	if (ndim != ptr_indexVector->size() * 3)
	{
		err("Found inconsistent dimension in path.dat and ind.");
		die("exiting...");
	}

	cout << "REMARK Dimension: " << ndim << '\n';
}

void ComputePathCV::calc_lambda()
{
	double sum = 0;
	for (int i = 0; i < nrep - 1; i++)
	{
		double var = 0;
		for (int j = 0; j < ndim; j++)
		{
			double del = images[i + 1][j] - images[i][j];

			var += del * del;
		}

		sum += sqrt(var);
	}

	lambda = sum / (nrep - 1);

	lambda = 2.3 / (lambda * lambda);
}

void ComputePathCV::calc_pathcv()
{
	double numer = 0;
	double denom = 0;

	for (int i = 0; i < nrep; i++)
	{
		double del = 0;
		for (int j = 0; j < ndim / 3; j++)
		{
			Atom& atom =
			ptr_atomVector->at( ptr_indexVector->at(j) - 1 );
			double dx = atom.position.x() - images[i][3 * j    ];
			double dy = atom.position.y() - images[i][3 * j + 1];
			double dz = atom.position.z() - images[i][3 * j + 2];

			del += dx * dx + dy * dy + dz * dz;
		}

		del = exp( -lambda * del );

		numer += (i + 1) * del;

		denom += del;
	}

	progress = numer / denom;

	distance = -1.0 / lambda * log(denom);
}
