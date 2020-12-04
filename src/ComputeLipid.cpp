#include <iostream>
#include <iomanip>
#include "common.hpp"
#include "ComputeLipid.hpp"

using namespace std;

ComputeLipid::ComputeLipid(double vmin, double vmax, int nbin, unsigned int mode, vector<Atom>* ptr_atomVector)
{
	this->vmin           = vmin;
	this->vmax           = vmax;
	this->nbin           = nbin;
	this->w              = (vmax - vmin) / nbin;
	this->mode           = mode;
	this->ptr_atomVector = ptr_atomVector;

	initialize();

	output_info();
}

void ComputeLipid::initialize()
{
	histogram.resize(nbin, 0.);
	coordinates.resize(nbin, 0.);

	for (int i = 0; i < nbin; i++)
		coordinates[i] = vmin + w / 2. * (2. * i + 1);

}

void ComputeLipid::output_info()
{
	switch (mode) {
		case 1:
			cout << "REMARK Compute number density profile.\n";
			break;
		case 2:
			cout << "REMARK Compute mass density profile.\n";
			break;
		default:
			die("Unknown mode is specified.");
	}
}

void ComputeLipid::set_selection()
{
	for (auto& atom: *ptr_atomVector)
	{
		if (atom.PDBAtomName == "P")
			tgt_list.push_back(atom.PSFIndex);
	}

	cout << "REMARK Found " << tgt_list.size() << " P atom(s).\n";
}

void ComputeLipid::calc_histogram(double boxx, double boxy)
{
	vector<double> vtmp(nbin, 0.);

	for (auto& id: tgt_list)
	{
		Atom& atom = ptr_atomVector->at(id - 1);

		for (int i = 0; i < nbin; i++)
		{
			double z = coordinates[i];

			if ( abs( atom.position.z() - z ) < 0.5 * w )
			{
				vtmp[i] += 1.0;
				break;
			}
		}
	}

	// accumulate density
	for (int i = 0; i < nbin; i++)
	{
		histogram[i] += vtmp[i] / (boxx * boxy * w);
	}
}

void ComputeLipid::calc_density_profile(int nsteps)
{
	double weight = 0;

	switch (mode) {
		case 1: weight =    1.0; break;
		case 2: weight = 30.974; break;
	}

	for (int i = 0; i < nbin; i++)
		histogram[i] = histogram[i] * weight / nsteps;
}

void ComputeLipid::output_density_profile()
{
	cout << setprecision(6) << scientific;

	for (int i = 0; i < nbin; i++)
	{
		cout
			<< setw(16) << coordinates[i]
			<< setw(16) << histogram[i]
			<< '\n';
	}
}
