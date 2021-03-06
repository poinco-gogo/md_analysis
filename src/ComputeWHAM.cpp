#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cmath>
#include "common.hpp"
#include "ComputeWHAM.hpp"

using namespace std;

ComputeWHAM::ComputeWHAM(string metafilename, double vmin, double vmax, int nbin, double tol, double temperature, string speriod)
{
	this->vmin        = vmin;
	this->vmax        = vmax;
	this->nbin        = nbin;
	this->tol         = tol;
	this->temperature = temperature;
	this->kbT         = temperature * BOLTZMAN;
	this->beta        = 1. / kbT;
	this->istep       = 0;
	this->dz          = (vmax - vmin) / nbin;

	this->is_periodic = true;
	if (speriod == "P")
	{
		period = N_TWOPI * RAD2DEG;
	}
	else if (speriod == "Ppi")
	{
		period = N_TWOPI;
	}
	else
	{
		this->is_periodic = false;
	}

	if (!load_metafile(metafilename))
		return;

	for (int i = 0; i < nbin; i++)
	{
		coordinates.push_back(vmin + dz / 2. * (2 * i + 1));

		prob_global.push_back(0.);
	}

	// turns histograms into probability.
	for (auto& hi: histograms)
	{
		hi.do_normalize();
	}
}

bool ComputeWHAM::load_metafile(string metafilename)
{
	ifstream fi(metafilename.c_str());
	if (!fi)
	{
		cerr << "\nerror: file \"" << metafilename
			<< "\" not exists.\n\n";
		return false;
	}

	string s;
	while (getline(fi, s))
	{
		if (s.empty() || is_comment(s))
			continue;

		istringstream is(s);
		string filename;
		double center, consk;
		is >> filename >> center >> consk;

		ComputeHistogram vhist(vmin, vmax, nbin, false);

		vhist.load_wham_data(filename, center, consk);

		vhist.calc_histogram();

		histograms.push_back(vhist);
	}

	return true;
}

double ComputeWHAM::wrap_delta(double diff)
{
	if (diff < -period / 2.) diff += period;
	if (diff >  period / 2.) diff -= period;

	return diff;
}

bool ComputeWHAM::check_convergence()
{
	for (auto& hi: histograms)
		if ( abs(hi.fene_new - hi.fene_old) > tol )
			return false;
	return true;
}

void ComputeWHAM::wham_iteration()
{
	for (auto& hi: histograms)
		hi.fene_old = hi.fene_new;

	// calc global probability density function
	for (int k = 0; k < coordinates.size(); k++)
	{
		double nP = 0.;
		double denom = 0.;
		for (auto& hi: histograms)
		{
			nP += hi._nsample() * hi.prob_hist[k];
			double dist = coordinates[k] - hi._center();
			if (is_periodic) dist = wrap_delta(dist);
			denom += hi._nsample() * 
			exp(
				beta * (hi.fene_old
				- 0.5 * hi._consk() *  dist * dist)
			);
		}
		prob_global[k] = nP / denom;
	}

	// calc free energy associated with adding bias
	for (auto& hi: histograms)
	{
		hi.fene_new = 0.;

		for (int l = 0; l < coordinates.size(); l++)
		{
			double d = coordinates[l] - hi._center();
			if (is_periodic) d = wrap_delta(d);
			hi.fene_new +=
			prob_global[l] * dz *
			exp(
				-beta * 0.5 * hi._consk()
				* d * d
			);
		}
	}

	for (auto& hi: histograms)
		hi.fene_new = -kbT * log( hi.fene_new );

	istep++;
}

void ComputeWHAM::output_results()
{
	vector<double> pmf(nbin, 0.);

	for (int i = 0; i < nbin; i++)
	{
		pmf[i] = -kbT * log( prob_global[i] );
	}

	double pivot = 9999.;
	for (int i = 0; i < nbin; i++)
	{
		if ( isfinite( pmf[i] ) && pmf[i] < pivot )
			pivot = pmf[i];
	}

	cout << setprecision(6) << fixed;
	for (int i = 0; i < nbin; i++)
	{
		cout
			<< setw(16) << coordinates[i]
			<< setw(16) << pmf[i] - pivot
			<< '\n';
	}
}
