#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cmath>
#include "common.hpp"
#include "ComputeMBAR.hpp"

using namespace std;

Bias::Bias(string filename, int ncvs, double center, double consk)
{
	this->ncvs     = ncvs;
	this->center   = center;
	this->consk    = consk;
	this->fene_new = 0;
	this->fene_old = 0;

	load_data(filename);
}

void Bias::load_data(string filename)
{
	ifstream fi(filename.c_str());
	if (!fi) die("Could not open file " + filename);

	string s;
	while (getline(fi, s))
	{
		if (s.empty() || is_comment(s)) continue;

		int itmp;
		double dtmp;
		vector<double> vtmp;
		istringstream is(s);
		is >> itmp;
		while (is >> dtmp)
			vtmp.push_back(dtmp);

		if (vtmp.size() != ncvs)
			die("error: inconsistency between ncvs and size of data");

		data.push_back(vtmp);
	}

	cout << "REMARK # of data from " << filename << ": " << data.size() << '\n';
}

ComputeMBAR::ComputeMBAR(string metafilename, int ncvs, double tol, double temperature, string speriod)
{
	this->ncvs        = ncvs;
	this->tol         = tol;
	this->temperature = temperature;
	this->kbT         = temperature * BOLTZMAN;
	this->beta        = 1. / kbT;
	this->istep       = 0;

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

	load_metafile(metafilename);
}

void ComputeMBAR::load_metafile(string metafilename)
{
	ifstream fi(metafilename.c_str());
	if (!fi)
	{
		die("Could not open file " + metafilename);
	}

	string s;
	vector<string> filenames;
	while (getline(fi, s))
	{
		if (s.empty() || is_comment(s))
			continue;

		istringstream is(s);
		string filename;
		double center, consk;
		is >> filename >> center >> consk;

		Bias btmp(filename, ncvs, center, consk);

		biases.push_back(btmp);
	}
}

double ComputeMBAR::wrap_delta(double diff)
{
	if (diff < -period / 2.) diff += period;
	if (diff >  period / 2.) diff -= period;

	return diff;
}

bool ComputeMBAR::check_convergence()
{
	for (auto& b: biases)
		if ( abs(b.fene_new - b.fene_old) > tol )
		{
			double sum = 0;
			for (auto& b: biases)
			{
				double del = b.fene_new - b.fene_old;
				sum += del * del;
			}
			sum = sqrt(sum/biases.size());
			cout << "REMARK rmsd = " << sum << '\n';
			return false;
		}

	cout << setprecision(6) << fixed;
	for (int i = 0; i < biases.size(); i++)
	{
		cout
			//<< setw(16) << coordinates[i]
			<< setw(16) << biases[i].fene_new - biases[0].fene_new
			<< '\n';
	}
	return true;
}

void ComputeMBAR::mbar_iteration()
{
	for (auto& b: biases)
		b.fene_old = b.fene_new;

	for (auto& bi: biases)
	{
		double sum = 0;

		for (auto& bj: biases)
		{
			for (auto& xjns: bj.data)
			{
				double numer = 0;
				for (auto& xjn: xjns)
				{
					double del = xjn - bi.center;
					numer += del * del;
				}
				numer = exp( -beta * 0.5 * bi.consk * numer );

				double denom = 0;
				for (auto& bk: biases)
				{
					double dtmp = 0;
					for (auto& xjn: xjns)
					{
						double del = xjn - bk.center;
						dtmp += del * del;
					}
					dtmp = beta * (bk.fene_old - 0.5 * bk.consk * dtmp);
					dtmp = exp ( dtmp );

					denom += bk.data.size() * dtmp;
				}

				sum += numer / denom;
			}
		}

		bi.fene_new = -kbT * log( sum );
	}

	if (!(++istep % 10))
		cout << "REMARK Iteration # " << istep << '\n';
}

void ComputeMBAR::output_results()
{
}
