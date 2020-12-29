#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cmath>
#include <Eigen/Core>
#include <Eigen/LU>
#include "common.hpp"
#include "ComputeMBAR.hpp"

using namespace std;

Bias::Bias(string filename, unsigned int ndim, double center, double consk)
{
	this->ndim     = ndim;
	this->center   = center;
	this->consk    = consk;
	this->fene_new = 0;
	this->fene_old = 0;
	this->ci       = 0;

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

		if (vtmp.size() != ndim)
			die("error: inconsistency between ndim and size of data");

		data.push_back(vtmp);

		Wna.push_back(0.0);
	}

	cout << "REMARK # of data from " << filename << ": " << data.size() << '\n';
}

ComputeMBAR::ComputeMBAR(string metafilename, unsigned int ndim, double vmin, double vmax, unsigned int nbin, double tol, double temperature, string ofilename, unsigned int nself,  string speriod)
{
	this->ndim         = ndim;
	this->vmin         = vmin;
	this->vmax         = vmax;
	this->nbin         = nbin;
	this->dz           = (vmax - vmin) / nbin;
	this->tol          = tol;
	this->temperature  = temperature;
	this->kbT          = temperature * BOLTZMAN;
	this->beta         = 1. / kbT;
	this->ofilename    = ofilename;
	this->istep        = 0;
	this->ndata        = 0;
	this->nself        = nself;
	this->metafilename = metafilename;

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

	initialize();
}

void ComputeMBAR::initialize()
{
	load_metafile();

	for (auto& b: biases)
		b.Fni.resize(ndata);

	Wni.resize(ndata, biases.size());
}

void ComputeMBAR::load_metafile()
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

		Bias btmp(filename, ndim, center, consk);

		ndata += btmp.data.size();

		biases.push_back(btmp);
	}

	cout << "REMARK total data size: " << ndata << '\n';
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
			cout << setprecision(6) << scientific;
			cout << "REMARK # iteration:" << istep << ", rmsd = " << sum << '\n';
			return false;
		}

	return true;
}

void ComputeMBAR::mbar_iteration()
{
	for (auto& b: biases)
		b.fene_old = b.fene_new;

	calc_Fni_and_ci();

	if (++istep <= nself)
	{
		if (istep == 1)
			cout <<
			"REMARK ===========================================\n"
			"REMARK Self-consistent iterations...\n"
			"REMARK ===========================================\n";

		// attempt self-consistent iteration for at least nself cycles

		mbar_self_consistent();
	}
	else
	{
		// try Newton-Raphson method

		if (istep == nself + 1)
			cout <<
			"REMARK ===========================================\n"
			"REMARK Switch to Newton-Raphson minimizer...\n"
			"REMARK ===========================================\n";

		mbar_newton_raphson();
	}
}

void ComputeMBAR::calc_Fni_and_ci()
{
	for (auto& bi: biases)
	{
		bi.ci = 0;

		int icnt = 0;
		for (auto& bj: biases)
		{
			for (auto& xjns: bj.data)
			{
				double numer = 0;
				for (auto& xjn: xjns)
				{
					double del = xjn - bi.center;
					if (is_periodic) del = wrap_delta(del);
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
						if (is_periodic)
							del = wrap_delta(del);
						dtmp += del * del;
					}
					dtmp = beta * (bk.fene_old - 0.5 * bk.consk * dtmp);
					dtmp = exp ( dtmp );

					denom += bk.data.size() * dtmp;
				}

				bi.Fni[icnt] = numer / denom;

				bi.ci += bi.Fni[icnt++];
			}
		}
	}
}

void ComputeMBAR::mbar_self_consistent()
{
	for (auto& b: biases)
		b.fene_new = -kbT * log( b.ci );

	// constrain f0 = 0
	double ftmp = biases[0].fene_new;
	for (auto& b: biases)
		b.fene_new -= ftmp;
}

void ComputeMBAR::mbar_newton_raphson()
{
	calc_weight_matrix();

	int l = biases.size() - 1;
	Eigen::VectorXd g(l);
	Eigen::MatrixXd h(l, l);

	for (int i = 1; i < biases.size(); i++)
	{
		Bias& bi = biases[i];

		double sum = 0;
		for (int n = 0; n < ndata; n++)
			sum += Wni(n,i);

		g(i-1) = bi.data.size() * ( 1.0 - sum );

		sum = 0;
		for (int n = 0; n < ndata; n++)
		{
			sum +=
			bi.data.size() * Wni(n,i) *
			(1.0 - bi.data.size() * Wni(n,i));
		}
		h(i-1, i-1) = sum;

		for (int j = i + 1; j < biases.size(); j++)
		{
			Bias& bj = biases[j];

			sum = 0;
			for (int n = 0; n < ndata; n++)
			{
				sum -=
				bi.data.size() * Wni(n,i) *
				bj.data.size() * Wni(n,j);
			}

			h(i-1, j-1) = sum;
			h(j-1, i-1) = h(i-1, j-1);
		}
	}

	h = h.inverse();

	Eigen::VectorXd obj = h * g;

	double r = istep == nself + 1 ? 0.1 : 1.0;

	for (int i = 1; i < biases.size(); i++)
	{
		Bias& b = biases[i];

		b.fene_new = b.fene_old + r * obj(i - 1);
	}
}

void ComputeMBAR::calc_weight_matrix()
{
	for (int i = 0; i < biases.size(); i++)
	{
		Bias& bi = biases[i];

		for (int n = 0; n < ndata; n++)
		{
			Wni(n, i)
			= bi.Fni[n] * exp( beta * bi.fene_old );
		}
	}
}

void ComputeMBAR::calc_unbiasing_weights()
{
	double ca = 0.;

	for (auto& b: biases)
	{
		int icnt = 0;
		for (auto& xjns: b.data)
		{
			double numer = 1.0;

			double denom = 0.0;
			for (auto& bk: biases)
			{
				double dtmp = 0.0;
				for (auto& xjn: xjns)
				{
					double del = xjn - bk.center;
					if (is_periodic) del = wrap_delta(del);
					dtmp += del * del;
				}
				dtmp = beta * (bk.fene_new - 0.5 * bk.consk * dtmp);
				dtmp = exp( dtmp );

				denom += bk.data.size() * dtmp;
			}

			b.Wna[icnt] = numer / denom;

			ca += b.Wna[icnt++];
		}
	}

	for (auto& b: biases)
		for (auto& w: b.Wna)
			w /= ca;
}

void ComputeMBAR::output_results()
{
	output_fene();

	output_unbiasing_weights();

	output_pmf();
}

void ComputeMBAR::output_fene()
{
	ofstream fo( ofilename + ".fene" );
	fo << setprecision(12) << fixed;
	fo << "REMARK Free energy of the biased systems (kcal/mol)\n";
	for (int i = 0; i < biases.size(); i++)
	{
		fo
			<< setw(12) << i + 1
			<< setw(20) << biases[i].fene_new - biases[0].fene_new
			<< '\n';
	}
}

void ComputeMBAR::output_unbiasing_weights()
{
	for (int i = 0; i < biases.size(); i++)
	{
		ostringstream os;
		os << ofilename << i + 1 << ".weight";
		ofstream fo(os.str().c_str());
		fo << setprecision(12) << scientific;

		int icnt = 0;
		for (auto& w: biases[i].Wna)
		{
			fo
			<< setw(12) << ++icnt
			<< setw(20) << w
			<< '\n';
		}
	}
}

void ComputeMBAR::output_pmf()
{
	calc_bincenters();

	vector<double> histogram(nbin, 0.);

	for (auto& b: biases)
	{
		int icnt = 0;
		for (auto& xjns: b.data)
		{
			for (auto& xjn: xjns)
			{
				for (int i = 0; i < nbin; i++)
				{
					if ( abs( xjn - bincenters[i] ) < dz * 0.5 )
						histogram[i] += b.Wna[icnt];
				}
			}

			++icnt;
		}
	}

	double pivot = 9999.;
	for (auto& d: histogram)
	{
		d = -kbT * log(d);

		if ( isfinite(d) && d < pivot) pivot = d;
	}

	for (auto& d: histogram)
		d -= pivot;

	ofstream fo( ofilename + ".pmf" );
	fo << setprecision(6) << fixed;
	for (int i = 0; i < nbin; i++)
	{
		fo
			<< setw(16) << bincenters[i]
			<< setw(16) << histogram[i]
			<< '\n';
	}
}

void ComputeMBAR::calc_bincenters()
{
	bincenters.resize(nbin);

	for (int i = 0; i < nbin; i++)
		bincenters[i] = vmin + dz / 2. * (2. * i + 1);
}
