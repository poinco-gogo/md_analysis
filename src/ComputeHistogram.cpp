#include <iostream>
#include <iomanip>
#include <fstream>
#include <numeric>
#include "common.hpp"
#include "ComputeHistgram.hpp"

using namespace std;

ComputeHistgram::ComputeHistgram(vector<double>* ptr_dataVector, double vmin, double vmax, int nbin, bool normalize)
{	
	this->ptr_dataVector = ptr_dataVector;
	this->vmin = vmin;
	this->vmax = vmax;
	this->nbin = nbin;
	this->normalize = normalize;

	this->w = (vmax - vmin) / nbin;

	this->histgram.resize(nbin, 0.);

	this->prob_hist.resize(nbin, 0.);
}

ComputeHistgram::ComputeHistgram(double vmin, double vmax, int nbin, bool normalize)
{
	this->vmin = vmin;
	this->vmax = vmax;
	this->nbin = nbin;
	this->normalize = normalize;

	this->w = (vmax - vmin) / nbin;

	this->histgram.resize(nbin, 0.);

	this->prob_hist.resize(nbin, 0.);
}

bool ComputeHistgram::load_wham_data(string filename, double center, double consk)
{
	this-> center   = center;
	this-> consk    = consk;
	this-> fene_old = 0;
	this-> fene_new = 0;

	ifstream fi(filename.c_str());

	if (!fi)
	{
		cerr << "\nerror: file \"" << filename << "\" not exists.\n\n";
		return false;
	}

	string s;
	while (getline(fi, s))
	{
		if (s.empty() || is_comment(s))
			continue;

		istringstream is(s);
		double val;
		// note that only 2nd column is used
		for (int i = 0; i < 2; i++)
			is >> val;

		dataVector.push_back( val );
	}

	this->ptr_dataVector = &dataVector;

	return true;
}

void ComputeHistgram::do_normalize()
{
	double dsum = static_cast<double>( nsample );

	for (int i = 0; i < nbin; i++)
	{
		prob_hist[i] = histgram[i] / (w * dsum);
	}
}

void ComputeHistgram::calc_histgram()
{
	for (auto& data: *ptr_dataVector)
	{
		for (int i = 0; i < nbin; i++)
		{
			if (vmin + w * i <= data && data < vmin + w * (i + 1))
			{
				histgram[i] += 1;
				break;
			}
		}
	}

	nsample = accumulate(histgram.begin(), histgram.end(), 0);
}

void ComputeHistgram::output()
{
	if (normalize) do_normalize();

	cout << setprecision(4) << scientific;

	for (int i = 0 ; i < nbin; i++)
	{
		cout 
			<< setw(12) << vmin + w * i
			<< setw(12) << vmin + w * (i + 1);
		
		if (normalize)
		{
			cout
			<< setw(16) << prob_hist[i];
		}
		else
		{
			cout
			<< setw(16) << histgram[i];
		}

		cout
			<< '\n';
	}
	cout << "REMARK " << nsample << "/" << ptr_dataVector->size()
		<< " SAMPLES COLLECTED.\n";
}
