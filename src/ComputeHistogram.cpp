#include <iostream>
#include <iomanip>
#include <fstream>
#include <numeric>
#include <cmath>
#include "FileIO.hpp"
#include "common.hpp"
#include "ComputeHistogram.hpp"

using namespace std;

ComputeHistogram::ComputeHistogram(vector<double>* ptr_dataVector, double vmin, double vmax, int nbin, bool normalize)
{	
	this->ptr_dataVector = ptr_dataVector;

	initialize( vmin, vmax, nbin, normalize );
}

ComputeHistogram::ComputeHistogram(double vmin, double vmax, int nbin, bool normalize)
{
	initialize( vmin, vmax, nbin, normalize );
}

void ComputeHistogram::initialize(double vmin, double vmax, int nbin, bool normalize)
{
	this->vmin = vmin;
	this->vmax = vmax;
	this->nbin = nbin;
	this->normalize = normalize;

	this->w = (vmax - vmin) / nbin;

	this->histogram.resize(nbin, 0.);
	this->w_histogram.resize(nbin, 0.);

	this->prob_hist.resize(nbin, 0.);

	this->coordinates.resize(nbin, 0.);

	calc_coordinates();
}

void ComputeHistogram::calc_coordinates()
{
	for (int i = 0; i < nbin; i++)
		coordinates[i] = vmin + w / 2. * (2. * i + 1);
}

bool ComputeHistogram::load_wham_data(string filename, double center, double consk)
{
	this-> center   = center;
	this-> consk    = consk;
	this-> fene_old = 0;
	this-> fene_new = 0;

	FileIO fi(filename);
	if (!fi.load_data(2, &dataVector)) return false;

	this->ptr_dataVector = &dataVector;

	return true;
}

bool ComputeHistogram::load_data(string filename)
{
	FileIO fi(filename);
	if (!fi.load_data(2, &dataVector)) return false;

	this->ptr_dataVector = &dataVector;

	return true;
}

bool ComputeHistogram::load_weight(string filename)
{
	FileIO fi(filename);
	if (!fi.load_data(2, &weightVector)) return false;

	return true;
}

void ComputeHistogram::do_normalize()
{
	double dsum = static_cast<double>( nsample );

	for (int i = 0; i < nbin; i++)
	{
		prob_hist[i] = histogram[i] / (w * dsum);
	}
}

void ComputeHistogram::calc_histogram()
{
	for (auto& data: *ptr_dataVector)
	{
		for (int i = 0; i < nbin; i++)
		{
			//if (vmin + w * i <= data && data < vmin + w * (i + 1))
			if ( abs( data - coordinates[i] ) < w * 0.5)
			{
				histogram[i] += 1;
				break;
			}
		}
	}

	nsample = accumulate(histogram.begin(), histogram.end(), 0);
}

void ComputeHistogram::calc_weighted_histogram()
{
	unsigned int icnt = 0;
	for (auto& data: *ptr_dataVector)
	{
		for (int i = 0; i < nbin; i++)
		{
			//if (vmin + w * i <= data && data < vmin + w * (i + 1))
			if ( abs( data - coordinates[i] ) < w * 0.5)
			{
				w_histogram[i] += weightVector[icnt];
				break;
			}
		}

		++icnt;
	}
}

void ComputeHistogram::output()
{
	if (normalize) do_normalize();

	cout << setprecision(4) << scientific;

	for (int i = 0 ; i < nbin; i++)
	{
		cout 
			<< setw(12) << coordinates[i];
		
		if (normalize)
		{
			cout
			<< setw(16) << prob_hist[i];
		}
		else
		{
			cout
			<< setw(16) << histogram[i];
		}

		cout
			<< '\n';
	}
	cout << "REMARK " << nsample << "/" << ptr_dataVector->size()
		<< " SAMPLES COLLECTED.\n";
}

void ComputeHistogram::output_pmf(double kbT)
{
	for (int i = 0; i < nbin; i++)
		w_histogram[i] = -kbT * log( w_histogram[i] );

	double pivot = 9999.;
	for (int i = 0; i < nbin; i++)
	{
		if ( isfinite( w_histogram[i] ) && w_histogram[i] < pivot )
			pivot = w_histogram[i];
	}

	cout << setprecision(4) << scientific;

	for (int i = 0 ; i < nbin; i++)
	{
		cout
			<< setw(12) << coordinates[i]
			<< setw(16) << w_histogram[i] - pivot
			<< '\n';
	}
}
