#include <iostream>
#include <cmath>
#include "ComputeKDE.hpp"
#include "common.hpp"
using namespace std;

ComputeKDE::ComputeKDE(vector<double>* ptr_dataVector1, double band_width)
{
	this->ptr_dataVector1    = ptr_dataVector1;

	this->band_width         = band_width;
}

ComputeKDE::ComputeKDE(vector<double>* ptr_dataVector1,	vector<double>* ptr_dataVector2, double band_width)
{
	this->ptr_dataVector1    = ptr_dataVector1;
	this->ptr_dataVector2    = ptr_dataVector2;

	this->band_width         = band_width;
}

ComputeKDE::ComputeKDE(vector<double>* ptr_dataVector1,	vector<double>* ptr_weightVector, vector<double>* ptr_dataVector2, double band_width, double cutoff)
{
	this->ptr_dataVector1    = ptr_dataVector1;
	this->ptr_weightVector   = ptr_weightVector;
	this->ptr_dataVector2    = ptr_dataVector2;

	this->band_width         = band_width;
	this->cutoff             = cutoff;
}


ComputeKDE::ComputeKDE(vector<double>* ptr_dataVector1,	vector<double>* ptr_dataVector2, vector<double>* ptr_weightVector, vector<double>* ptr_distanceVector, double min1, double max1, int nbin1, double band_width1, double min2, double max2, int nbin2, double band_width2, double cutoff)
{
	this->ptr_dataVector1    = ptr_dataVector1;
	this->ptr_dataVector2    = ptr_dataVector2;

	this->ptr_weightVector   = ptr_weightVector;

	this->ptr_distanceVector = ptr_distanceVector;

	this->min1               = min1;
	this->max1               = max1;
	this->nbin1              = nbin1;
	this->band_width1        = band_width1;

	this->min2               = min2;
	this->max2               = max2;
	this->nbin2              = nbin2;
	this->band_width2        = band_width2;

	this->cutoff             = cutoff;
}

double ComputeKDE::estimate_gauss(double x)
{
	double sum = 0;
	for (auto& data: *ptr_dataVector1)
	{
		double xx = ( x - data ) / band_width;
		xx *= xx;

		sum += exp( - xx * 0.5 );
	}
	sum /= sqrt( 2.0 * N_PI ) * ptr_dataVector1->size() * band_width;

	return sum;
}

double ComputeKDE::estimate_gauss(double x, double y)
{
	double sum = 0;
	for (int i = 0; i < ptr_dataVector1->size(); i++)
	{
		double xdata = ptr_dataVector1->at(i);
		double ydata = ptr_dataVector2->at(i);

		double xx = ( x - xdata ) / band_width;
		double yy = ( y - ydata ) / band_width;
		xx *= xx;
		yy *= yy;

		sum += exp( - ( xx + yy ) * 0.5 );
	}
	sum /= 2.0 * N_PI * ptr_dataVector1->size() * band_width * band_width;

	return sum;
}

double ComputeKDE::estimate_gauss_weight(double x)
{
	double sum = 0;
	int icnt = 0;
	for (auto& data: *ptr_dataVector1)
	{
		if (ptr_dataVector2->at(icnt++) > cutoff)
			continue;

		double xx = ( x - data ) / band_width;
		xx *= xx;

		sum += ptr_dataVector2->at(icnt++) * exp( - xx * 0.5 );
	}
	sum /= sqrt( 2.0 * N_PI ) * ptr_dataVector1->size() * band_width;

	return sum;
}

double ComputeKDE::estimate_gauss_weight(double x, double y)
{
	double sum = 0;
	for (int i = 0; i < ptr_dataVector1->size(); i++)
	{
		if (ptr_distanceVector->at(i) > cutoff)
			continue;

		double xdata = ptr_dataVector1->at(i);
		double ydata = ptr_dataVector2->at(i);

		double xx = ( x - xdata ) / band_width1;
		double yy = ( y - ydata ) / band_width2;
		xx *= xx;
		yy *= yy;

		sum += ptr_weightVector->at(i) * exp( - ( xx + yy ) * 0.5 );
	}
	sum /= 2.0 * N_PI * ptr_dataVector1->size() * band_width1 * band_width2;

	return sum;
}
