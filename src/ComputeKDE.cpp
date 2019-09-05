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
