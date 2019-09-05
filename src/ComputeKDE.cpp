#include <iostream>
#include <cmath>
#include "ComputeKDE.hpp"
#include "common.hpp"
using namespace std;

ComputeKDE::ComputeKDE(vector<double>* ptr_dataVector, double band_width)
{
	this->ptr_dataVector     = ptr_dataVector;

	this->band_width         = band_width;
}

double ComputeKDE::estimate_gauss(double x)
{
	double sum = 0;
	for (auto& data: *ptr_dataVector)
	{
		double xx = ( x - data ) / band_width;
		xx *= xx;

		sum += exp( - xx * 0.5 );
	}
	sum /= sqrt( 2.0 * N_PI ) * ptr_dataVector->size() * band_width;

	return sum;
}
