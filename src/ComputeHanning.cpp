#include <iostream>
#include <cmath>
#include <vector>
#include "ComputeHanning.hpp"
using namespace std;

ComputeHanning::ComputeHanning(vector<double>* ptr_HannVector)
{
	this -> ptr_HannVector = ptr_HannVector;
	this -> nsample = ptr_HannVector->size();
	this -> indenom = 0.;
	generate_Hanning();
}

void ComputeHanning::generate_Hanning()
{
	double inN = 1. / (nsample - 1);
	for (int k = 0; k < nsample; k++)
	{
		ptr_HannVector->at(k) = 0.5 - 0.5 * cos(2. * N_PI * k * inN);
		indenom += ptr_HannVector->at(k);
	}
	indenom = 1.0 / indenom;
}

double ComputeHanning::get_Hanning_ave(vector<double>* ptr_DataVector)
{
	double average = 0.;
	vector<double> w_i(ptr_DataVector->size(), 0.);
	
	for (int i = 0; i < nsample; i++)
	{
		average += ptr_HannVector->at(i) * ptr_DataVector->at(i);
	}
	return average * indenom;
}
