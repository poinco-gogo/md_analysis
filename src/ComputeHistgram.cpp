#include <iostream>
#include <iomanip>
#include <numeric>
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
}

void ComputeHistgram::output()
{
	int    isum = accumulate(histgram.begin(), histgram.end(), 0);

	double dsum = static_cast<double>( isum );

	cout << setprecision(4) << scientific;

	for (int i = 0 ; i < nbin; i++)
	{
		cout 
			<< setw(12) << vmin + w * i
			<< setw(12) << vmin + w * (i + 1);
		
		if (normalize)
		{
			cout
			<< setw(16) << histgram[i] / (w * dsum);
		}
		else
		{
			cout
			<< setw(16) << histgram[i];
		}

		cout
			<< '\n';
	}
	cout << "REMARK " << isum << "/" << ptr_dataVector->size()
		<< " SAMPLES COLLECTED.\n";
}
