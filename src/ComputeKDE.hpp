#ifndef ___CLASS_COMPUTEKDE
#define ___CLASS_COMPUTEKDE

#include <vector>

class ComputeKDE
{
	private:
	std::vector<double>* ptr_dataVector;
	double band_width;

	public:
	ComputeKDE(std::vector<double>* ptr_dataVector, double band_width);

	double estimate_gauss(double x);
};
#endif
