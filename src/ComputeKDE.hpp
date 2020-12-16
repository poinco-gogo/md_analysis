#ifndef ___CLASS_COMPUTEKDE
#define ___CLASS_COMPUTEKDE

#include <vector>

class ComputeKDE
{
	private:
	std::vector<double>* ptr_dataVector1;
	std::vector<double>* ptr_dataVector2;
	double band_width;

	public:
	ComputeKDE(std::vector<double>* ptr_dataVector1, double band_width);
	ComputeKDE(
			std::vector<double>* ptr_dataVector1,
			std::vector<double>* ptr_dataVector2,
			double band_width
		);

	double estimate_gauss(double x);
	double estimate_gauss(double x, double y);
	double estimate_gauss_weight(double x);
};
#endif
