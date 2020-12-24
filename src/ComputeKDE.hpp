#ifndef ___CLASS_COMPUTEKDE
#define ___CLASS_COMPUTEKDE

#include <vector>

class ComputeKDE
{
	private:
	std::vector<double>* ptr_dataVector1;
	std::vector<double>* ptr_dataVector2;
	std::vector<double>* ptr_weightVector;
	double band_width;
	double min1, max1, band_width1;
	double min2, max2, band_width2;
	int nbin1, nbin2;

	public:
	ComputeKDE(std::vector<double>* ptr_dataVector1, double band_width);
	ComputeKDE(
			std::vector<double>* ptr_dataVector1,
			std::vector<double>* ptr_dataVector2,
			double band_width
		);
	ComputeKDE(
			std::vector<double>* ptr_dataVector1,
			std::vector<double>* ptr_dataVector2,
			std::vector<double>* ptr_weightVector,
			double min1, double max1, int nbin1, double band_width1,
			double min2, double max2, int nbin2, double band_width2
		);

	double estimate_gauss(double x);
	double estimate_gauss(double x, double y);
	double estimate_gauss_weight(double x);
	double estimate_gauss_weight(double x, double y);
};
#endif
