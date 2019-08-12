#ifndef ___CLASS_COMPUTESPLINE
#define ___CLASS_COMPUTESPLINE
#include <string>
#include <vector>
class ComputeSpline
{
	private:
	std::string filename;
	std::vector<double> x;
	std::vector<double> y;
	std::vector<double> b;
	std::vector<double> c;
	std::vector<double> d;

	public:
	
	ComputeSpline() {};
	ComputeSpline(std::vector<double>& x, std::vector<double>& y);
	ComputeSpline(std::string filename);
	ComputeSpline(int numdata, double x_origin, double stepSize, double* data_value);
	
	bool load_data();

	double calc_cubic_spline(double s);
	double calc_gradient(double s);

	void calc_1d_coefficients();
	void calc_2d_coefficients(double* f, double* dfd1, double* dfd2, double* dfd1d2, double stepSize, double (*w)[4]);
};
#endif
