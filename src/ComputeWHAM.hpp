#ifndef ___CLASS_COMPUTEWHAM
#define ___CLASS_COMPUTEWHAM

#include <vector>
#include "ComputeHistogram.hpp"

class ComputeWHAM
{
	private:
	double vmin, vmax;
	int    nbin;
	double dz;
	double tol;
	double temperature;
	double kbT;
	double beta;
	int istep;
	std::string speriod;
	bool is_periodic;
	double period;

	public:
	std::vector<ComputeHistogram> histograms;
	std::vector<double> coordinates;
	std::vector<double> prob_global;

	// constructor
	ComputeWHAM(std::string metafilename, double vmin, double vmax, int nbin, double tol, double temperature, std::string speriod);

	private:
	bool load_metafile(std::string filename);
	double wrap_delta(double diff);

	public:
	bool check_convergence();
	void wham_iteration();
	void output_results();
};
#endif
