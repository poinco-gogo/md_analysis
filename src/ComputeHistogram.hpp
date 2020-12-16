#ifndef ___COMPUTEHISTOGRAM
#define ___COMPUTEHISTOGRAM

#include <vector>

class ComputeHistogram
{
	private:
	std::vector<double>  dataVector, weightVector;
	std::vector<double>* ptr_dataVector;
	std::vector<double> coordinates;
	double vmin, vmax, w;
	int    nbin, nsample;
	bool   normalize;

	// for WHAM
	double center, consk;

	public:
	double fene_old, fene_new;
	std::vector<unsigned long int> histogram;
	std::vector<double> w_histogram;
	std::vector<double> prob_hist;

	public:
	ComputeHistogram(std::vector<double>* ptr_dataVector,
			double vmin, 
			double vmax, 
			int    nbin, 
			bool   normalize);
	ComputeHistogram(
			double vmin,
			double vmax,
			int    nbin,
			bool   normalize);

	bool load_wham_data(std::string filename, double center, double consk);
	bool load_data(std::string filename);
	bool load_weight(std::string filename);
	void do_normalize();
	void calc_histogram();
	void calc_weighted_histogram();
	void output();
	void output_pmf(double kbT);

	int    _nsample() { return nsample; }
	double _center()  { return center; }
	double _consk()   { return consk; }

	private:
	void initialize(
			double vmin,
			double vmax,
			int    nbin,
			bool   normalize);
	void calc_coordinates();
};
#endif
