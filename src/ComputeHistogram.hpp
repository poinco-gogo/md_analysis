#ifndef ___COMPUTEHISTOGRAM
#define ___COMPUTEHISTOGRAM

#include <vector>

class ComputeHistogram
{
	private:
	std::vector<double>  dataVector;
	std::vector<double>* ptr_dataVector;
	double vmin, vmax, w;
	int    nbin, nsample;
	bool   normalize;

	// for WHAM
	double center, consk;

	public:
	double fene_old, fene_new;
	std::vector<unsigned long int> histogram;
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
	void do_normalize();
	void calc_histogram();
	void output();

	int    _nsample() { return nsample; }
	double _center()  { return center; }
	double _consk()   { return consk; }
};
#endif
