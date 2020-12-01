#ifndef ___COMPUTEHISTGRAM
#define ___COMPUTEHISTGRAM

#include <vector>

class ComputeHistgram
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
	std::vector<unsigned long int> histgram;
	std::vector<double> prob_hist;

	public:
	ComputeHistgram(std::vector<double>* ptr_dataVector,
			double vmin, 
			double vmax, 
			int    nbin, 
			bool   normalize);
	ComputeHistgram(
			double vmin,
			double vmax,
			int    nbin,
			bool   normalize);

	bool load_wham_data(std::string filename, double center, double consk);
	void do_normalize();
	void calc_histgram();
	void output();

	int    _nsample() { return nsample; }
	double _center()  { return center; }
	double _consk()   { return consk; }
};
#endif
