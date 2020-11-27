#ifndef ___COMPUTEHISTGRAM
#define ___COMPUTEHISTGRAM

#include <vector>

class ComputeHistgram
{
	private:
	std::vector<double>* ptr_dataVector;
	double vmin, vmax, w;
	int    nbin;
	bool   normalize;
	std::vector<unsigned long int> histgram;

	public:
	ComputeHistgram(std::vector<double>* ptr_dataVector,
			double vmin, 
			double vmax, 
			int    nbin, 
			bool   normalize);

	void calc_histgram();
	void output();
};
#endif
