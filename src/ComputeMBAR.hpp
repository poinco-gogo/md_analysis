#ifndef ___CLASS_COMPUTEMBAR
#define ___CLASS_COMPUTEMBAR

#include <vector>

class Bias
{
	private:
	unsigned int ndim;

	public:
	std::vector< std::vector<double> > data;
	double center;
	double consk;
	double fene_new, fene_old;
	std::vector<double> Wni;

	Bias(std::string filename, unsigned int ndim, double center, double consk);

	void load_data(std::string filename);
};

class ComputeMBAR
{
	private:
	unsigned int ndim, nbin;
	double vmin, vmax, dz, tol;
	std::vector<double> bincenters;
	double temperature;
	double kbT;
	double beta;
	std::string ofilename;
	int istep;
	std::string speriod;
	bool is_periodic;
	double period;
	std::vector<Bias> biases;
	double ci;

	public:
	ComputeMBAR(std::string metafilename, unsigned int ndim, double vmin, double vmax, unsigned int nbin, double tol, double temperature, std::string ofilename, std::string speriod);

	private:
	void load_metafile(std::string metafilename);
	double wrap_delta(double diff);
	void output_fene();
	void output_unbiasing_weights();
	void output_pmf();
	void calc_bincenters();

	public:
	bool check_convergence();
	void mbar_iteration();
	void calc_unbiasing_weights();
	void output_results();
};
#endif
