#ifndef ___CLASS_COMPUTEMBAR
#define ___CLASS_COMPUTEMBAR

#include <vector>

class Bias
{
	private:
	int ndim;

	public:
	std::vector< std::vector<double> > data;
	double center;
	double consk;
	double fene_new, fene_old;
	std::vector<double> Wni;

	Bias(std::string filename, int ndim, double center, double consk);

	void load_data(std::string filename);
};

class ComputeMBAR
{
	private:
	int    ndim, nbin;
	double vmin, vmax, dz, tol;
	std::vector<double> coordinates;
	double temperature;
	double kbT;
	double beta;
	int istep;
	std::string speriod;
	bool is_periodic;
	double period;
	std::vector<Bias> biases;
	double ci;

	public:
	ComputeMBAR(std::string metafilename, int ndim, double vmin, double vmax, double nbin, double tol, double temperature, std::string speriod);

	private:
	void load_metafile(std::string metafilename);
	double wrap_delta(double diff);
	void output_weights();
	void output_pmf();
	void calc_coordinates();

	public:
	bool check_convergence();
	void mbar_iteration();
	void calc_weights();
	void output_results();
};
#endif
