#ifndef ___CLASS_COMPUTEMBAR
#define ___CLASS_COMPUTEMBAR

#include <vector>

class Bias
{
	private:
	int ncvs;

	public:
	std::vector< std::vector<double> > data;
	double center;
	double consk;
	double fene_new, fene_old;

	Bias(std::string filename, int ncvs, double center, double consk);

	void load_data(std::string filename);
};

class ComputeMBAR
{
	private:
	int    ncvs;
	double tol;
	double temperature;
	double kbT;
	double beta;
	int istep;
	std::string speriod;
	bool is_periodic;
	double period;
	std::vector<Bias> biases;

	public:
	ComputeMBAR(std::string metafilename, int ncvs, double tol, double temperature, std::string speriod);

	private:
	void load_metafile(std::string metafilename);
	double wrap_delta(double diff);

	public:
	bool check_convergence();
	void mbar_iteration();
	void output_results();
};
#endif
