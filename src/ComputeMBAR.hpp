#ifndef ___CLASS_COMPUTEMBAR
#define ___CLASS_COMPUTEMBAR

#include <vector>
#include <Eigen/Core>

class Bias
{
	private:
	unsigned int ndim;

	public:
	std::vector< std::vector<double> > data;
	Eigen::VectorXd center;
	double consk;
	double fene_new, fene_old, ci;
	std::vector<double> Wna, Fni, qni;

	Bias(std::string filename, unsigned int ndim, double center, double consk);
	Bias(std::string filename, unsigned int ndim, double consk, int icnt, Eigen::MatrixXd& Rdi);

	void load_data(std::string filename);
};

class ComputeMBAR
{
	private:
	unsigned int ndim, nbin, ndata, nself, nbias;
	double vmin, vmax, dz, tol;
	std::vector<double> bincenters;
	double temperature;
	double kbT;
	double beta;
	std::string metafilename;
	std::string ofilename;
	int istep;
	std::string speriod;
	bool is_periodic;
	double period;
	std::vector<Bias> biases;
	Eigen::MatrixXd Wni;
	Eigen::MatrixXd Rdi;

	public:
	ComputeMBAR(std::string metafilename, unsigned int ndim, double vmin, double vmax, unsigned int nbin, double tol, double temperature, std::string ofilename, unsigned int nself, std::string speriod);
	ComputeMBAR(std::string metafilename, unsigned int ndim, unsigned int nbias, double tol, double temperature, std::string ofilename, unsigned int nself, std::string speriod);

	private:
	void initialize();
	void load_metafile();
	void load_metafile_posi();
	void load_references(std::string filename);
	void calc_qni();
	double wrap_delta(double diff);
	void calc_Fni_and_ci();
	void mbar_self_consistent();
	void mbar_newton_raphson();
	void calc_weight_matrix();
	void output_fene();
	void output_unbiasing_weights();
	void output_pmf();
	void calc_bincenters();
	void output_rmsd();
	double calc_rmsd();

	public:
	bool check_convergence();
	void mbar_iteration();
	void calc_unbiasing_weights();
	void output_results();
};
#endif
