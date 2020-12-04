#ifndef ___COMPUTE_LIPID
#define ___COMPUTE_LIPID

#include <vector>
#include "Atom.hpp"

class ComputeLipid
{
	private:
	double vmin, vmax, w;
	int nbin;
	unsigned int mode;
	std::vector<Atom>* ptr_atomVector;
	std::vector<double> coordinates;
	std::vector<int> tgt_list;

	public:
	std::vector<double> histogram;

	public:
	ComputeLipid(double vmin, double vmax, int nbin, unsigned int mode, std::vector<Atom>* ptr_atomVector);

	void set_selection();
	void calc_histogram(double boxx, double boxy);
	void calc_density_profile(int nsteps);
	void output_density_profile();

	private:
	void initialize();
	void output_info();
};
#endif
