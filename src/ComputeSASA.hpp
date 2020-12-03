#ifndef ___CLASS_COMPUTESASA
#define ___CLASS_COMPUTESASA

#include <vector>
#include "Atom.hpp"

class ComputeSASA
{
	private:
	unsigned int npoints, mode;
	std::vector<Atom>* ptr_atomVector;
	std::vector<Eigen::Vector3d> points;

	public:
	ComputeSASA(unsigned int npoints, unsigned int mode, std::vector<Atom>* ptr_atomVector);
	double calc_sasa(unsigned int resid);

	private:
	void assign_radius();
	void generate_points();
	void generate_points_gaussian();
	void generate_points_fibonacci();
};
#endif
