#ifndef ___COMPUTEPATHCV
#define ___COMPUTEPATHCV

#include <vector>
#include "Atom.hpp"

class ComputePathCV
{
	private:
	unsigned int nrep, ndim;
	std::string filename;
	double progress, distance, lambda;
	std::vector<int>* ptr_indexVector;
	std::vector<Atom>* ptr_atomVector;
	std::vector< std::vector<double> > images;

	public:
	ComputePathCV(std::string filename, std::vector<int>* ptr_indexVector, std::vector<Atom>* ptr_atomVector);

	void calc_pathcv();
	double _progress() { return progress; }
	double _distance() { return distance; }

	private:
	void open_images();
	void calc_lambda();
};
#endif
