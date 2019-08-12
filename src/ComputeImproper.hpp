#ifndef ___CLASS_COMPUTE_IMPROPER
#define ___CLASS_COMPUTE_IMPROPER

#include <vector>
#include "Improper.hpp"

class ComputeImproper
{
	private:

	double sum_vimproper;

	std::vector<Improper>* ptr_improperVector;
	
	public:

	ComputeImproper(){};
	ComputeImproper(std::vector<Improper>& improperVector);

	void reset();

	double compute_force();
};
#endif
