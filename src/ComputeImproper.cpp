#include <iostream>
#include <iomanip>
#include "ComputeImproper.hpp"
using namespace std;

ComputeImproper::ComputeImproper(vector<Improper>& improperVector)
{
	this -> ptr_improperVector = &improperVector;

	reset();
}

void ComputeImproper::reset()
{
	this->sum_vimproper = 0;
}

double ComputeImproper::compute_force()
{
	reset();

	for (Improper& imp: *ptr_improperVector)
	{
		imp.calc_force();

		imp.ptr_atom1->force += imp.force1;
		imp.ptr_atom2->force += imp.force2;
		imp.ptr_atom3->force += imp.force3;
		imp.ptr_atom4->force += imp.force4;

		sum_vimproper += imp.energy;
	}

	return sum_vimproper;
}
