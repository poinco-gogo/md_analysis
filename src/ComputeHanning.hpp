#ifndef ___CLASS_COMPUTEHANNING
#define ___CLASS_COMPUTEHANNING
#include <vector>
#include "common.hpp"

class ComputeHanning
{
	private:
	int nsample;
	double indenom;
	std::vector<double>* ptr_HannVector;

	public:

	ComputeHanning(std::vector<double>* ptr_HannVector);
	
	int _nsample() { return nsample; }
	double _indenom() { return indenom; }
	void generate_Hanning();
	double get_Hanning_ave(std::vector<double>* ptr_DataVector);
	
};
#endif
