#ifndef ___CLASS_COMPUTE_CMAP
#define ___CLASS_COMPUTE_CMAP
#include "LoadParm.hpp"
class ComputeCmap
{
	private:

	static const int res               = 24;
	static const int half_res          = res / 2;
	static const int two_res           = res * 2;
	static constexpr double stepSize   = 15.0;
	static constexpr double gridOrigin = -180.0;

	LoadParm* ptr_All22;
	std::vector<Cmap>* ptr_cmapVector;
	

	public:

	ComputeCmap(){};
	ComputeCmap(LoadParm& All22, std::vector<Cmap>& cmapVector);

	double compute_force();

	double calc_cmap(const double phi, const double psi, const int cmap_type, double& dEdPhi, double& dEdPsi);

	private:
	void generate_cmap_grid_derivatives();
};
#endif
