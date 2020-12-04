#include <iostream>
#include <iomanip>
#include <cstdlib>
#include "PSF.hpp"
#include "ReadDCD.hpp"
#include "ComputeLipid.hpp"

using namespace std;

int main(int argc, char** argv)
{
	if (argc < 7)
	{
		cout << 
		"\nCompute Lipid Density Profile\n"
		"\nusage: ./a.out psf dcd min max bin mode\n\n"
		"      mode = 1: number density\n"
		"      mode = 2:   mass density\n\n";
		return 1;
	}

	output_args(argc, argv);

	vector<Atom> atomVector;
	PSF PSFFile(argv[1], &atomVector);
	ReadDCD DCD(&atomVector);
	if (!DCD.open_dcd_read(argv[2]))
		return 0;

	double vmin = atof(argv[3]);
	double vmax = atof(argv[4]);
	int    nbin = atoi(argv[5]);
	unsigned int mode = atoi(argv[6]);

	ComputeLipid JOB(vmin, vmax, nbin, mode, &atomVector);
	JOB.set_selection();

	cout << "REMARK ============\n";
	
	cout << setprecision(8) << scientific;
	while (DCD.read_1step())
	{
		JOB.calc_histogram(DCD._boxx(), DCD._boxy());
	}

	JOB.calc_density_profile(DCD._nsets());

	JOB.output_density_profile();
}
