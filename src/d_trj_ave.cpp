#include <iostream>
#include <fstream>
#include <sstream>
#include "PSF.hpp"
#include "ReadDCD.hpp"
#include "NAMDBin.hpp"
using namespace std;
int main(int argc, char** argv)
{
	if (argc != 4)
	{
		cout << 
		"\nD_TRJ_AVE\n"
		"\nCalculate average structure\n"
		"\nusage: ./a.out psf[natom] dcd ofilename\n\n";
		return 1;
	}

	string ofilename(argv[3]);
	cout << "REMARK ";
	for (int i = 0; i < argc; i++)
		cout << argv[i] << " ";
	cout <<
	"\nREMARK D_TRJ_AVE\n"
	"REMARK Average structure calculation\n"
	"REMARK Output to file " << ofilename << '\n';

	int natom = 0;
	vector<Atom> atomVector;

	string spsf(argv[1]);
	string spsfex = spsf.substr(spsf.find_last_of(".") + 1);
	if (spsfex == "psf")
	{
		PSF PSFFile(argv[1]);
		atomVector = PSFFile.atomVector;
		natom = atomVector.size();
	}
	else
	{
		natom = atoi(argv[1]);
		atomVector.resize(natom);
	}

	//  average structure coordinates  
	vector<Eigen::Vector3d> aveCoords;

	for (int i = 0; i < natom; i++)
	{
		Eigen::Vector3d v(0., 0., 0.);
		aveCoords.push_back(v);
	}
	
	ReadDCD TRJFile(&atomVector);
	if (!TRJFile.open_dcd_read(argv[2]))
		return 0;
	
	while(TRJFile.read_1step())
	{
		for (int i = 0; i < natom; i++)
		{
			aveCoords[i] += atomVector[i].position;
		}
	}
	
	double dstep = static_cast<double>(TRJFile._nsteps());

	//  average over step
	for (int i = 0; i < natom; i++)
	{
		atomVector[i].position = aveCoords[i] / dstep;
	}

	// set some valuables to 0.
	for (int i = 0; i < natom; i++)
	{
		atomVector[i].PDBo = 0.;
		atomVector[i].PDBb = 0.;
	}

	string header = "generated by d_trj_ave.x \nREMARK AVERAGED STRUCTURE FROM : ";
	header += TRJFile._ifilename();
	ostringstream os;
	os << TRJFile._nsteps();
	header += " -> total " + os.str() + " steps.";

	if (spsfex == "psf")
	{
		PSF PSFFile(argv[1], &atomVector);
		PSFFile.writePDB(argv[3], header);
	}
	else
	{
		NAMDBin JOB(argv[3], "coor");
		JOB.write_fo(atomVector);
	}
}
