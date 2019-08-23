#include <iostream>
#include "PSF.hpp"
#include "PDB.hpp"
#include "ForceUtil.hpp"
using namespace std;
int main (int argc, char** argv)
{
	if (argc < 3)
	{
		cout << "\nusage: ./a.out psf ave\n\n";
		return 1;
	}

	PSF PSFFile(argv[1]);
	PDB PDBFile(argv[2]);
	if (!PDBFile.LoadCoords(PSFFile.atomVector))
		return 1;

	int iD92, iK95, iD195, iK198, iD291, iK294;

	for (int i = 0; i < PSFFile.atomVector.size(); i++)
	{
		Atom& at = PSFFile.atomVector[i];

		if (at.PDBAtomName != "CA") continue;

		switch (at.PSFResID) {
			case  92: iD92  = i; break;
			case  95: iK95  = i; break;
			case 195: iD195 = i; break;
			case 198: iK198 = i; break;
			case 291: iD291 = i; break;
			case 294: iK294 = i; break;
		}
	}

	Atom& D92 = PSFFile.atomVector[iD92];
	Atom& K95 = PSFFile.atomVector[iK95];
	Atom& D195 = PSFFile.atomVector[iD195];
	Atom& K198 = PSFFile.atomVector[iK198];
	Atom& D291 = PSFFile.atomVector[iD291];
	Atom& K294 = PSFFile.atomVector[iK294];

	vector<Eigen::Vector3d> vfrc;
	for (int i = 0; i < PSFFile.atomVector.size(); i++)
	{
		Eigen::Vector3d vtmp(0., 0., 0.);
		vfrc.push_back(vtmp);
	}

	Eigen::Vector3d d1 = K294.position - D92.position;
	Eigen::Vector3d d2 =  K95.position - D195.position;
	Eigen::Vector3d d3 = K198.position - D291.position;

	vfrc[iD92]  =  d1;
	vfrc[iK294] = -d1;

	vfrc[iD195] =  d2;
	vfrc[iK95]  = -d2;
	
	vfrc[iD291] =  d3;
	vfrc[iK198] = -d3;

	ForceUtil FRC;
	for (int i = 0; i < PSFFile.atomVector.size(); i++)
	{
		if (PSFFile.atomVector[i].PDBAtomName == "CA")
		FRC.writeForce("stdout", PSFFile.atomVector[i], vfrc[i]);
	}

}
