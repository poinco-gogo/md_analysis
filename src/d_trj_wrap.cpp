#include <iostream>
#include <cstdlib>
#include <iomanip>
#include <sstream>
#include "PSF.hpp"
#include "ReadDCD.hpp"
#include "Lattice.hpp"
using namespace std;
void wrap_molecule(const Lattice& lattice, const vector<int>& ind, const Eigen::Vector3d& com, PSF& psf);
int main (int argc, char** argv)
{
	if (argc < 5)
	{
		cout << "\nusage: ./a.out psf dcd stride center\n\n";
		return 0;
	}
	cout << "REMARK ";
	for (int i = 0; i < argc; i++)
		cout << argv[i] << ' ';
	cout << '\n';


	PSF PSFFile(argv[1]);
	ReadDCD DCD(&PSFFile.atomVector);

	if (!DCD.open_dcd_read(argv[2]))
		return 1;
	
	int nstride = atoi(argv[3]);
	cout << "REMARK Calculation will be done every " << nstride << " step(s).\n";
	string center(argv[4]);
	if (center == "origin" || center == "pcom")
	{
		cout << "REMARK COODINATES WILL BE WRAPPED AROUND "
			<< center << ".\n";
	}
	else
	{
		cerr << "\nerror: unknown keyword \"" << center
			<< "\"\n\n";
		return 1;
	}

	DCD.open_dcd_write("outdcd");

	vector<int> iprt, iprt_hv, iADP, iGDP, iCAL, iBON, iBEN, iATP;
	vector<int> iATPA, iATPC, iATPD, iHEME, iADPC, iADPD;


	vector<int> iPL, iOH2;
	for (int i = 0; i < PSFFile.atomVector.size(); i++)
	{	
		if (PSFFile.atomVector[i].PSFAtomName == "PL")
			iPL.push_back(i);
		if (PSFFile.atomVector[i].PDBAtomName == "OH2")
			iOH2.push_back(i);

	}
	for (int i = 0; i < PSFFile.atomVector.size(); i++)
	{
		Atom& atom = PSFFile.atomVector[i];

		const string res = atom.PSFResName;
		const string seg = atom.PSFSegmentName;

		if (seg == "PRT" 
			|| seg == "PRA"
			|| seg == "PRB"
			|| seg == "PRC"
			|| seg == "PRD"
			|| seg == "PRI"
			|| seg == "PRJ"
			|| seg == "PRK"
			|| seg == "PRL")
		{
			iprt.push_back(i);

			if (atom.mass < 2.)
				continue;
			else
				iprt_hv.push_back(i);
		}
		else if (res == "CAL")
		{
			iCAL.push_back(i);
		}
		else if (seg == "ADP")
		{
			iADP.push_back(i);
		}
		else if (seg == "ADPC")
		{
			iADPC.push_back(i);
		}
		else if (seg == "ADPD")
		{
			iADPD.push_back(i);
		}
		else if (seg == "ATP")
		{
			iATP.push_back(i);
		}
		else if (seg == "ATPA")
		{
			iATPA.push_back(i);
		}
		else if (seg == "ATPC")
		{
			iATPC.push_back(i);
		}
		else if (seg == "ATPD")
		{
			iATPD.push_back(i);
		}
		else if (seg == "GDP")
		{
			iGDP.push_back(i);
		}
		else if (seg == "BKA" || seg == "BKC")
		{
			iBON.push_back(i);
		}
		else if (seg == "BEN")
		{
			iBEN.push_back(i);
		}
		else if (seg == "HEME")
		{
			iHEME.push_back(i);
		}
	}

	cout << "REMARK " << iprt_hv.size() << " protein heavy atoms found.\n";
	cout << "REMARK " << iprt.size() << " protein       atoms found.\n";
	cout << "REMARK " << iADP.size() << " ADP       atoms found.\n";
	cout << "REMARK " << iADPC.size() << " ADPC      atoms found.\n";
	cout << "REMARK " << iADPD.size() << " ADPD      atoms found.\n";
	cout << "REMARK " << iATP.size() << " ATP       atoms found.\n";
	cout << "REMARK " << iATPA.size() << " ATPA      atoms found.\n";
	cout << "REMARK " << iATPC.size() << " ATPC      atoms found.\n";
	cout << "REMARK " << iATPD.size() << " ATPD      atoms found.\n";
	cout << "REMARK " << iGDP.size() << " GDP       atoms found.\n";
	cout << "REMARK " << iCAL.size() << " CAL       atoms found.\n";
	cout << "REMARK " << iBON.size() << " BON       atoms found.\n";
	cout << "REMARK " << iBEN.size() << " BEN       atoms found.\n";
	cout << "REMARK " << iHEME.size() << " HEME      atoms found.\n";

	const int atoms_per_heme = 73;
	const int num_heme = iHEME.size() / atoms_per_heme;

	int icnt = 0;
	while (DCD.read_1step())
	{
		if (DCD._nsteps() % nstride)
			continue;

		// get box info
		const double x = DCD._boxx();
		const double y = DCD._boxy();
		const double z = DCD._boxz();

		int icnt_heme = 0;

		if (!x || !y || !z)
		{
			cerr << "error: cannot read box info from dcd.\n";
			return 1;
		}

		Lattice lattice(x, y, z);

		Eigen::Vector3d com(0., 0., 0.);
		if (center == "pcom")
			com = PSFFile.getcomof_zerobased(iprt_hv);
		else
			// wrap protein around origin
			wrap_molecule(lattice, iprt, com, PSFFile);

		for (int i = 0; i < PSFFile.atomVector.size(); i++)
		{
			Atom& atom = PSFFile.atomVector[i];

			const string res = atom.PSFResName;
			const string seg = atom.PSFSegmentName;

			if (seg == "PRT")
				continue;

			if (atom.PDBAtomName == "OH2")
			{
				Eigen::Vector3d del = atom.position - com;
				atom.position += lattice.wrap_delta(del);
				Atom& H1 = PSFFile.atomVector[i + 1];
				Atom& H2 = PSFFile.atomVector[i + 2];
				H1.position += lattice.wrap_delta(del);
				H2.position += lattice.wrap_delta(del);
			}
			else if (res == "POPG" && atom.PDBAtomName == "C13")
			{
				int nlipid = 127; // # per lipid
				Eigen::Vector3d lipid_com(0., 0., 0.);
				for (int j = 0; j < nlipid; j++)
				{
					lipid_com += PSFFile.atomVector[i + j].position;
				}
				lipid_com = lipid_com / nlipid;
				Eigen::Vector3d del = lipid_com - com;
				for (int j = 0; j < nlipid; j++)
				{
					Atom& lp = PSFFile.atomVector[i + j];
					lp.position += lattice.wrap_delta(del);
				}
			}
			else if (res == "POPE" && atom.PSFAtomName == "NH3L")
			{
				int nlipid = 125; // # per lipid
				Eigen::Vector3d lipid_com(0., 0., 0.);
				for (int j = 0; j < nlipid; j++)
				{
					lipid_com += PSFFile.atomVector[i + j].position;
				}
				lipid_com = lipid_com / nlipid;
				Eigen::Vector3d del = lipid_com - com;
				for (int j = 0; j < nlipid; j++)
				{
					Atom& lp = PSFFile.atomVector[i + j];
					lp.position += lattice.wrap_delta(del);
				}
			}
			else if (res == "POPC" && atom.PSFAtomName == "NTL")
			{
				int nlipid = 134; // # per lipid
				Eigen::Vector3d lipid_com(0., 0., 0.);
				for (int j = 0; j < nlipid; j++)
				{
					lipid_com += PSFFile.atomVector[i + j].position;
				}
				lipid_com = lipid_com / nlipid;
				Eigen::Vector3d del = lipid_com - com;
				for (int j = 0; j < nlipid; j++)
				{
					Atom& lp = PSFFile.atomVector[i + j];
					lp.position += lattice.wrap_delta(del);
				}
			}
			else if (res == "DPPE" && atom.PSFAtomName == "NH3L")
			{
				int nlipid = 121; // # per lipid
				Eigen::Vector3d lipid_com(0., 0., 0.);
				for (int j = 0; j < nlipid; j++)
				{
					lipid_com += PSFFile.atomVector[i + j].position;
				}
				lipid_com = lipid_com / nlipid;
				Eigen::Vector3d del = lipid_com - com;
				for (int j = 0; j < nlipid; j++)
				{
					Atom& lp = PSFFile.atomVector[i + j];
					lp.position += lattice.wrap_delta(del);
				}
			}
			else if (res == "DPPG" && atom.PDBAtomName == "C13")
			{
				int nlipid = 123; // # per lipid
				Eigen::Vector3d lipid_com(0., 0., 0.);
				for (int j = 0; j < nlipid; j++)
				{
					lipid_com += PSFFile.atomVector[i + j].position;
				}
				lipid_com = lipid_com / nlipid;
				Eigen::Vector3d del = lipid_com - com;
				for (int j = 0; j < nlipid; j++)
				{
					Atom& lp = PSFFile.atomVector[i + j];
					lp.position += lattice.wrap_delta(del);
				}
			}
			else if (res == "PYPE" && atom.PSFAtomName == "NH3L")
			{
				int nlipid = 119; // # per lipid
				Eigen::Vector3d lipid_com(0., 0., 0.);
				for (int j = 0; j < nlipid; j++)
				{
					lipid_com += PSFFile.atomVector[i + j].position;
				}
				lipid_com = lipid_com / nlipid;
				Eigen::Vector3d del = lipid_com - com;
				for (int j = 0; j < nlipid; j++)
				{
					Atom& lp = PSFFile.atomVector[i + j];
					lp.position += lattice.wrap_delta(del);
				}
			}
			else if (res == "PYPG" && atom.PDBAtomName == "C13")
			{
				int nlipid = 121; // # per lipid
				Eigen::Vector3d lipid_com(0., 0., 0.);
				for (int j = 0; j < nlipid; j++)
				{
					lipid_com += PSFFile.atomVector[i + j].position;
				}
				lipid_com = lipid_com / nlipid;
				Eigen::Vector3d del = lipid_com - com;
				for (int j = 0; j < nlipid; j++)
				{
					Atom& lp = PSFFile.atomVector[i + j];
					lp.position += lattice.wrap_delta(del);
				}
			}
			else if (atom.PDBAtomName == "C4'" && res == "ADP")
			{
				wrap_molecule(lattice, iADP, com, PSFFile);
			}
			else if (atom.PDBAtomName == "C4'" && res == "ADPC")
			{
				wrap_molecule(lattice, iADPC, com, PSFFile);
			}
			else if (atom.PDBAtomName == "C4'" && res == "ADPD")
			{
				wrap_molecule(lattice, iADPD, com, PSFFile);
			}
			else if (atom.PDBAtomName == "C4'" && seg == "ATP")
			{
				wrap_molecule(lattice, iATP, com, PSFFile);
			}
			else if (atom.PDBAtomName == "C4'" && seg == "ATPA")
			{
				wrap_molecule(lattice, iATPA, com, PSFFile);
			}
			else if (atom.PDBAtomName == "C4'" && seg == "ATPC")
			{
				wrap_molecule(lattice, iATPC, com, PSFFile);
			}
			else if (atom.PDBAtomName == "C4'" && seg == "ATPD")
			{
				wrap_molecule(lattice, iATPD, com, PSFFile);
			}
			else if (atom.PDBAtomName == "N9" && res == "GDP")
			{
				wrap_molecule(lattice, iGDP, com, PSFFile);
			}
			else if ((seg == "BKA" || seg == "BKC")
					&& atom.PDBAtomName == "C1")
			{
				wrap_molecule(lattice, iBON, com, PSFFile);
			}
			else if (seg == "BEN" && atom.PDBAtomName == "C1")
			{
				wrap_molecule(lattice, iBEN, com, PSFFile);
			}
			else if (seg == "HEME" && atom.PDBAtomName == "FE")
			{
				Eigen::Vector3d bcom(0., 0., 0.);
				for (int j = 0; j < atoms_per_heme; j++)
				{
					int k = j + atoms_per_heme * icnt_heme;
					bcom += PSFFile.atomVector[iHEME[k]].position;
				}
				bcom = bcom / atoms_per_heme;
				Eigen::Vector3d dr = bcom - com;
				for (int j = 0; j < atoms_per_heme; j++)
				{
					int k = j + atoms_per_heme * icnt_heme;
					PSFFile.atomVector[iHEME[k]].position
						+= lattice.wrap_delta(dr);
				}

				++icnt_heme;
			}
			else if (seg == "ION")
			{
				Eigen::Vector3d dr = atom.position - com;
				atom.position += lattice.wrap_delta(dr);
			}
			else if (atom.PDBAtomName == "CLA")
			{
				Eigen::Vector3d dr = atom.position - com;
				atom.position += lattice.wrap_delta(dr);
			}
			else if (atom.PDBAtomName == "MG")
			{
				Eigen::Vector3d dr = atom.position - com;
				atom.position += lattice.wrap_delta(dr);
			}
		}
		// このスナップショットの巻き直し完了

		DCD.write_1step();
		
		++icnt;
	}
	cout << "REMARK " << icnt << " calculations done.\n";
}

void wrap_molecule(const Lattice& lattice, const vector<int>& ind, const Eigen::Vector3d& com, PSF& psf)
{
	Eigen::Vector3d vcom = psf.getcomof_zerobased(ind);

	Eigen::Vector3d del  = vcom - com;

	for (auto& i: ind)
	{
		Atom& at = psf.atomVector[i];

		at.position += lattice.wrap_delta(del);
	}
}
