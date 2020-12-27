#include <iostream>
#include <fstream>
#include <sstream>
#include "LoadAmberPrmtop.hpp"
using namespace std;
void show_sm(string s, int m)
{
	cout << "REMARK " << s << m << '\n';
}

template <class X> bool show_err(string s, int n, vector<X>& vec)
{
	if (vec.size() != n)
	{
		cerr << "error: " << s << '\n';
		return false;
	}
	return true;
}

bool show_err_3(int n, vector<int>& v1, vector<int>& v2, vector<int>& v3)
{
	if (v1.size() != n || v2.size() != n || v3.size() != n)
		return false;
	else
		return true;
}

void load_20a4(ifstream& fi, int numtot, vector<string>& vec)
{
	string s;
	for (int i = 0; i < numtot / 20; i++)
	{
		getline(fi, s);
		for (int ii = 0; ii < 20; ii++)
		{
			istringstream is(s.substr(4 * ii, 4));
			vec.push_back(is.str());
		}
	}
	if (numtot % 20)
	{
		getline(fi, s);
		for (int i = 0; i < numtot % 20; i++)
		{
			istringstream is(s.substr(4 * i, 4));
			vec.push_back(is.str());
		}
	}
}

void load_5e16(ifstream& fi, int numtot, vector<double>& vec)
{
	string s;
	for (int i = 0; i < numtot / 5; i++)
	{
		getline(fi, s);
		istringstream is(s);
		for (int ii = 0; ii < 5; ii++)
		{
			double dtmp;
			is >> dtmp;
			vec.push_back(dtmp);
		}
	}
	if (numtot % 5)
	{
		getline(fi, s);
		istringstream is(s);
		for (int i = 0; i < numtot % 5; i++)
		{
			double dtmp;
			is >> dtmp;
			vec.push_back(dtmp);
		}
	}
}

void load_10i8(ifstream& fi, int numtot, vector<int>& vec)
{
	string s;
	for (int i = 0; i < numtot / 10; i++)
	{
		getline(fi, s);
		istringstream is(s);
		for (int ii = 0; ii < 10; ii++)
		{
			int itmp;
			is >> itmp;
			vec.push_back(itmp);
		}
	}
	if (numtot % 10)
	{
		getline(fi, s);
		istringstream is(s);
		for (int i = 0; i < numtot % 10; i++)
		{
			int itmp;
			is >> itmp;
			vec.push_back(itmp);
		}
	}
}

void load_10i8_n(ifstream& fi, int numtot, vector<int>& vtmp, int n)
{
	string s;
	for (int i = 0; i < numtot * n / 10; i++)
	{
		getline(fi, s);
		istringstream is(s);
		for (int ii = 0; ii < 10; ii++)
		{
			int itmp;
			is >> itmp;
			vtmp.push_back(itmp);
		}
	}
	
	if ((numtot * n) % 10)
	{
		getline(fi, s);
		istringstream is(s);
		for (int i = 0; i < (numtot * n) % 10; i++)
		{
			int itmp;
			is >> itmp;
			vtmp.push_back(itmp);
		}
	}
}

void load_10i8_3(ifstream& fi, int numtot, vector<int>& v1, vector<int>& v2, vector<int>& v3)
{
	vector<int> vtmp;

	load_10i8_n(fi, numtot, vtmp, 3);
	
	for (int i = 0; i < numtot; i++)
	{
		v1.push_back(vtmp[3 * i    ]);
		v2.push_back(vtmp[3 * i + 1]);
		v3.push_back(vtmp[3 * i + 2]);
	}
}

void load_10i8_4(ifstream& fi, int numtot, vector<int>& v1, vector<int>& v2, vector<int>& v3, vector<int>& v4)
{
	vector<int> vtmp;

	load_10i8_n(fi, numtot, vtmp, 4);

	for (int i = 0; i < numtot; i++)
	{
		v1.push_back(vtmp[4 * i    ]);
		v2.push_back(vtmp[4 * i + 1]);
		v3.push_back(vtmp[4 * i + 2]);
		v4.push_back(vtmp[4 * i + 3]);
	}
}

void load_10i8_5(ifstream& fi, int numtot, vector<int>& v1, vector<int>& v2, vector<int>& v3, vector<int>& v4, vector<int>& v5)
{
	vector<int> vtmp;

	load_10i8_n(fi, numtot, vtmp, 5);

	for (int i = 0; i < numtot; i++)
	{
		v1.push_back(vtmp[5 * i    ]);
		v2.push_back(vtmp[5 * i + 1]);
		v3.push_back(vtmp[5 * i + 2]);
		v4.push_back(vtmp[5 * i + 3]);
		v5.push_back(vtmp[5 * i + 4]);
	}
}

LoadAmberPrmtop::LoadAmberPrmtop(string filename)
{
	ifstream fi(filename.c_str());
	cout << "REMARK Opening prmtop: " << filename << '\n';
	load_prmtop(fi);
}

void LoadAmberPrmtop::load_prmtop(ifstream& fi)
{
	string s;
	int i = 0;
	while (getline(fi, s))
	{
		if (s.find("FLAG TITLE", 0) != string::npos)
		{
			getline(fi, s);
			getline(fi, s);
			cout << "REMARK TITLE: " << 
				s.substr(0, s.find(" ", 0)) << '\n';
		}
		if (s.find("FLAG POINTERS", 0) != string::npos)
		{
			getline(fi, s);
			getline(fi, s);
			istringstream is1(s);
			is1 >> natom >> ntypes >> nbonh >> mbona >> ntheth >> mtheta >> nphih >> mphia >> nhparm >> nparm;
			getline(fi, s);
			istringstream is2(s);
			is2 >> nnb >> nres >> nbona >> ntheta >> nphia >> numbnd >> numang >> nptra >> natyp >> nphb;
			getline(fi, s);
			istringstream is3(s);
			is3 >> ifpert >> nbper >> ngper >> ndper >> mbper >> mgper >> mdper >> ifbox >> nmxrs >> ifcap;
			getline(fi, s);
			istringstream is4(s);
			is4 >> numextra;
			
			show_sm("Total Number of Atoms: ", natom);
			show_sm("Total Number of Distinct Atom Types: ", ntypes);
			show_sm("Number of Bonds w/  hydrogen: ", nbonh);
			show_sm("Number of Bonds w/o hydrogen: ", mbona);
			show_sm("Number of Angles w/  hydrogen: ", ntheth);
			show_sm("Number of Angles w/o hydrogen: ", mtheta);
			show_sm("Number of Dihedrals w/  hydrogen: ", nphih);
			show_sm("Number of Dihedrals w/o hydrogen: ", mphia);
			show_sm("Currently not used: ", nhparm);
			show_sm("This prmtop was created by ADDLES: ", nparm);
			show_sm("Number of Excluded Atoms: ", nnb);
			show_sm("Number of Residues: ", nres);
			show_sm("MBONA + Number of Constraint Bonds: ", nbona);
			show_sm("MTHETA + Number of Constraint Angles: ", ntheta);
			show_sm("MPHIA + Number of Constraint Dihedrals: ", nphia);
			show_sm("Number of Unique Bond Types: ", numbnd);
			show_sm("Number of Unique Angle Types: ", numang);
			show_sm("Number of Unique Dihedral Types: ", nptra);
			show_sm("Number of Atom Types in Parameter File: ", natyp);
			show_sm("Number of Distinct 10-12 Hydrogen Bond Pair Types: ", nphb);
			show_sm("Perturbation Info to be read: ", ifpert);
			show_sm("Number of Bonds to be Perturbed: ", nbper);
			show_sm("Number of Angles to be Perturbed: ", ngper);
			show_sm("Number of Dihedrals to be Perturbed: ", ndper);
			show_sm("Number of Bonds with Atoms Completely in Perturbed Group: ", mbper);
			show_sm("Number of Angles with Atoms Completely in Perturbed Group: ", mgper);
			show_sm("Number of Dihedrals with Atoms Completely in Perturbed Group: ", mdper);
			show_sm("1 for Standard Periodic Box: ", ifbox);
			show_sm("Number of Atoms in the Largest Residue: ", nmxrs);
			show_sm("CAP Option from Edit was Specified: ", ifcap);
			show_sm("Number of Extra Points Found in Topology: ", numextra);

		}
		if (s.find("FLAG ATOM_NAME", 0) != string::npos)
		{
			getline(fi, s);

			load_20a4(fi, natom, igraph);
			
			if (!show_err("igraph.size() != natom", natom, igraph))
			       	break;
		}
		if (s.find("FLAG CHARGE", 0) != string::npos)
		{
			getline(fi, s);

			load_5e16(fi, natom, charge);
			
			if (!show_err("charge.size() != natom", natom, charge))
				break;	
		}
		if (s.find("FLAG MASS", 0) != string::npos)
		{
			getline(fi, s);

			load_5e16(fi, natom, amass);
			
			if (!show_err("amass.size() != natom", natom, amass))
				break;	
		}
		if (s.find("FLAG ATOM_TYPE_INDEX", 0) != string::npos)
		{
			getline(fi, s);

			load_10i8(fi, natom, iac);
			
			if (!show_err("iac.size() != natom", natom, iac))
				break;	
		}
		if (s.find("FLAG NUMBER_EXCLUDED_ATOMS", 0) != string::npos)
		{
			getline(fi, s);

			load_10i8(fi, natom, numex);
			
			if (!show_err("numex.size() != natom", natom, numex))
				break;	
		}
		if (s.find("FLAG NONBONDED_PARM_INDEX", 0) != string::npos)
		{
			getline(fi, s);

			load_10i8(fi, ntypes * ntypes, ico);
			
			if (!show_err("ico.size() != ntypes**2", ntypes*ntypes, ico))
				break;	
		}
		if (s.find("FLAG RESIDUE_LABEL", 0) != string::npos)
		{
			getline(fi, s);

			load_20a4(fi, nres, lbres);
			
			if (!show_err("lbres.size() != nres", nres, lbres))
				break;	
		}
		if (s.find("FLAG RESIDUE_POINTER", 0) != string::npos)
		{
			getline(fi, s);

			load_10i8(fi, nres, ipres);
			
			if (!show_err("ipres.size() != nres", nres, ipres))
				break;	
		}
		if (s.find("FLAG BOND_FORCE_CONSTANT", 0) != string::npos)
		{
			getline(fi, s);

			load_5e16(fi, numbnd, rk);
			
			if (!show_err("rk.size() != numbnd", numbnd, rk))
				break;	
		}
		if (s.find("FLAG BOND_EQUIL_VALUE", 0) != string::npos)
		{
			getline(fi, s);

			load_5e16(fi, numbnd, req);
			
			if (!show_err("req.size() != numbnd", numbnd, req))
				break;	
		}
		if (s.find("FLAG ANGLE_FORCE_CONSTANT", 0) != string::npos)
		{
			getline(fi, s);

			load_5e16(fi, numang, tk);
			
			if (!show_err("tk.size() != numang", numang, tk))
				break;	
		}
		if (s.find("FLAG ANGLE_EQUIL_VALUE", 0) != string::npos)
		{
			getline(fi, s);

			load_5e16(fi, numang, teq);
			
			if (!show_err("teq.size() != numang", numang, teq))
				break;	
		}
		if (s.find("FLAG DIHEDRAL_FORCE_CONSTANT", 0) != string::npos)
		{
			getline(fi, s);

			load_5e16(fi, nptra, pk);
			
			if (!show_err("pk.size() != nptra", nptra, pk))
				break;	
		}
		if (s.find("FLAG DIHEDRAL_PERIODICITY", 0) != string::npos)
		{
			getline(fi, s);

			load_5e16(fi, nptra, pn);
			
			if (!show_err("pn.size() != nptra", nptra, pn))
				break;	
		}
		if (s.find("FLAG DIHEDRAL_PHASE", 0) != string::npos)
		{
			getline(fi, s);

			load_5e16(fi, nptra, phase);
			
			if (!show_err("phase.size() != nptra", nptra, phase))
				break;	
		}
		if (s.find("FLAG SCEE_SCALE_FACTOR", 0) != string::npos)
		{
			getline(fi, s);

			load_5e16(fi, nptra, one_scee);
			
			if (!show_err("one_scee.size() != nptra", nptra, one_scee))
				break;	
		}
		if (s.find("FLAG SCNB_SCALE_FACTOR", 0) != string::npos)
		{
			getline(fi, s);

			load_5e16(fi, nptra, one_scnb);
			
			if (!show_err("one_scnb.size() != nptra", nptra, one_scnb))
				break;	
		}
		if (s.find("FLAG SOLTY", 0) != string::npos)
		{
			getline(fi, s);

			load_5e16(fi, natyp, solty);
			
			if (!show_err("solty.size() != nptra", natyp, solty))
				break;	
		}
		if (s.find("FLAG LENNARD_JONES_ACOEF", 0) != string::npos)
		{
			getline(fi, s);

			int itmp = ntypes * (ntypes + 1) / 2;
			load_5e16(fi, itmp, cn1);
			
			if (!show_err("cn1.size() != ntypes*(ntypes+1)/2", itmp, cn1))
				break;	
		}
		if (s.find("FLAG LENNARD_JONES_BCOEF", 0) != string::npos)
		{
			getline(fi, s);

			int itmp = ntypes * (ntypes + 1) / 2;
			load_5e16(fi, itmp, cn2);
			
			if (!show_err("cn2.size() != ntypes*(ntypes+1)/2", itmp, cn2))
				break;	
		}
		if (s.find("FLAG BONDS_INC_HYDROGEN", 0) != string::npos)
		{
			getline(fi, s);

			load_10i8_3(fi, nbonh, ibh, jbh, icbh);
			
			if (!show_err_3(nbonh, ibh, jbh, icbh))
			{
				cerr << "error: in FLAG BONDS_INC_HYDROGEN\n";
				break;	
			}
		}
		if (s.find("FLAG BONDS_WITHOUT_HYDROGEN", 0) != string::npos)
		{
			getline(fi, s);

			load_10i8_3(fi, nbona, ib, jb, icb);
			
			if (!show_err_3(nbona, ib, jb, icb))
			{
				cerr << "error: in FLAG BONDS_WITHOUT_HYDROGEN\n";
				break;	
			}
		}
		if (s.find("FLAG ANGLES_INC_HYDROGEN", 0) != string::npos)
		{
			getline(fi, s);

			load_10i8_4(fi, ntheth, ith, jth, kth, icth);
			
			if (ith.size() != ntheth || jth.size() != ntheth || kth.size() != ntheth || icth.size() != ntheth)
			{
				cerr << "error: in FLAG ANGLES_INC_HYDROGEN\n";
				break;	
			}
		}
		if (s.find("FLAG ANGLES_WITHOUT_HYDROGEN", 0) != string::npos)
		{
			getline(fi, s);

			load_10i8_4(fi, ntheta, it, jt, kt, ict);
			
			if (it.size() != ntheta || jt.size() != ntheta || kt.size() != ntheta || ict.size() != ntheta)
			{
				cerr << "error: in FLAG ANGLES_WITHOUT_HYDROGEN\n";
				break;	
			}
		}
		if (s.find("FLAG DIHEDRALS_INC_HYDROGEN", 0) != string::npos)
		{
			getline(fi, s);

			load_10i8_5(fi, nphih, iph, jph, kph, lph, icph);
			
			if (iph.size() != nphih || jph.size() != nphih || kph.size() != nphih || lph.size() != nphih || icph.size() != nphih)
			{
				cerr << "error: in FLAG DIHEDRALS_INC_HYDROGEN\n";
				break;	
			}
		}
		if (s.find("FLAG DIHEDRALS_WITHOUT_HYDROGEN", 0) != string::npos)
		{
			getline(fi, s);

			load_10i8_5(fi, nphia, ip, jp, kp, lp, icp);
			
			if (ip.size() != nphia || jp.size() != nphia || kp.size() != nphia || lp.size() != nphia || icp.size() != nphia)
			{
				cerr << "error: in FLAG DIHEDRALS_WITHOUT_HYDROGEN\n";
				break;	
			}
		}
		if (s.find("FLAG EXCLUDED_ATOMS_LIST", 0) != string::npos)
		{
			getline(fi, s);

			load_10i8(fi, nnb, inb);
			
			if (!show_err("inb.size() != nnb", nnb, inb))
				break;	
		}
		if (s.find("FLAG HBOND_ACOEF", 0) != string::npos)
		{
			getline(fi, s);

			load_5e16(fi, nphb, asol);
			
			if (!show_err("asol.size() != nphb", nphb, asol))
				break;	
		}
		if (s.find("FLAG HBOND_BCOEF", 0) != string::npos)
		{
			getline(fi, s);

			load_5e16(fi, nphb, bsol);
			
			if (!show_err("bsol.size() != nphb", nphb, bsol))
				break;	
		}
		if (s.find("FLAG HBCUT", 0) != string::npos)
		{
			getline(fi, s);

			load_5e16(fi, nphb, hbcut);
			
			if (!show_err("hbcut.size() != nphb", nphb, hbcut))
				break;	
		}
		if (s.find("FLAG AMBER_ATOM_TYPE", 0) != string::npos)
		{
			getline(fi, s);

			load_20a4(fi, natom, isymbl);
			
			if (!show_err("isymbl.size() != natom", natom, isymbl))
				break;	
		}
		if (s.find("FLAG TREE_CHAIN_CLASSIFICATION", 0) != string::npos)
		{
			getline(fi, s);

			load_20a4(fi, natom, itree);
			
			if (!show_err("itree.size() != natom", natom, itree))
				break;	
		}
		if (s.find("FLAG JOIN_ARRAY", 0) != string::npos)
		{
			getline(fi, s);

			load_10i8(fi, natom, join);
			
			if (!show_err("join.size() != natom", natom, join))
				break;	
		}
		if (s.find("FLAG IROTAT", 0) != string::npos)
		{
			getline(fi, s);

			load_10i8(fi, natom, irotat);
			
			if (!show_err("irotat.size() != natom", natom, irotat))
				break;	
		}
		if (s.find("FLAG SOLVENT_POINTERS", 0) != string::npos)
		{
			getline(fi, s);
			getline(fi, s);
			istringstream is(s);
			is >> iptres >> nspm >> nspsol;
		}
		if (s.find("FLAG ATOMS_PER_MOLECULE", 0) != string::npos)
		{
			getline(fi, s);

			load_10i8(fi, nspm, nsp);
			
			if (!show_err("nsp.size() != nspm", nspm, nsp))
				break;	
		}
		if (s.find("FLAG BOX_DIMENSIONS", 0) != string::npos)
		{
			getline(fi, s);
			getline(fi, s);
			istringstream is(s);
			box.resize(3);
			is >> oldbeta >> box[0] >> box[1] >> box[2];
		}
		
	}
}

int LoadAmberPrmtop::get_num_nwt()
{
	int numwat = 0;
	int ipcnt  = 0;

        for (int i = 0; i < natom; i++)
	{
		// get residue name for atom with index i
		if (lbres[ipcnt] == "WAT ") ++numwat;
		if (i + 1 == ipres[ipcnt + 1] - 1) ++ipcnt;
	}
	return natom - numwat;
}

void LoadAmberPrmtop::make_atomVector(vector<Atom>& atomVector)
{
	atomVector.resize(natom);

	for (int i = 0; i < natom; i++)
	{
		atomVector[i].PDBIndex    = i + 1;
		atomVector[i].PDBAtomName = igraph[i];
		atomVector[i].charge      = charge[i];
		atomVector[i].mass        =  amass[i];
		atomVector[i].AMBAtomName = isymbl[i];
		atomVector[i].idximm      = 4;
	}

	for (int i = 0; i < nres - 1; i++)
	{
		for (int j = ipres[i] - 1; j < ipres[i + 1] - 1; j++)
		{
			atomVector[j].PDBResName = lbres[i];
			atomVector[j].PDBResID   = i + 1;
		}
	}

	//cout << ipres[nres - 1] - 1 << '\n';
	//cout << natom - ipres[nres - 1] + 1 << '\n';
	for (int i = 0; i < natom - ipres[nres - 1] + 1; i++)
	{
		int j = i + ipres[nres - 1] - 1;
		atomVector[j].PDBResName = lbres[nres - 1];
		atomVector[j].PDBResID   = nres;
	}
/*
   	     
	合っているような気がするんだけど、うまくいかない
	原因不明　とりあえず使うな >> 未来の自分

	for (int i = 0; i < natom; i++)
	{
		Atom& atom = atomVector[i];

		int iexcl = 0;

		for (int j = 0; j < i; j++)
			iexcl += numex[j];

		for (int j = iexcl; j < iexcl + numex[i]; j++)
		{
			atom.exclusionVector.push_back(inb[j]);
		}

	}
*/
}
