#include <iostream>
#include <iomanip>
#include <ctime>
#include <cmath>
#include <fstream>
#include <sstream>
#include "PSF.hpp"
#include "PDBUtil.hpp"

using namespace std;

PSF::PSF(string filename) 
{
	ptr_atomVector = &atomVector;
	this -> filename = filename;
	read_PSF();
}

PSF::PSF(string filename, vector<Atom>* ptr_atomVector) 
{
	this -> ptr_atomVector = ptr_atomVector;
	this -> filename = filename;
	read_PSF();
}

void PSF::read_PSF()
{
	ifstream fi(filename.c_str());

	if ( !fi )
	{
		cerr << "エラー：ファイル" << filename << "を開けません。\n";

		return;
	}

	int loc, ncnt;
	string s;
	
	while ( getline(fi, s) )
	{
		loc = s.find("NATOM", 0);
		if (loc != string::npos) 
		{
			istringstream is(s.substr(0, 8));
			is >> natom;
			get_atom_list(natom, fi); 
		}

		loc = s.find("NBOND", 0);
		if (loc != string::npos)
		{
			istringstream is(s.substr(0, 8));
			is >> ncnt;
			get_bond_list(ncnt, fi);
		}

		loc = s.find("NTHETA", 0);
		if (loc != string::npos)
		{
			istringstream is(s.substr(0, 8));
			is >> ncnt;
			get_angle_list(ncnt, fi);
		}

		loc = s.find("NPHI", 0);
		if (loc != string::npos)
		{
			istringstream is(s.substr(0, 8));
			is >> ncnt;
			get_dihedral_list(ncnt, fi);
		}

		loc = s.find("NIMPHI", 0);
		if (loc != string::npos)
		{
			istringstream is(s.substr(0, 8));
			is >> ncnt;
			get_improper_list(ncnt, fi);
		}

		loc = s.find("NCRTERM", 0);
		if (loc != string::npos)
		{
			istringstream is(s.substr(0, 8));
			is >> ncnt;
			get_cmap_list(ncnt, fi);
		}
	}

	fi.close();
}

void PSF::get_atom_list(const int natom, ifstream& fi)
{
	cout << "REMARK TOTAL ATOM : " << natom << '\n';
	ptr_atomVector -> resize(natom);
	string s;
	this->nwater = 0;
	for (int i = 0; i < natom; i++)
	{
		Atom& atom = ptr_atomVector->at(i);
		getline(fi, s);
		istringstream is(s);
		is >> atom.PSFIndex 
		   >> atom.PSFSegmentName
		   >> atom.PSFResID
		   >> atom.PSFResName
		   
		   /*  overwrite PDB value */
		   >> atom.PDBAtomName
		   
		   >> atom.PSFAtomName 
		   >> atom.charge
		   >> atom.mass;

		atom.invmass = 1. / atom.mass;

		if (atom.PDBAtomName == "OH2") ++nwater;
	}

	cout << "REMARK Number of water molecules: " << nwater << '\n';
}

void PSF::get_bond_list(const int nbond, ifstream& fi)
{
	cout << "REMARK NUMBER OF BOND FROM PSF : " << nbond << '\n';

	string s;

	while ( getline(fi, s) && s != "" )
	{
		istringstream is(s);
		int b;
		while (is >> b)
		{
			bondArray.push_back(b);
		}
	}
	if (bondArray.size() / 2 != nbond)
	{
		die("somethig is wrong with PSF::get_bond_list()");
	}
/*
	// Debug      list all bond

	for (int i = 0; i < nbond; i++)
	{
		cout << "REMARK  "
		     << "bond no. " << i+1 << " / " << nbond << " "
		     << bondArray[i*2] << " - " << bondArray[i*2 + 1] << " "
		     << ptr_atomVector -> at(bondArray[i*2] - 1).PSFAtomName
		     << " - "
		     << ptr_atomVector -> at(bondArray[i*2 + 1] - 1).PSFAtomName << '\n';
		if (ptr_atomVector -> at(bondArray[i * 2] - 1).PSFSegmentName == "WAT")
			break;
		num_nwt_bond++;

	}
	cout << "REMARK NUMBER OF BOND IN NON-WAT : " << num_nwt_bond << '\n';
*/
}  /*  end of func() get_bond_list  */

bool PSF::set_bond_parm(vector<Bond>& bondParmVector)
{
	string iat1;
	string iat2;
	string jat1;
	string jat2;

	// now make bondVector
	for (int i = 0; i < bondArray.size() / 2; i++)
	{
		iat1 = ptr_atomVector -> at(bondArray[i * 2    ] - 1).PSFAtomName;
		iat2 = ptr_atomVector -> at(bondArray[i * 2 + 1] - 1).PSFAtomName;
		for (int j = 0; j < bondParmVector.size(); j++)
		{
			jat1 = bondParmVector[j].at1;
			jat2 = bondParmVector[j].at2;

			if ((iat1 == jat1 && iat2 == jat2) || (iat1 == jat2 && iat2 == jat1))
			{
				Bond bond(
				&ptr_atomVector -> at(bondArray[i * 2] - 1),
				&ptr_atomVector -> at(bondArray[i * 2 + 1] - 1),
				bondParmVector[j].Kb,
				bondParmVector[j].b0);

				// set flag
				bond.set_flag_XH();

				bondVector.push_back(bond);

				break;
			}

			if (j == bondParmVector.size() - 1)
			{
				cerr
				<< "error: bond parameter missing for:\n"
				<< iat1 << " - " << iat2 << endl;
				return false;
			}
		}
	}

	// remove all elements in bondArray
	bondArray.clear();
	// remove all parameters in BondParmArray
	bondParmVector.clear();

	cout << "REMARK Bond Parameters are successfully assigned.\n";
	return true;

} // end of func() set_bond_parm

void PSF::get_angle_list(const int nangle, ifstream& fi)
{
	cout << "REMARK NUMBER OF ANGLE FROM PSF : " << nangle << '\n';

	string s;

	while ( getline(fi, s) && s != "" )
	{
		istringstream is(s);
		int b;
		while (is >> b)
		{
			angleArray.push_back(b);
		}
	}

	/*  Debug      list all angle  */
/*
	for (int i = 0; i < nangle; i++)
	{
		cout << "REMARK  "
		     << "angle no. " << i+1 << " / " << nangle << " "
		     << angleArray[i*3] << " - "
		     << angleArray[i*3 + 1] << " - "
		     << angleArray[i*3 + 2] <<  " "
		     << ptr_atomVector -> at(angleArray[i*3] - 1).PDBAtomName
		     << " - "
		     << ptr_atomVector -> at(angleArray[i*3 + 1] - 1).PDBAtomName
		     << " - "
		     << ptr_atomVector -> at(angleArray[i*3 + 2] - 1).PDBAtomName << '\n';
	}
*/

} /* end of func() get_angle_list() */

bool PSF::set_angle_parm(vector<Angle>& angleParmVector)
{
	string iat1;
	string iat2;
	string iat3;
	string jat1;
	string jat2;
	string jat3;

	// make AngleVector from parameters
	for (int i = 0; i < angleArray.size() / 3; i++)
	{
		iat1 = ptr_atomVector -> at(angleArray[i * 3    ] - 1).PSFAtomName;
		iat2 = ptr_atomVector -> at(angleArray[i * 3 + 1] - 1).PSFAtomName;
		iat3 = ptr_atomVector -> at(angleArray[i * 3 + 2] - 1).PSFAtomName;
		for (int j = 0; j < angleParmVector.size(); j++)
		{
			jat1 = angleParmVector[j].at1;
			jat2 = angleParmVector[j].at2;
			jat3 = angleParmVector[j].at3;

			if ((iat2 == jat2) &&
			((iat1 == jat1 && iat3 == jat3) || (iat1 == jat3 && iat3 == jat1)))
			{
				Angle angle(
				&ptr_atomVector->at(angleArray[i * 3] - 1),
				&ptr_atomVector->at(angleArray[i * 3 + 1] - 1),
				&ptr_atomVector->at(angleArray[i * 3 + 2] - 1),
				angleParmVector[j].Kt,
				angleParmVector[j].t0,
				angleParmVector[j].Kub,
				angleParmVector[j].s0);

				// set flag
				//angle.set_flag_XH();

				angleVector.push_back(angle);

				break;
			}

			if (j == angleParmVector.size() - 1)
			{
				cerr << "error: angle parameter missing for: "
					<< iat1 << " - " << iat2
					<< " - " << iat3 << '\n';
				return false;
			}
		}
	}
/*
	cout << "DEBUG angle " << angleArray.size() / 3 << " should be "
	<< AngleVector.size() << '\n';
*/
	// remove all elements in angleArray
	angleArray.clear();
	// remove all elements in AngleParmVector
	angleParmVector.clear();

	cout << "REMARK Angle Parameters are successfully assigned.\n";
	return true;

} // end of func() set_angle_parm

void PSF::get_dihedral_list(const int ndihedral, ifstream& fi)
{
	cout
	<< "REMARK NUMBER OF UNIQUE DIHEDRAL FROM PSF : " << ndihedral << '\n';

	string s;
	int b;
	while ( getline(fi, s) && s != "" )
	{
		istringstream is(s);
		while (is >> b)
			dihedralArray.push_back(b);
	}

	/*  Debug      list all dihedral  */
/*
	for (int i = 0; i < ndihedral; i++)
	{
		cout << "REMARK  "
		     << "dihedral no. " << i+1 << " / " << ndihedral << " "
		     << dihedralArray[i*4    ] << " - "
		     << dihedralArray[i*4 + 1] << " - "
		     << dihedralArray[i*4 + 2] << " - "
		     << dihedralArray[i*4 + 3] << " "
		     << ptr_atomVector -> at(dihedralArray[i*4] - 1).PSFAtomName
		     << " - "
		     << ptr_atomVector -> at(dihedralArray[i*4 + 1] - 1).PSFAtomName
		     << " - "
		     << ptr_atomVector -> at(dihedralArray[i*4 + 2] - 1).PSFAtomName
		     << " - "
		     << ptr_atomVector -> at(dihedralArray[i*4 + 3] - 1).PSFAtomName << '\n';
	}*/

}

bool PSF::set_dihedral_parm(vector<Dihedral>& dihedralParmVector)
{
	string iat1, iat2, iat3, iat4;
	string jat1, jat2, jat3, jat4;

	// now make dihedralVector
	for (int i = 0; i < dihedralArray.size() / 4; i++)
	{
		bool assigned = false;

		int i4 = i * 4;

		int id1 = dihedralArray[i4    ] - 1;
		int id2 = dihedralArray[i4 + 1] - 1;
		int id3 = dihedralArray[i4 + 2] - 1;
		int id4 = dihedralArray[i4 + 3] - 1;

		iat1 = ptr_atomVector -> at( id1 ).PSFAtomName;
		iat2 = ptr_atomVector -> at( id2 ).PSFAtomName;
		iat3 = ptr_atomVector -> at( id3 ).PSFAtomName;
		iat4 = ptr_atomVector -> at( id4 ).PSFAtomName;

		for (Dihedral& dparm: dihedralParmVector)
		{
			// skip wildcard parameters
			if (dparm.is_wildcard) continue;

			jat1 = dparm.at1; jat2 = dparm.at2;
			jat3 = dparm.at3; jat4 = dparm.at4;

			if ((iat1 == jat1 && iat2 == jat2 && iat3 == jat3 && iat4 == jat4) || (iat1 == jat4 && iat2 == jat3 && iat3 == jat2 && iat4 == jat1))
			{
				assigned = true;

				Dihedral dihed(
				&ptr_atomVector -> at( id1 ),
				&ptr_atomVector -> at( id2 ),
				&ptr_atomVector -> at( id3 ),
				&ptr_atomVector -> at( id4 ),
				dparm.Kchi,
				dparm.n,
				dparm.delta);

				dihedralVector.push_back(dihed);

				// ここでbreakすると
			        // O    C    N    CP1 みたいにパラメータを
				// 複数持つ二面角にパラメータがアサインされない.
				//break;
			}
		}

		// forward to next dihedral if assigned
		if (assigned) continue;

		// paramaeres for this dihedral is missing.
		// search for wildcard parameters
		for (Dihedral& dparm: dihedralParmVector)
		{
			// skip non-wildcard parameters
			if (!dparm.is_wildcard) continue;

			jat1 = dparm.at1; jat2 = dparm.at2;
			jat3 = dparm.at3; jat4 = dparm.at4;

			// assumed X-A-B-X, not A-X-X-B... is this OK?
			if ( (iat2 == jat2 && iat3 == jat3) ||
					(iat2 == jat3 && iat3 == jat2))
			{
				assigned = true;

				Dihedral dihed(
				&ptr_atomVector -> at( id1 ),
				&ptr_atomVector -> at( id2 ),
				&ptr_atomVector -> at( id3 ),
				&ptr_atomVector -> at( id4 ),
				dparm.Kchi,
				dparm.n,
				dparm.delta);

				dihedralVector.push_back(dihed);
			}
		}

		if (!assigned)
		{
			cerr << "\nerror: Dihedral parameter missing for: "
					<< iat1 << " - " << iat2 << " - "
					<< iat3 << " - " << iat4 << "\n\n";
			return false;
		}
	} // for dihedralArray
/*
	// DEBUG: list all dihedral
	for (int i = 0; i < DihedralVector.size(); i++)
	{
		cout << "DEBUG " << i + 1 << " / " << DihedralVector.size()
		<< " "
		<< DihedralVector[i].ptr_atom1 -> PSFAtomName << " - "
		<< DihedralVector[i].ptr_atom2 -> PSFAtomName << " - "
		<< DihedralVector[i].ptr_atom3 -> PSFAtomName << " - "
		<< DihedralVector[i].ptr_atom4 -> PSFAtomName << "   "
		<< DihedralVector[i].Kchi << "  "
		<< DihedralVector[i].n    << "  "
		<< DihedralVector[i].delta << endl;
	}

	cout << "DEBUG dihedral " << dihedralArray.size() / 4 << " should be "
	<< DihedralVector.size() << '\n';
*/
	cout << "REMARK Dihedral parameters successfully assigned.\n";
	cout << "REMARK NUMBER OF DUPLICATED DIHEDRAL = "
		<< dihedralVector.size() << '\n';

	dihedralArray.clear();
	dihedralParmVector.clear();

	return true;

} // end of func() set_dihedral_parm


void PSF::get_improper_list(int nimproper, ifstream& fi)
{
	cout << "REMARK NUMBER OF IMPROPER FROM PSF : " << nimproper << '\n';

	string s;
	int b;
	while ( getline(fi, s) && s != "" )
	{
		istringstream is(s);
		while (is >> b)
			improperArray.push_back(b);
	}
/*
	for (int i = 0; i < nimproper; i++)
	{
		cout << "REMARK  "
		     << "improper no. " << i+1 << " / " << nimproper << " "
		     << improperArray[i*4    ] << " - "
		     << improperArray[i*4 + 1] << " - "
		     << improperArray[i*4 + 2] << " - "
		     << improperArray[i*4 + 3] << " "
		     << ptr_atomVector -> at(improperArray[i*4] - 1).PDBAtomName
		     << " - "
		     << ptr_atomVector -> at(improperArray[i*4 + 1] - 1).PDBAtomName
		     << " - "
		     << ptr_atomVector -> at(improperArray[i*4 + 2] - 1).PDBAtomName
		     << " - "
		     << ptr_atomVector -> at(improperArray[i*4 + 3] - 1).PDBAtomName << '\n';
	}
*/
} // end of func() get_improper_list

bool PSF::set_improper_parm(vector<Improper>& improperParmVector)
{
	for (int i = 0; i < improperArray.size() / 4; i++)
	{
		bool assigned = false;

		int id1 = improperArray[i * 4    ];
		int id2 = improperArray[i * 4 + 1];
		int id3 = improperArray[i * 4 + 2];
		int id4 = improperArray[i * 4 + 3];

		string iat1 = ptr_atomVector -> at(id1 - 1).PSFAtomName;
		string iat2 = ptr_atomVector -> at(id2 - 1).PSFAtomName;
		string iat3 = ptr_atomVector -> at(id3 - 1).PSFAtomName;
		string iat4 = ptr_atomVector -> at(id4 - 1).PSFAtomName;

		for (Improper& iparm: improperParmVector)
		{
			// skip wildcard parameters
			if (iparm.is_wildcard) continue;

			string jat1 = iparm.at1;
			string jat2 = iparm.at2;
			string jat3 = iparm.at3;
			string jat4 = iparm.at4;

			if ((iat1 == jat1 && iat2 == jat2 && iat3 == jat3 && iat4 == jat4) || (iat1 == jat4 && iat2 == jat3 && iat3 == jat2 && iat4 == jat1))
			{
				assigned = true;

				Improper impro(
				&ptr_atomVector -> at(id1-1),
				&ptr_atomVector -> at(id2-1),
				&ptr_atomVector -> at(id3-1),
				&ptr_atomVector -> at(id4-1),
				iparm.Kpsi,
				iparm.psi0);

				improperVector.push_back(impro);

				break;
			}
		}

		if (assigned) continue;

		// then try wildcard parameters
		for (Improper& iparm: improperParmVector)
		{
			// skip non-wildcard parameters
			if (!iparm.is_wildcard) continue;

			string jat1 = iparm.at1;
			string jat2 = iparm.at2;
			string jat3 = iparm.at3;
			string jat4 = iparm.at4;

			if ((iat1 == jat1 && iat4 == jat4) || (iat1 == jat4 && iat4 == jat1))
			{
				assigned = true;

				Improper impro(
				&ptr_atomVector -> at(id1-1),
				&ptr_atomVector -> at(id2-1),
				&ptr_atomVector -> at(id3-1),
				&ptr_atomVector -> at(id4-1),
				iparm.Kpsi,
				iparm.psi0);

				improperVector.push_back(impro);

				break;
			}
		}

		if (!assigned)
		{
			cerr << "error: improper parameter missing for: "
				<< iat1 << " - " << iat2 << " - "
				<< iat3 << " - " << iat4 << '\n';
			return false;
		}
	}

	improperArray.clear();
	improperParmVector.clear();
	cout << "REMARK Improper parameters successfully assigned.\n";
	return true;
/*
	cout << "DEBUG IMPROPER ARRAY.size() / 4 = " << improperArray.size() / 4
	<< " should be \"assigned\" " << assigned << '\n';
	for (int i = 0; i < ImproperVector.size(); i++)
	{
		Improper* impro = &ImproperVector[i];
		cout << i + 1 << "   "
		<< impro->ptr_atom1->PSFIndex << " "
		<< impro->ptr_atom1->PSFAtomName << " "
		<< impro->ptr_atom2->PSFIndex << " "
		<< impro->ptr_atom2->PSFAtomName << " "
		<< impro->ptr_atom3->PSFIndex << " "
		<< impro->ptr_atom3->PSFAtomName << " "
		<< impro->ptr_atom4->PSFIndex << " "
		<< impro->ptr_atom4->PSFAtomName << " "
		<< impro->Kpsi << " "
		<< impro->psi0 << endl;
	}
*/
}// end of func() set_improper_parm

void PSF::get_cmap_list(int ncmap, ifstream& fi)
{
	cout << "REMARK NUMBER OF CROSS-TERM : "  << ncmap << '\n';
	string s;
	while(getline(fi, s) && !s.empty())
	{
		istringstream is(s);
		int itmp;
		while (is >> itmp) cmapArray.push_back(itmp);
	}
/*	cout << "cmapArray.size() / 4 = " << cmapArray.size() / 4 << endl;
	cout << "DEBUG>>>>>>>>>>>>>>>>>>>>>>>\n";
	for (int i = 0; i < ncrossterm * 8; i = i + 4)
	{
		//showAtom(&ptr_atomVector->at(cmapArray[i] - 1));
		int i1 = cmapArray[i    ]; int i2 = cmapArray[i + 1];
		int i3 = cmapArray[i + 2]; int i4 = cmapArray[i + 3];
		cout << i1 << " " << i2 << " " << i3 << " " << i4 << '\n';
	}*/
}

void PSF::set_cmap_parm(vector<Cmap>& cmapParmVector)
{
	for (int i = 0; i < cmapArray.size() / 8; i++)
	{
		vector<Atom*>  p_at(8);
		vector<string> attype(8);
		for (int j = 0; j < 8; j++)
		{
			p_at[j] =
			&ptr_atomVector->at(cmapArray[8 * i + j] - 1);

			attype[j] = p_at[j]->PSFAtomName;
		}

		Dihedral dphi(p_at[0], p_at[1], p_at[2], p_at[3], 0., 0., 0.);
		Dihedral dpsi(p_at[4], p_at[5], p_at[6], p_at[7], 0., 0., 0.);

		Cmap cmap(dphi, dpsi);

		for (int j = 0; j < cmapParmVector.size(); j++)
		{
			bool match = true;
			for (int k = 0; k < 8; k++)
			{
				if (cmapParmVector[j].attype[k]
					!= attype[k])
				{
					match = false;
					break;
				}
			}
			if (match)
			{
				cmap.attype = cmapParmVector[j].attype;

				cmap.cmap_type_index = j;

				break;
			}
		}

		cmapVector.push_back(cmap);
	}

	cout << "REMARK " << cmapVector.size() << " cmap created.\n";
}

bool PSF::set_lj_parm(vector<Atom>& LJParmVector)
{
	if (!LJParmVector.size())
	{
		cerr << "error: please set LJ parm\n";
		return false;
	}

	for (int i = 0; i < ptr_atomVector -> size(); i++)
	{
		Atom& atom = ptr_atomVector -> at(i);
		for (int j = 0; j < LJParmVector.size(); j++)
		{
			Atom& LJatom = LJParmVector.at(j);
			//cout << atom->PSFAtomName << "  "
			//	<< LJatom->PSFAtomName << '\n';
			if (atom.PSFAtomName == LJatom.PSFAtomName)
			{
				atom.epsilon = LJatom.epsilon;
				atom.Rmin_div2 = LJatom.Rmin_div2;
				atom.eps1_4 = LJatom.eps1_4;
				atom.Rmin1_4 = LJatom.Rmin1_4;

				break;
			}
			if (j == LJParmVector.size() - 1)
			{
				cerr << "error: LJ parameter missing for: "
					<< atom.PSFAtomName << '\n';
				return false;
			}
		}
	}

	// remove all elements
	LJParmVector.clear();
	cout << "REMARK LJ parameters are successfully assigned.\n";
	return true;

} // end of function() set_LJ_parm

void PSF::make_exclusion_vector()
{
	if (bondVector.empty() || angleVector.empty())
	{
		cerr << "error in PSF::make_exclusion_vector()\n"
		     << "please set parameter of bond and angle\n";
		return;
	}

	// exclusion of 1-2 pair
	for (auto& bond: bondVector)
	{
		if (bond.Kb < 0.1) continue;

		Atom* atom1 = bond.ptr_atom1;
		Atom* atom2 = bond.ptr_atom2;

		atom1 -> exclusionVector.push_back(atom2 -> PSFIndex);
		atom2 -> exclusionVector.push_back(atom1 -> PSFIndex);
	}

	// exclusion of 1-3 pair
	for (auto& angle: angleVector)
	{
		Atom* atom1 = angle.ptr_atom1;
		Atom* atom3 = angle.ptr_atom3;

		atom1 -> exclusionVector.push_back(atom3 -> PSFIndex);
		atom3 -> exclusionVector.push_back(atom1 -> PSFIndex);
	}

} // end of function PSF::make_exclusion_vector()

void PSF::make_scaled1_4_vector()
{
	for (auto& dihed: dihedralVector)
	{
		Atom* atom1 = dihed.ptr_atom1;
		Atom* atom4 = dihed.ptr_atom4;

		atom1 -> scaled1_4Vector.push_back(atom4 -> PSFIndex);
		atom4 -> scaled1_4Vector.push_back(atom1 -> PSFIndex);
	}
}

void PSF::writePDB(string filename, string header, const vector<int>& indexVector)
{
	PDBUtil PDBJOBS;
	
	ofstream fo(filename.c_str());

	time_t now = time(NULL);
	struct tm* local = localtime(&now);
	
	fo << "REMARK " << asctime(local);
	fo << "REMARK " << header;
	fo << '\n';
	
	for (int i = 0; i < indexVector.size(); i++)
	{
		int ii = indexVector[i] - 1;
		PDBJOBS.writePDBline(fo, ptr_atomVector -> at(ii));
	}
}// end of function writePDB()

void PSF::writePDB(string filename, string header)
{
	PDBUtil PDBJOBS;
	
	ofstream fo(filename.c_str());

	time_t now = time(NULL);
	struct tm* local = localtime(&now);
	
	fo << "REMARK " << asctime(local);
	fo << "REMARK " << header;
	fo << '\n';
	
	for (int i = 0; i < ptr_atomVector -> size(); i++)
	{
		PDBJOBS.writePDBline(fo, ptr_atomVector -> at(i));
	}
}// end of function writePDB()

void PSF::writePDB(string filename, vector<Atom>& atomVector, string header)
{
	PDBUtil PDBJOBS;
	
	ofstream fo(filename.c_str());

	time_t now = time(NULL);
	struct tm* local = localtime(&now);
	
	fo << "REMARK  " << asctime(local);
	fo << "REMARK  " << header;
	fo << '\n';
	
	for (int i = 0; i < atomVector.size(); i++)
	{
		Atom& base = ptr_atomVector->at(i);
		Atom& atom = atomVector[i];

		atom.PSFIndex = base.PSFIndex;
		atom.PDBAtomName = base.PDBAtomName;
		atom.PSFResName = base.PSFResName;
		atom.PSFSegmentName = base.PSFSegmentName;
		atom.PSFResID = base.PSFResID;
//		atom.PDBo = base.PDBo;
//		atom.PDBb = base.PDBb;
		atom.PDBo = 0;
		atom.PDBb = 0;

		PDBJOBS.writePDBline(fo, atomVector[i]);
	}
}// end of function writePDB()

void PSF::showAtom(const Atom& atom)
{
	cout
	<< atom.PSFIndex << ' ' << atom.PDBAtomName
	<< " @" << atom.PSFResName << atom.PSFResID 
	<< "  " << atom.charge << 'e' 
	<< "  " << atom.mass   << "amu\n";
}

void PSF::showResidue(int resID)
{
	bool hit = false;
	for (int i = 0; i < ptr_atomVector -> size(); i++)
	{
		if (ptr_atomVector -> at(i).PSFResID == resID)
		{
			cout << ptr_atomVector -> at(i).PSFResName << '\n';
			hit = true;
			break;
		}
	}
	if (!hit)
		cout << "REMARK warning: no residue hit for resID: "
		<< resID << '\n';
}

void PSF::make_water_shake_bond_array()
{
	for (int i = 0; i < ptr_atomVector->size(); i++)
	{
		Atom& atom = ptr_atomVector->at(i);

		if (atom.PDBAtomName == "OH2")
		{
			// store 1-based index
			bondArray.push_back(i + 1 + 1);
			bondArray.push_back(i + 2 + 1);
		}
	}
}

Eigen::Vector3d PSF::getcomof_zerobased(const vector<int>& ind)
{
	Eigen::Vector3d com(0., 0., 0.);

	for (auto& i: ind)
		com += ptr_atomVector->at(i).position;
	com /= ind.size();

	return com;
}
