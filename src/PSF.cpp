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
