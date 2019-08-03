#include <iostream>
#include <fstream>
#include <sstream>
#include "LoadParm.hpp"
using namespace std;

// constructor
LoadParm::LoadParm(string filename)
{
	ifstream fi(filename.c_str());
	if (!fi)
	{
		cerr << "エラー：ファイル" << filename << "は存在しません。\n";
		return;
	}
	
	this -> filename = filename;

	open_fi(fi);
}

void LoadParm::open_fi(ifstream& fi)
{
	string s;

	cout << "REMARK ---------------------\n"
		"REMARK parameter summary (" << filename << ")\n";

	while (getline(fi, s))
	{
		if (s.substr(0, 5) == "BONDS")
		{	
			get_bond_parameters(fi);
			get_angle_parameters(fi);
			get_cmap_parameters(fi);
			get_LJ_parameters(fi);
		}
	}
	cout << "REMARK ---------------------\n";
}

void LoadParm::get_bond_parameters(ifstream& fi)
{
	string s;
	while (getline(fi, s))
	{
		if (s.substr(0, 6) == "ANGLES")
			break;

		istringstream is(s);
		string at1;
		string at2;
		double Kb;
		double b0;
		is >> at1 >> at2 >> Kb >> b0;
		if (!is)
			continue;

		Bond bond(at1, at2, Kb, b0);

		bondParmVector.push_back(bond);
	}

	cout << "REMARK " << bondParmVector.size() << " bonds\n";

	//for (int i = 0; i < BondParmVector.size(); i++)
	//	cout << BondParmVector[i]._at1() << '\n';
}

void LoadParm::get_angle_parameters(ifstream& fi)
{
	string s;
	int iub = 0;
	while (getline(fi, s))
	{
		if (s.substr(0, 5) == "DIHED")
			break;

		istringstream is(s);
		string at1;
		string at2;
		string at3;
		double Kt;
		double t0;
		double Kub;
		double s0;
		is >> at1 >> at2 >> at3 >> Kt >> t0 >> Kub >> s0;
		if (!is)
		{
			istringstream is2(s);
			// probabry not having Urey-Bradley parameters
			is2 >> at1 >> at2 >> at3 >> Kt >> t0;

			if (!is2) continue;

			Kub = 0;
			s0  = 0;
			++iub;
		}
		Angle angle(at1, at2, at3, Kt, t0, Kub, s0);

		angleParmVector.push_back(angle);
	}

	cout << "REMARK " << angleParmVector.size() << " angles (with ub terms)\n";
	cout << "REMARK " << iub << " angles (wo ub terms)\n";

/*	itmp = 0;
	for (int i = 0; i < AngleParmVector.size(); i++)
	{
		cout << AngleParmVector[i]._at1() << " - ";
		cout << AngleParmVector[i]._at2() << " - ";
		cout << AngleParmVector[i]._at3() << "   ";
		if (AngleParmVector[i]._Kub() > 0.0)
			cout << "1-3ub pair  no. " << ++itmp;
		cout << '\n';
	}*/
}

void LoadParm::get_cmap_parameters(ifstream& fi)
{
	string s;
	while (getline(fi, s))
	{
		if (s.substr(0, 1) == "!" || s.empty()) continue;
		if (s.find("NONBONDED", 0) != string::npos) break;
	}
}

void LoadParm::get_LJ_parameters(ifstream& fi)
{
	string s;
	while(getline(fi, s))
	{
		istringstream is(s);
		
		string atomName;
		double ignored;
		double epsilon;
		double Rmin_div2;
		double eps1_4;
		double Rmin1_4;

		is >> atomName >> ignored >> epsilon >> Rmin_div2
		>> ignored >> eps1_4 >> Rmin1_4;
		
		if (!is)
		{
			// probably not having 1-4 parameter...
			istringstream is2(s);
			is2 >> atomName >> ignored >> epsilon >> Rmin_div2;

			if (!is2) continue;
			
			// set normal value
			eps1_4  = epsilon;
			Rmin1_4 = Rmin_div2;
		}

		// assign parameters
		Atom atom;
		atom.PSFAtomName = atomName;
		atom.epsilon     = epsilon;
		atom.Rmin_div2   = Rmin_div2;
		atom.eps1_4      = eps1_4;
		atom.Rmin1_4     = Rmin1_4;
		LJParmVector.push_back(atom);
	}
	cout << "REMARK " << LJParmVector.size() << " lj parms\n";
}

bool LoadParm::merge(string filename)
{
	ifstream fi(filename.c_str());

	if (!fi)
	{
		cerr << "error: file not exist.\n";
		return false;
	}

	open_fi(fi);

	duplication_check();
}

void LoadParm::duplication_check()
{
	for (int i = 0; i < angleParmVector.size(); i++)
	{
		Angle& a1 = angleParmVector[i];

		string s1 = a1.at1;
		string s2 = a1.at2;
		string s3 = a1.at3;

		for (int j = i + 1; j < angleParmVector.size(); j++)
		{
			Angle& a2 = angleParmVector[j];

			string t1 = a2.at1;
			string t2 = a2.at2;
			string t3 = a2.at3;

			if (s1 == t1 && s2 == t2 && s3 == t3)
				die("angle same parameters!");
		}
	}
}
