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
	bool skipCmap;
	string s;

	cout << "REMARK ---------------------\n"
		"REMARK parameter summary (" << filename << ")\n";

	while (getline(fi, s))
	{
		if (s.substr(0, 5) == "BONDS")
		{	
			get_bond_parameters(fi);
			get_angle_parameters(fi);
			get_dihedral_parameters(fi);
			skipCmap = get_improper_parameters(fi);
			if (!skipCmap) get_cmap_parameters(fi);
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

void LoadParm::get_dihedral_parameters(ifstream& fi)
{
	string s;
	while (getline(fi, s))
	{
		if (s.substr(0, 5) == "IMPRO")
			break;

		istringstream is(s);
		string at1, at2, at3, at4;
		double Kchi;
		double n;
		double delta;

		is >> at1 >> at2 >> at3 >> at4 >> Kchi >> n >> delta;
		if (!is) continue;

		Dihedral dihed(at1, at2, at3, at4, Kchi, n, delta);

		dihed.is_wildcard = false;
		if (at1 == "X" || at2 == "X" || at3 == "X" || at4 == "X")
			dihed.is_wildcard = true;

		dihedralParmVector.push_back(dihed);

	}

	cout << "REMARK " << dihedralParmVector.size() << " dihedrals\n";

	// multiplicity check
	/*int multiplicity = 1;
	string sbef = DihedralParmVector[0].at1
			+ DihedralParmVector[0].at2
			+ DihedralParmVector[0].at3
			+ DihedralParmVector[0].at4;

	for (int i = 1; i < DihedralParmVector.size(); i++)
	{
		string snow = DihedralParmVector[i].at1
				+ DihedralParmVector[i].at2
				+ DihedralParmVector[i].at3
				+ DihedralParmVector[i].at4;
		if (snow == sbef)
			++multiplicity;
		else
		{
			DihedralParmVector[i - 1].multiplicity = multiplicity;
			multiplicity = 1;
		}

		// reserve previous values
		sbef = snow;
	}

	for (int i = 0; i < DihedralParmVector.size(); i++)
		cout << "DEBUG " << i + 1 << " " << DihedralParmVector[i].multiplicity << '\n';*/
}

bool LoadParm::get_improper_parameters(ifstream& fi)
{
	string at1, at2, at3, at4, s;
	double Kpsi, dtmp, psi0;

	while (getline(fi, s))
	{
		if (s.substr(0, 4) == "CMAP" || s.substr(0, 9) == "NONBONDED")
			break;

		istringstream is(s);
		is >> at1 >> at2 >> at3 >> at4 >> Kpsi >> dtmp >> psi0;
		if (!is) continue;

		Improper impro(at1, at2, at3, at4, Kpsi, psi0);

		impro.is_wildcard = false;
		if (at1 == "X" || at2 == "X" || at3 == "X" || at4 == "X")
			impro.is_wildcard = true;

		improperParmVector.push_back(impro);
	}

	cout << "REMARK " << improperParmVector.size() << " impropers\n";
/*
	for (int i = 0; i < ImproperParmVector.size(); i++)
		cout << i + 1 << " "
		<< ImproperParmVector[i].at1 << " "
		<< ImproperParmVector[i].at2 << " "
		<< ImproperParmVector[i].at3 << " "
		<< ImproperParmVector[i].at4 << " "
		<< ImproperParmVector[i].Kpsi << " "
		<< ImproperParmVector[i].psi0 << endl;
*/
	return (s.substr(0, 4) == "CMAP") ? false : true;
}


void LoadParm::get_cmap_parameters(ifstream& fi)
{
	string s;
	while (getline(fi, s))
	{
		if (s.substr(0, 1) == "!" || s.empty()) continue;
		if (s.find("NONBONDED", 0) != string::npos) break;

		vector<string> vs(8);
		int resolution = 0;;
		istringstream is1(s);
		is1
			>> vs[0] >> vs[1] >> vs[2] >> vs[3] >> vs[4]
			>> vs[5] >> vs[6] >> vs[7] >> resolution;

		if (resolution != 24)
		{
			cerr << "\n\n  ERROR ERROR ERROR ERROR \n\n\n"
				"\n\n  Unexpected Cmap Resolution of "
				<< resolution << "    \n\n\n\n";
			die("program terminated.");
		}

		Cmap cmap(vs[0], vs[1], vs[2], vs[3],
			  vs[4], vs[5], vs[6], vs[7]);
		int icnt = 0;

		for (int i = 0; i < resolution; i++)
		{
			getline(fi, s); // brank line
			getline(fi, s); // line like "! phi = -180.0"
			//cout << s << endl;
			for (int j = 0; j < resolution; j++)
			{
				fi >> cmap.cmap_grid_data( j, i );
			}
			getline(fi, s); // need this.
		}

		cmapParmVector.push_back(cmap);
	}

	cout << "REMARK " << cmapParmVector.size() << " cmaps\n";
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

	return true;
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
