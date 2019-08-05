#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cmath>
#include "ForceUtil.hpp"
using namespace std;
ForceUtil::ForceUtil()
{
	// nothing to do currently...	
}

bool ForceUtil::LoadForce(string filename, vector<Eigen::Vector3d>& v)
{
	ifstream fi(filename.c_str());
	if (!fi)
	{
		cerr << "error: file \"" << filename << "\" not exists.\n";
		return false;
	}

	this -> filename = filename;

	string s;
	while (getline(fi, s))
	{
		if (s.empty() || s.find("REMARK", 0) != string::npos)
			continue;
		Eigen::Vector3d vtmp;
		if (!ReadLine(s, vtmp))
		{
			cerr << "error: bad force data type.\n";
			return false;
		}
		v.push_back(vtmp);
	}
	return true;
}// end of func ForceUtil::LoadForce()

bool ForceUtil::ReadLine(string s, Eigen::Vector3d& vtmp)
{
	istringstream is(s);
	int itmp, icnt;
	string stmp;
	is >> icnt >> itmp >> stmp >> vtmp.x() >> vtmp.y() >> vtmp.z();
	if (!is)
		return false;
	else
	{
		atomIndex.push_back(icnt);
		return true;
	}
}

void ForceUtil::writeForce(string filename, Atom& atom, Eigen::Vector3d& f)
{
	if (filename == "stdout")
		cout << setw(5) << atom.PSFIndex
		<< setw(4) << atom.PSFResID
		<< setw(5) << atom.PDBAtomName
		<< setprecision(8) << scientific
		<< setw(16) << f.x()
		<< setw(16) << f.y()
		<< setw(16) << f.z()
		<< setw(16) << f.norm()
		<< '\n';
}

void ForceUtil::normalize(vector<Eigen::Vector3d>& v)
{
	double denom = 0;
	for (auto& f: v)
		denom += f.squaredNorm();
	denom = 1. / sqrt(denom);

	for (auto& f: v)
		f *= denom;
}
