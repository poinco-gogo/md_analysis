#include <Eigen/Core>
#include <iostream>
#include <iomanip>
#include <sstream>
#include "PDBUtil.hpp"
using namespace std;
PDBUtil::PDBUtil()
{
	// currently nothing to do.
}

void PDBUtil::writePDBline(ofstream& fo, Atom& atom)
{
	ostringstream os;
	os << "ATOM" << setw(7) << atom.PSFIndex << ' ';

	string at = atom.PDBAtomName;
	switch (at.size())
	{
		case  1: os << setw(2) << at << "  "; break;
		case  2: os << setw(3) << at <<  ' '; break;
		default: os << setw(4) << at        ; break;
	}
	
	os << ' ' << setw(3) << atom.PSFResName;

	if (atom.PSFResName == "TIP3")
		os << 'W';
	else if (atom.PSFResName == "POPC")
		os << 'L';
	else if (atom.PSFResName == "ALAD")
		os << ' ';
	else
		os << ' ' << atom.PSFSegmentName[0];

	
	os << setw(4) << atom.PSFResID
	<< ' ' << setprecision(3) << fixed
	<< setw(11) << atom.position.x()
	<< setw(8)  << atom.position.y()
	<< setw(8)  << atom.position.z()
	<< setprecision(2)
	<< setw(6) << atom.PDBo
	<< setw(6) << atom.PDBb
	<< setw(9) << atom.PSFSegmentName
	<< '\n';

	string s = os.str();
	fo.write(s.c_str(), s.size());
}

void PDBUtil::readPDBline(string line, Atom& atom)
{
	istringstream is1(line.substr(4, 7)); // index
	istringstream is2(line.substr(11, 6)); // atom name
	istringstream is3(line.substr(17, 4)); // residue name
	istringstream is4(line.substr(21, 1)); // chain identifer
	istringstream is5(line.substr(22, 4)); // residue index
	istringstream is6(line.substr(30, 8)); // x
	istringstream is7(line.substr(38, 8)); // y
	istringstream is8(line.substr(46, 8)); // z
	istringstream is9(line.substr(54, 6)); // occupancy
	istringstream isa(line.substr(60, 6)); // beta
	
	is1 >> atom.PDBIndex;
	is2 >> atom.PDBAtomName;
	is3 >> atom.PDBResName;
	is5 >> atom.PDBResID;
	is6 >> atom.position.x();
	is7 >> atom.position.y();
	is8 >> atom.position.z();
	is9 >> atom.PDBo;
	isa >> atom.PDBb;
}
