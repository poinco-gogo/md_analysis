#include <sstream>
#include <iostream>
#include "common.hpp"
#include "Atom.hpp"
using namespace std;

Atom::Atom()
{
	reset_parameters();
}

Atom::Atom(string PDBfield)
{
	reset_parameters();
	parse_PDBfield(PDBfield);
}

void Atom::parse_PDBfield(string PDBField)
{
	if (PDBField.substr(0, 4) != "ATOM")
	{
		cerr << "エラー：in <Atom::parse_field>...\n"
		     << "This is not ATOM.\n";
		return;
	}

/*	COLUMNS        DATA  TYPE    FIELD       DEFINITION
	---------------------------------------------------
	 1 -  6        Record name   "ATOM  "
	 7 - 11        Integer       serial      Atom  serial number.
	 13 - 16       Atom          name        Atom name.
	 17            Character     altLoc      Alternate location indicator.
	 18 - 20       Residue name  resName     Residue name.
	 22            Character     chainID     Chain identifier.
	 23 - 26       Integer       resSeq      Residue sequence number.
	 27            AChar         iCode       Code for insertion of residues.
	 31 - 38       Real(8.3)     x           
	 39 - 46       Real(8.3)     y          
	 47 - 54       Real(8.3)     z   
	 55 - 60       Real(6.2)     occupancy   Occupancy.
	 61 - 66       Real(6.2)     tempFactor  Temperature  factor.
	 77 - 78       LString(2)    element     Element symbol,right-justified.
	 79 - 80       LString(2)    charge      Charge  on the atom.
*/
	istringstream is;
	is.str(PDBField.substr(6, 5));
	is >> PDBIndex;
	
	PDBAtomName = PDBField.substr(12, 4);
	is.clear();
	is.str(PDBAtomName);
	is >> PDBAtomName;

	PDBResName = PDBField.substr(17, 3);
	
	is.clear();
	is.str(PDBField.substr(22, 4));
	is >> PDBResID;
	
	is.clear();
	is.str(PDBField.substr(30, 8));
	is >> position.x();
	
	is.clear();
	is.str(PDBField.substr(38, 8));
	is >> position.y();
	
	is.clear();
	is.str(PDBField.substr(46, 8));
	is >> position.z();
}

void Atom::reset_parameters()
{
	position = V3ZERO;
	force    = V3ZERO;
	constfrc = V3ZERO;
	velocity = V3ZERO;
	vnew     = V3ZERO;
	rnew     = V3ZERO;
	rold     = V3ZERO;

	// Langevin parameteres
	A     = 0;
	B     = 0;
	invB  = 0;
	R     = 0;
	gamma = 0;
	
	PDBo  = 0.;
	PDBb  = 0.;

	charge    = 0.;
	mass      = 0.;
	epsilon   = 0.;
	Rmin_div2 = 0.;
	eps1_4    = 0.;
	Rmin1_4   = 0.;

	is_rigid = false;
}

bool Atom::checkExclusionPair(const Atom& atom)
{
	// return true if exclusion pair.
	for (const auto& i: exclusionVector)
        {
		if (i == atom.PSFIndex)
		{
			return true;
		}
	}
	return false;
}

bool Atom::checkScaled1_4Pair(const Atom& atom)
{
	// return true if 1-4 pair.
	for (const auto& i: scaled1_4Vector)
	{
		if (i == atom.PSFIndex)
		{
			return true;
		}
	}
	return false;
}

bool Atom::is_mainchain()
{
	string s = this->PDBAtomName;

	if (s == "CA" || s == "N" || s == "C" || s == "O" || s == "HN"
			|| s == "HA")
		return true;
	else
		return false;
}
