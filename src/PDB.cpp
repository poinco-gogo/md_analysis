#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include "PDB.hpp"
#include "PDBUtil.hpp"

using namespace std;

PDB::PDB(string filename)
{
	this -> filename = filename;
}

void PDB::read_PDB(vector<Atom>& atomVector)
{
	ifstream fi(filename.c_str());

	if ( !fi )
	{
		cerr << "error: file " << filename << " not exists.\n";
		return;
	}

	if (atomVector.size() == 0)
	{
		cerr << "PDB error: topology not loaded?\n";
		return;
	}

	string s;
	int icnt = 0;
	
	while (getline(fi, s))
	{
		if (s.substr(0, 4) != "ATOM")
			continue;

		// i am "ATOM"
		++icnt;
		
		if (icnt > atomVector.size())
		{	
			cout << "REMARK PDB contains water?\n";
			break;
		}
		
		// read coordinates and PDBo, PDBb.
		istringstream is(s.substr(30, 36));
		Atom& at = atomVector[icnt - 1];
		Eigen::Vector3d& pv = at.position;
		is
		>> pv.x() 
		>> pv.y() 
		>> pv.z()
		>> at.PDBo
		>> at.PDBb;
	}
}

bool PDB::LoadCoords(vector<Atom>& atomVector)
{
	if (!atomVector.size())
	{
		cerr << "error: psf not loaded?\n";
		return false;
	}
	ifstream fi(filename.c_str());
	if ( !fi )
	{
		cerr << "error: file " << filename << " not exists.\n";
		return false;
	}
	string s;
	int icnt = 0;
	while(getline(fi, s))
	{
		if (s.substr(0, 4) != "ATOM")
			continue;
		else if (atomVector.size() < icnt + 1)
		{
			cerr << "warning: PDB contains water? ... anyway, "
				<< icnt << " coordinates loaded from PDB.\n";
			return true;
		}

		// i`m ATOM.
		istringstream is(s.substr(30, 24));
		Atom& atom = atomVector[icnt++];
		is
		>> atom.position.x()
		>> atom.position.y()
		>> atom.position.z();
	}
	if (icnt != atomVector.size())
	{
		cerr << "error: PDB has wrong atom size.\n";
		return false;
	}
	else
	{
		cout 
		<< "REMARK " << icnt << " coordinates loaded from PDB.\n";
	}
	return true;
}

bool PDB::LoadNotOnlyCoords(vector<Atom>& atomVector)
{
	if (!atomVector.size())
	{
		cerr << "error: psf not loaded?\n";
		return false;
	}
	ifstream fi(filename.c_str());
	if ( !fi )
	{
		cerr << "error: file " << filename << " not exists.\n";
		return false;
	}

	PDBUtil JOB;
	string s;
	int icnt = 0;
	while(getline(fi, s))
	{
		if (s.substr(0, 4) != "ATOM")
		{
			continue;
		}

		JOB.readPDBline(s, atomVector[icnt++]);
	}

	if (icnt != atomVector.size())
	{
		cerr << "error: PDB has wrong atom size.\n";
		return false;
	}
	else
	{
		cout 
		<< "REMARK " << icnt << " coordinates loaded from PDB.\n";
	}
	return true;
}
