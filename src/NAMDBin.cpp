#include <iostream>
#include <fstream>
#include "NAMDBin.hpp"
using namespace std;
NAMDBin::NAMDBin(string filename, std::string type)
{
	this->filename = filename;
	this->type     = type;
}

bool NAMDBin::read_fi(vector<Atom>& atomVector)
{
	ifstream fi;
	
	fi.open(this->filename.c_str(), ios::in | ios::binary);

	if (!fi)
	{
		cerr << "\nerror: file " << filename << " does not exist.\n\n";
		return false;
	}

	int itmp;
	fi.read( (char*) &itmp, sizeof(int)); // # of atoms

	if (itmp != atomVector.size())
	{
		cerr << "\nerror: # of atoms specified in the namd bin file ("
		       << itmp << ") does not coincide with that of psf file ("
		       << atomVector.size() << ").\n\n";
		return false;
	}

	if (this->type == "coor")
	{
		read_coor(fi, atomVector);
		return true;
	}
	else if (this->type == "vel")
	{
		read_vel(fi, atomVector);
		return true;
	}
	else
	{
		cerr << "\nerror: unknown namd bin type.\n\n";
		return false;
	}

	return true;
}

void NAMDBin::read_coor(ifstream& fi, vector<Atom>& atomVector)
{
	double dtmp[3];
	for (auto& atom: atomVector)
	{
		fi.read( (char*) &dtmp, sizeof(double) * 3 );

		atom.position.x() = dtmp[0];
		atom.position.y() = dtmp[1];
		atom.position.z() = dtmp[2];
	}
}

void NAMDBin::read_vel(ifstream& fi, vector<Atom>& atomVector)
{
	double dtmp[3];
	for (auto& atom: atomVector)
	{
		fi.read( (char*) &dtmp, sizeof(double) * 3 );

		atom.velocity.x() = dtmp[0];
		atom.velocity.y() = dtmp[1];
		atom.velocity.z() = dtmp[2];
	}
}

bool NAMDBin::write_fo(vector<Atom>& atomVector)
{
	ofstream fo(this->filename.c_str(), ios::out | ios::binary);

	if (!fo)
	{
		cerr << "\nerror: cannot open namd output binary file \n\n";
		return false;
	}

	int itmp = atomVector.size();

	// write # of atoms
	fo.write( (char*) &itmp, sizeof(int) );

	if (this->type == "coor")
	{
		write_coor(fo, atomVector);
		return true;
	}
	else if (this->type == "vel")
	{
		write_vel(fo, atomVector);
		return true;
	}
	else
	{
		cerr << "\nerror: unknown namd bin type\n\n";
		return false;
	}

	return true;
}

void NAMDBin::write_coor(ofstream& fo, vector<Atom>& atomVector)
{
	double dtmp[3];
	for (auto& atom: atomVector)
	{
		dtmp[0] = atom.position.x();
		dtmp[1] = atom.position.y();
		dtmp[2] = atom.position.z();
		fo.write( (char*) dtmp, sizeof(double) * 3 );
	}
}

void NAMDBin::write_vel(ofstream& fo, vector<Atom>& atomVector)
{
	double dtmp[3];
	for (auto& atom: atomVector)
	{
		dtmp[0] = atom.velocity.x();
		dtmp[1] = atom.velocity.y();
		dtmp[2] = atom.velocity.z();
		fo.write( (char*) dtmp, sizeof(double) * 3 );
	}
}
