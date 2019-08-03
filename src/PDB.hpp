#ifndef ___CLASS_PDB
#define ___CLASS_PDB

#include <vector>
#include <string>
#include "Atom.hpp"

class PDB
{
	private:
	std::string filename;
	
	public:
	// constructor
	PDB(std::string filename);

	std::string _filename() { return filename; }
	bool LoadCoords(std::vector<Atom>& atomVector);
	bool LoadNotOnlyCoords(std::vector<Atom>& atomVector);

	// member function
	private:
	void read_PDB(std::vector<Atom>& atomVector);
};
#endif
