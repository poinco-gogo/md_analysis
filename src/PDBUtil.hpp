#ifndef ___CLASS_PDBUTIL
#define ___CLASS_PDBUTIL
#include <fstream>
#include <string>
#include "Atom.hpp"
class PDBUtil
{
	public:
	
	// constructor
	PDBUtil();

	// member funtion
	void writePDBline(std::ofstream& fo, Atom& atom);
	void  readPDBline(std::string line, Atom& atom);
};
#endif
