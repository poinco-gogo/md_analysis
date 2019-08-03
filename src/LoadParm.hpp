#ifndef ___CLASS_LOADPARM
#define ___CLASS_LOADPARM

#include <vector>
#include "Atom.hpp"
#include "Bond.hpp"
#include "Angle.hpp"

class LoadParm
{
	private:
	std::string filename;
	
	public:
	
	std::vector<Bond> bondParmVector;
	std::vector<Angle> angleParmVector;
	std::vector<Atom> LJParmVector;

	// constructor
	LoadParm(std::string filename);
	bool merge(std::string filename);

	//member
	inline std::string _filename() { return filename; }

	
	private:
	// for internal use
	void open_fi(std::ifstream& fi);
	void duplication_check();
	void get_bond_parameters(std::ifstream& fi);
	void get_angle_parameters(std::ifstream& fi);
	void get_cmap_parameters(std::ifstream& fi);
	void get_LJ_parameters(std::ifstream& fi);
};
#endif
