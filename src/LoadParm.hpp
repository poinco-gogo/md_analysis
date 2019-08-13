#ifndef ___CLASS_LOADPARM
#define ___CLASS_LOADPARM

#include <vector>
#include "Atom.hpp"
#include "Bond.hpp"
#include "Angle.hpp"
#include "Dihedral.hpp"
#include "Improper.hpp"
#include "Cmap.hpp"

class LoadParm
{
	private:
	std::string filename;
	
	public:
	
	std::vector<Bond>        bondParmVector;
	std::vector<Angle>       angleParmVector;
	std::vector<Dihedral>    dihedralParmVector;
	std::vector<Improper>    improperParmVector;
	std::vector<Cmap>        cmapParmVector;
	std::vector<Atom>        LJParmVector;

	// constructor
	LoadParm(std::string filename);

	bool merge(std::string filename);

	inline std::string _filename() { return filename; }
	
	private:
	void open_fi(std::ifstream& fi);
	void duplication_check();
	void get_bond_parameters(std::ifstream& fi);
	void get_angle_parameters(std::ifstream& fi);
	void get_dihedral_parameters(std::ifstream& fi);
	bool get_improper_parameters(std::ifstream& fi);
	void get_cmap_parameters(std::ifstream& fi);
	void get_LJ_parameters(std::ifstream& fi);
};
#endif
