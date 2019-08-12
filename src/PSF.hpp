#ifndef ___CLASS_PSF
#define ___CLASS_PSF

#include <Eigen/Core>
#include <string>
#include <fstream>
#include <vector>
#include "Atom.hpp"
#include "Bond.hpp"
#include "Angle.hpp"
#include "Dihedral.hpp"
#include "Improper.hpp"
#include "Cmap.hpp"

class PSF 
{
	private:
	int natom, nwater;
	std::string filename;
	std::vector<Atom>* ptr_atomVector;
	
	public:
	std::vector<Atom>      atomVector;
	std::vector<Bond>      bondVector;
	std::vector<Angle>     angleVector;
	std::vector<Dihedral>  dihedralVector;
	std::vector<Improper>  improperVector;
	std::vector<Cmap>      cmapVector;

	std::vector<int> bondArray;
	std::vector<int> angleArray;
	std::vector<int> dihedralArray;
	std::vector<int> improperArray;
	std::vector<int> cmapArray;
	
	// constructor
	PSF(std::string filename);
	PSF(std::string filename, std::vector<Atom>* ptr_atomVector);

	// member function
	public:

	bool set_bond_parm(std::vector<Bond>& bondParmVector);
	bool set_angle_parm(std::vector<Angle>& angleParmVector);
	bool set_dihedral_parm(std::vector<Dihedral>& dihedralParmVector);
	bool set_improper_parm(std::vector<Improper>& improperParmVector);
	void set_cmap_parm(std::vector<Cmap>& cmapParmVector);
	bool set_lj_parm(std::vector<Atom>& LJParmVector);

	void make_exclusion_vector();
	void make_scaled1_4_vector();

	void showAtom(const Atom& atom);
	void showResidue(int resID);

	void writePDB(std::string filename, std::string header);
	void writePDB(std::string filename, std::string header, const std::vector<int>& indexVector);
	void writePDB(std::string filename, std::vector<Atom>& atomVector, std::string header);

	void make_water_shake_bond_array();

	Eigen::Vector3d getcomof_zerobased(const std::vector<int>& ind);

	double _nwater() { return nwater; }

	// for internal use
	private:
	void read_PSF();
	void get_atom_list(const int natom, std::ifstream& fi);
	void get_bond_list(const int nbond, std::ifstream& fi);
	void get_angle_list(const int nangle, std::ifstream& fi);
	void get_dihedral_list(const int ndihedral, std::ifstream& fi);
	void get_improper_list(const int nimproper, std::ifstream& fi);
	void get_cmap_list(const int ncmap, std::ifstream& fi);
};
#endif
