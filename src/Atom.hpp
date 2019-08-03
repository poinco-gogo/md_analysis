#ifndef ___CLASS_ATOM
#define ___CLASS_ATOM

#include <Eigen/Core>
#include <string>
#include <vector>

class Atom
{
	public:
	
	std::vector<int> neighborVector;      // neighbor atom indices
	std::vector<double> overlapVector;    // overlap area

	int PDBIndex;
	int PSFIndex;

	std::string PSFSegmentName;
	std::string PDBAtomName;
	std::string PSFAtomName;
	std::string PDBResName;
	std::string PSFResName;

	int PDBResID;
	int PSFResID;
	
	Eigen::Vector3d position;
	Eigen::Vector3d velocity;
	Eigen::Vector3d    force;
	Eigen::Vector3d constfrc;
	Eigen::Vector3d     fold;
	Eigen::Vector3d     rnew;
	Eigen::Vector3d     rold;
	Eigen::Vector3d     vnew;
	Eigen::Vector3d     vold;
	Eigen::Vector3d  r_local;

	double A; // Langevin A parameter -i.e., 1 - gamma * dt * 0.5
	double B; // Langevin B parameter -i.e., 1 + gamma * dt * 0.5
	double invB; // 1. / B
	double R; // Langevin noise parameter
	double gamma; // Langevin Damping parameter
	Eigen::Vector3d random_f;

	double PDBo; // PDB occupancy
	double PDBb; // PDB beta-coupling

	double charge;
	double mass;
	double invmass;

	double epsilon;
	double Rmin_div2;
	double eps1_4;
	double Rmin1_4;

	bool is_rigid;

	std::vector<int> exclusionVector;
	std::vector<int> scaled1_4Vector;
	std::vector<int> ex_pair_list;

//	default constructor
	Atom();
//      constructor
	Atom(std::string s);
	// some utilities
	bool checkExclusionPair(const Atom& atom);
	bool checkScaled1_4Pair(const Atom& atom);
	bool is_mainchain();

	private:
	void reset_parameters();
	void parse_PDBfield(std::string PDBField);
};
#endif
