#ifndef ___CLASS_COMPUTERMSD
#define ___CLASS_COMPUTERMSD

#include <Eigen/Core>
#include <vector>
#include "Atom.hpp"

class ComputeRMSD
{
	private:
	double bef_rmsd;
	double aft_rmsd;
	std::vector<Atom>* ptr_tgtAtomVector;
	std::vector<int>*     ptr_indexVector;
	std::vector<Atom>* ptr_refAtomVector;
	Eigen::Vector3d refcom, tgtcom;
	Eigen::Matrix3d R;

	public:

	ComputeRMSD(std::vector<Atom>* ptr_tgtAtomVector);
	ComputeRMSD(std::vector<int>* ptr_indexVector, std::vector<Atom>* ptr_refAtomVector);
	ComputeRMSD(std::vector<Atom>* ptr_tgtAtomVector, std::vector<int>* ptr_indexVector, std::vector<Atom>* ptr_refAtomVector);
	
	void resetParm();

	void set_tgtAtomVector(std::vector<Atom>* ptr_tgtAtomVector);
	void set_refAtomVector(std::vector<Atom>* ptr_refAtomVector);
	void set_indexVector(std::vector<int>* ptr_indexVector);

	void add_ref_com(std::vector<Atom>& atomVector);
	void remove_ref_com();
	void remove_tgt_com();

	void get_rotation_matrix();

	void compute_rmsd();

	double _bef_rmsd() { return bef_rmsd; }
	double _aft_rmsd() { return aft_rmsd; }
};
#endif
