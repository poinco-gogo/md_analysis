#ifndef ___CLASS_CMAP
#define ___CLASS_CMAP
#include <string>
#include <vector>
#include "Dihedral.hpp"
class Cmap
{
	private:

	static const int res = 24;	

	public:
	Cmap(std::string at1, std::string at2, std::string at3, 
		std::string at4, std::string at5,
		std::string at6, std::string at7, std::string at8)
	{
		attype.push_back(at1); attype.push_back(at2);
		attype.push_back(at3); attype.push_back(at4);
		attype.push_back(at5); attype.push_back(at6);
		attype.push_back(at7); attype.push_back(at8);

		resize_matrix();
	}
	Cmap(Dihedral dihed_phi, Dihedral dihed_psi)
	{
		this -> dihed_phi = dihed_phi;
		this -> dihed_psi = dihed_psi;

		resize_matrix();
	}
	
	Dihedral dihed_phi;
	Dihedral dihed_psi;

	double phi, psi;

	Eigen::Vector3d dPhid1, dPhid2, dPhid3, dPhid4;
	Eigen::Vector3d dPsid1, dPsid2, dPsid3, dPsid4;

	inline void calc_angle_and_derivatives()
	{
		dihed_phi.calc_derivatives_for_cmap();
		dihed_psi.calc_derivatives_for_cmap();

		phi = dihed_phi.chi;
		psi = dihed_psi.chi;

		dPhid1 = dihed_phi.force1;
		dPhid2 = dihed_phi.force2;
		dPhid3 = dihed_phi.force3;
		dPhid4 = dihed_phi.force4;
		
		dPsid1 = dihed_psi.force1;
		dPsid2 = dihed_psi.force2;
		dPsid3 = dihed_psi.force3;
		dPsid4 = dihed_psi.force4;
	}
	
	void show_profile()
	{
		dihed_phi.show_profile();
		dihed_psi.show_profile();
		std::cout << "  cmap_type_index " << cmap_type_index << '\n';
	}
	void show_residue_profile()
	{
		std::cout 
			<< dihed_phi.ptr_atom3->PSFResName
			<< dihed_phi.ptr_atom3->PSFResID 
			<< '\n';
	}

	std::vector<std::string> attype;
	Eigen::MatrixXd cmap_grid_data, cmap_dPhi, cmap_dPsi, cmap_dPhi_dPsi;
	int cmap_type_index;
	Eigen::Vector3d force_phi1, force_phi2, force_phi3, force_phi4;
	Eigen::Vector3d force_psi1, force_psi2, force_psi3, force_psi4;

	private:

	inline void resize_matrix()
	{
		cmap_grid_data.resize( res, res );
		cmap_dPhi.resize( res, res );
		cmap_dPsi.resize( res, res );
		cmap_dPhi_dPsi.resize( res, res );
	}
};
#endif
