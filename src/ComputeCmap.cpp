// most codes from AMBER cmap.f90
#include <iostream>
#include <cmath>
#include "Cmap.hpp"
#include "ComputeSpline.hpp"
#include "ComputeCmap.hpp"
using namespace std;
ComputeCmap::ComputeCmap(LoadParm& All22, vector<Cmap>& cmapVector)
{
	this -> ptr_All22        = &All22;
	this -> ptr_cmapVector   = &cmapVector;

	generate_cmap_grid_derivatives();
}

void ComputeCmap::generate_cmap_grid_derivatives()
{
	// 各グリッド点における 1 階微分と 2 階微分を求める

	for (Cmap& cmap_parm: ptr_All22->cmapParmVector)
	{
		//  calculate dE/dPhi
		for (int row = 0; row < this->res; row++)
		{
			double tmpy[two_res];
			for (int k = 0; k < res; k++)
			{
				tmpy[k + half_res] 
					= cmap_parm.cmap_grid_data( row, k );
			}
			for (int k = 0; k < half_res; k++)
			{
				tmpy[k] = tmpy[k + res];
				tmpy[k + half_res + res] = tmpy[k + half_res];
			}

			ComputeSpline JOB(two_res, gridOrigin, stepSize, tmpy);
			JOB.calc_1d_coefficients();

			for (int j = half_res; j < res+half_res; j++)
			{
				cmap_parm.cmap_dPhi( row, j-half_res )
				= JOB.calc_gradient(gridOrigin + j * stepSize);
			}
		}
		
		// calculate dE/dPsi	
		for (int col = 0; col < res; col++)
		{
			double tmpy[two_res];
			for (int k = 0; k < res; k++)
			{
				tmpy[k + half_res]
					= cmap_parm.cmap_grid_data( k, col );
			}
			for (int k = 0; k < half_res; k++)
			{
				tmpy[k] = tmpy[k + res];
				tmpy[k + half_res + res] = tmpy[k + half_res];
			}
			
			ComputeSpline JOB(two_res, gridOrigin, stepSize, tmpy);
			JOB.calc_1d_coefficients();
				
			for (int j = half_res; j < res + half_res; j++)
			{
				cmap_parm.cmap_dPsi( j-half_res, col )
				= JOB.calc_gradient(gridOrigin + j * stepSize);
			}
		}

		// calculate d^2E/dPhidPsi
		for (int col = 0; col < res; col++)
		{
			double tmpy[two_res];
			for (int k = 0; k < res; k++)
			{
				tmpy[k + half_res]
				       = cmap_parm.cmap_dPhi( k, col );
			}
			for (int k = 0; k < half_res; k++)
			{
				tmpy[k] = tmpy[k + res];
				tmpy[k + half_res + res] = tmpy[k + half_res];
			}

			ComputeSpline JOB(two_res, gridOrigin, stepSize, tmpy);
			JOB.calc_1d_coefficients();
				
			for (int j = half_res; j < res + half_res; j++)
			{
				cmap_parm.cmap_dPhi_dPsi( j-half_res, col )
				= JOB.calc_gradient(gridOrigin + j * stepSize);
			}
		}
	}
}

double ComputeCmap::compute_force()
{
	double sum_energy = 0;

	for (Cmap& cmap: *ptr_cmapVector)
	{
		cmap.calc_angle_and_derivatives();

		double phi = cmap.phi * RAD2DEG;
		double psi = cmap.psi * RAD2DEG;
		double dEdPhi = 9999;
		double dEdPsi = 9999;
		
		int id = cmap.cmap_type_index;

		sum_energy += calc_cmap(phi, psi, id, dEdPhi, dEdPsi);

		dEdPhi *= RAD2DEG;
		dEdPsi *= RAD2DEG;

		cmap.force_phi1 = -dEdPhi * cmap.dPhid1;
		cmap.force_phi2 = -dEdPhi * cmap.dPhid2;
		cmap.force_phi3 = -dEdPhi * cmap.dPhid3;
		cmap.force_phi4 = -dEdPhi * cmap.dPhid4;

		cmap.force_psi1 = -dEdPsi * cmap.dPsid1;
		cmap.force_psi2 = -dEdPsi * cmap.dPsid2;
		cmap.force_psi3 = -dEdPsi * cmap.dPsid3;
		cmap.force_psi4 = -dEdPsi * cmap.dPsid4;

		cmap.dihed_phi.ptr_atom1->force += cmap.force_phi1;
		cmap.dihed_phi.ptr_atom2->force += cmap.force_phi2;
		cmap.dihed_phi.ptr_atom3->force += cmap.force_phi3;
		cmap.dihed_phi.ptr_atom4->force += cmap.force_phi4;
		
		cmap.dihed_psi.ptr_atom1->force += cmap.force_psi1;
		cmap.dihed_psi.ptr_atom2->force += cmap.force_psi2;
		cmap.dihed_psi.ptr_atom3->force += cmap.force_psi3;
		cmap.dihed_psi.ptr_atom4->force += cmap.force_psi4;
	}

	return sum_energy;
}

double ComputeCmap::calc_cmap(const double phi, const double psi, const int cmap_type, double& dEdPhi, double& dEdPsi)
{
	// 与えられた点 ( phi, psi ) を中に含む 4 つのグリッド点を求める

	// 四隅のうち、左下の点の X, Y インテックス（0-based）
	int x0 = floor( ( phi - gridOrigin ) / stepSize );
	int y0 = floor( ( psi - gridOrigin ) / stepSize );

	// 四隅のうち、右上の点の X, Y インデックス（0-based）
	int x1 = (x0 + 1) % res; // x0 = 23 のとき x1 = 0 である（PBC 的に）
	int y1 = (y0 + 1) % res;

	Cmap& cmap_parm = ptr_All22->cmapParmVector[cmap_type];

	// 四隅のグリッド点における CMAP エネルギーの値
	double grid_E[4] = {
				cmap_parm.cmap_grid_data( y0, x0 ),
				cmap_parm.cmap_grid_data( y0, x1 ),
				cmap_parm.cmap_grid_data( y1, x1 ),
				cmap_parm.cmap_grid_data( y1, x0 ),
	};

	double grid_dPhi[4] = {
				cmap_parm.cmap_dPhi( y0, x0 ),
				cmap_parm.cmap_dPhi( y0, x1 ),
				cmap_parm.cmap_dPhi( y1, x1 ),
				cmap_parm.cmap_dPhi( y1, x0 ),
	};

	double grid_dPsi[4] = {
				cmap_parm.cmap_dPsi( y0, x0 ),
				cmap_parm.cmap_dPsi( y0, x1 ),
				cmap_parm.cmap_dPsi( y1, x1 ),
				cmap_parm.cmap_dPsi( y1, x0 ),
	};

	double grid_dPhi_dPsi[4] = {
				cmap_parm.cmap_dPhi_dPsi( y0, x0 ),
				cmap_parm.cmap_dPhi_dPsi( y0, x1 ),
				cmap_parm.cmap_dPhi_dPsi( y1, x1 ),
				cmap_parm.cmap_dPhi_dPsi( y1, x0 ),
	};

	// spline coefficients
	double a[4][4] = {9999.};

	ComputeSpline JOB;
	JOB.calc_2d_coefficients(grid_E, grid_dPhi, grid_dPsi, grid_dPhi_dPsi, stepSize, a);

	double E = 0.;
	dEdPhi = 0;
	dEdPsi = 0;

	double t = ( phi - (gridOrigin + stepSize * x0) ) / stepSize;
	double u = ( psi - (gridOrigin + stepSize * y0) ) / stepSize;

	for (int i = 3; i >= 0; i--)
	{
		E = E * t + ( (a[i][3]*u+a[i][2])*u+a[i][1] )*u+a[i][0];
		dEdPhi = dEdPhi * u + ( 3.*a[3][i]*t+2.*a[2][i] )*t + a[1][i];
		dEdPsi = dEdPsi * t + ( 3.*a[i][3]*u+2.*a[i][2] )*u + a[i][1];
	}

	dEdPhi /= stepSize;
	dEdPsi /= stepSize;

	return E;
}
