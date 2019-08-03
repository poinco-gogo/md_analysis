#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <Eigen/Geometry>
#include <Eigen/Dense>
#include "ComputeRMSD.hpp"

using namespace std;

ComputeRMSD::ComputeRMSD(vector<Atom>* ptr_tgtAtomVector)
{
	this->ptr_tgtAtomVector  = ptr_tgtAtomVector;

	resetParm();
}

ComputeRMSD::ComputeRMSD(vector<int>* ptr_indexVector, vector<Atom>* ptr_refAtomVector)
{
	this-> ptr_indexVector   = ptr_indexVector;
	this-> ptr_refAtomVector = ptr_refAtomVector;

	resetParm();
}

ComputeRMSD::ComputeRMSD(vector<Atom>* ptr_tgtAtomVector, vector<int>* ptr_indexVector, vector<Atom>* ptr_refAtomVector)
{
	this->ptr_tgtAtomVector  = ptr_tgtAtomVector;
	this->ptr_indexVector    = ptr_indexVector;
	this->ptr_refAtomVector  = ptr_refAtomVector;

	resetParm();
}

void ComputeRMSD::resetParm()
{
	bef_rmsd = 0.0;
	aft_rmsd = 0.0;
	R      = M3ZERO;
	refcom = V3ZERO;
	tgtcom = V3ZERO;
}

void ComputeRMSD::set_tgtAtomVector(vector<Atom>* ptr_tgtAtomVector)
{
	this->ptr_tgtAtomVector  = ptr_tgtAtomVector;
}

void ComputeRMSD::set_refAtomVector(vector<Atom>* ptr_refAtomVector)
{
	this->ptr_refAtomVector  = ptr_refAtomVector;
}

void ComputeRMSD::set_indexVector(vector<int>* ptr_indexVector)
{
	this->ptr_indexVector    = ptr_indexVector;
}

void ComputeRMSD::remove_ref_com()
{
	refcom = V3ZERO;

	for (auto& i: *ptr_indexVector)
		refcom += ptr_refAtomVector -> at( i - 1 ).position;
	refcom = refcom / ptr_indexVector->size();

	for (auto& at: *ptr_refAtomVector)
		at.position -= refcom;
}

void ComputeRMSD::remove_tgt_com()
{
	tgtcom = V3ZERO;

	for (auto& i: *ptr_indexVector)
		tgtcom += ptr_tgtAtomVector -> at( i - 1 ).position;
	tgtcom = tgtcom / ptr_indexVector->size();

	for (auto& at: *ptr_tgtAtomVector)
		at.position -= tgtcom;
}

void ComputeRMSD::get_rotation_matrix()
{
	Eigen::Matrix3d C = M3ZERO;

	for (auto& i: *ptr_indexVector)
	{
		Atom& tgt = ptr_tgtAtomVector->at( i - 1 );
		Atom& ref = ptr_refAtomVector->at( i - 1 );

		C(0,0) += tgt.position.x() * ref.position.x();
		C(0,1) += tgt.position.x() * ref.position.y();
		C(0,2) += tgt.position.x() * ref.position.z();
		C(1,0) += tgt.position.y() * ref.position.x();
		C(1,1) += tgt.position.y() * ref.position.y();
		C(1,2) += tgt.position.y() * ref.position.z();
		C(2,0) += tgt.position.z() * ref.position.x();
		C(2,1) += tgt.position.z() * ref.position.y();
		C(2,2) += tgt.position.z() * ref.position.z();
	}

	double det = C(0,0)*(C(1,1)*C(2,2)-C(1,2)*C(2,1))
		   - C(0,1)*(C(1,0)*C(2,2)-C(1,2)*C(2,0))
		   + C(0,2)*(C(1,0)*C(2,1)-C(1,1)*C(2,0));
	double sign = (det > 0 ? 1.0 : -1.0);

	Eigen::JacobiSVD<Eigen::Matrix3d> svd(C, Eigen::ComputeFullU | 
			Eigen::ComputeFullV);

	Eigen::Matrix3d U  = svd.matrixU();
	Eigen::Matrix3d VT = svd.matrixV().transpose();

	R = M3ZERO;

	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			R(i,j) += VT(0,i) * U(j,0) + VT(1,i) * U(j,1)
				+ sign * VT(2,i) * U(j,2);
		}
	}
}

void ComputeRMSD::compute_rmsd()
{
	double vsqrm0 = 0.;

	remove_tgt_com();

	for (auto& i: *ptr_indexVector)
	{
		vsqrm0 +=
		(ptr_tgtAtomVector -> at(i - 1).position - ptr_refAtomVector -> at(i - 1).position).squaredNorm();
	}

	bef_rmsd = sqrt( vsqrm0 / ptr_indexVector->size() );

	get_rotation_matrix();
	
	for (auto& at: *ptr_tgtAtomVector)
		at.position = R * at.position;
		
	vsqrm0 = 0.;
	for (auto& i: *ptr_indexVector)
	{
		vsqrm0 += 
		(ptr_tgtAtomVector -> at(i - 1).position - ptr_refAtomVector -> at(i - 1).position).squaredNorm();
	}

	aft_rmsd = sqrt( vsqrm0 / ptr_indexVector->size() );
}

void ComputeRMSD::add_ref_com(vector<Atom>& atomVector)
{
	for (auto& at: atomVector)
		at.position += refcom;
}
