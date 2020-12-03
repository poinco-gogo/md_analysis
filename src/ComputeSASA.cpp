#include <iostream>
#include <cstdlib>
#include <cmath>
#include <ctime>
#include <random>
#include "common.hpp"
#include "ComputeSASA.hpp"

using namespace std;

ComputeSASA::ComputeSASA(
		unsigned int npoints,
		unsigned int mode,
		vector<Atom>* ptr_atomVector
		)
{
	this -> npoints        = npoints;
	this -> mode           = mode;
	this -> ptr_atomVector = ptr_atomVector;

	assign_radius();
}

void ComputeSASA::assign_radius()
{
	// assign Bondi radius
	// J. Phys. Chem. 1964, 68(3), 441-451.
	for (auto& atom: *ptr_atomVector)
	{
		char ctmp = atom.PDBAtomName[0];

		if      (ctmp == 'C') atom.rvdw = 1.70;
		else if (ctmp == 'O') atom.rvdw = 1.52;
		else if (ctmp == 'N') atom.rvdw = 1.55;
		else if (ctmp == 'H') atom.rvdw = 1.20;
		else
		{
			die("Detected unsupported atom type.");
		}
		if (atom.PDBAtomName == "CL") atom.rvdw = 1.75;
	}
}

double ComputeSASA::calc_sasa(unsigned int resid)
{
	// see implementation by Bosco Ho (https://github.com/boscoh/pdbremix)

	vector<int> tgt_list;
	for (auto& atom: *ptr_atomVector)
	{
		if (atom.PSFResID == resid)
		{
			tgt_list.push_back(atom.PSFIndex);
		}
	}

	generate_points();

	double rprobe = 1.4;
	double sasa_pre = 4. * N_PI / npoints;
	double sasa = 0;
	for (auto& tgt_id: tgt_list)
	{
		Atom& iatom = ptr_atomVector->at(tgt_id - 1);

		double ri = iatom.rvdw + rprobe;
		
		// find neighbors 
		vector<int> neighbors;
		for (auto& jatom: *ptr_atomVector)
		{
			if (iatom.PSFIndex == jatom.PSFIndex)
				continue;

			double rj = jatom.rvdw + rprobe;

			double d = (iatom.position - jatom.position).norm();

			if ( d < (ri + rj) )
				neighbors.push_back(jatom.PSFIndex);
		}

		int n_accessible_point = 0;
		for (auto& point: points)
		{
			Eigen::Vector3d test_point
			= point * (iatom.rvdw + rprobe) + iatom.position;

			bool is_accessible = true;
			for (auto& j: neighbors)
			{
				Atom& jatom = ptr_atomVector->at(j - 1);

				double d = (jatom.position - test_point).norm();

				double r = jatom.rvdw + rprobe;

				if (d < r)
				{
					is_accessible = false;
					break;
				}
			}
			if (is_accessible)
				++n_accessible_point;
		}

		sasa += ri * ri * n_accessible_point;
	}

	sasa *= sasa_pre;
	
	return sasa;
}

void ComputeSASA::generate_points()
{
	switch (mode) {
		case 1: generate_points_gaussian();  break;
		case 2: generate_points_fibonacci(); break;
		default: die("Unknown mode is specified."); break;
	}
}

void ComputeSASA::generate_points_gaussian()
{
	// uniformly distributed points on the surface of the 3-d unit sphere
	// https://stats.stackexchange.com/questions/7977/how-to-generate-uniformly-distributed-points-on-the-surface-of-the-3-d-unit-sphe

	//unsigned int iseed = 1234;
	mt19937 engine(static_cast<unsigned int>(time(NULL)));
	normal_distribution<double> dist(0., 1.);

        for (int i = 0; i < npoints; i++)
        {
		Eigen::Vector3d vtmp(dist(engine), dist(engine), dist(engine));
		double len = vtmp.norm();
		points.push_back( vtmp / len );
        }
}

void ComputeSASA::generate_points_fibonacci()
{
	// https://stackoverflow.com/questions/9600801/evenly-distributing-n-points-on-a-sphere/26127012

	double inc = N_PI * (3. - sqrt(5));
        double off = 2. / npoints;
        for (int i = 0; i < npoints; i++)
        {
		double y = i * off - 1. + (off / 2.);
                double r = sqrt(1. - y * y);
                double phi = i * inc;
		Eigen::Vector3d vtmp( cos(phi) * r, y, sin(phi) * r );
		points.push_back(vtmp);
        }
}
