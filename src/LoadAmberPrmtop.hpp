#ifndef ___CLASS_LOADAMBERPRMTOP
#define ___CLASS_LOADAMBERPRMTOP

#include <string>
#include <vector>
#include <fstream>
#include "Atom.hpp"

class LoadAmberPrmtop
{
	public:
	int natom, ntypes, nbonh, mbona, ntheth, mtheta, nphih, mphia, nhparm,
	    nparm, nnb, nres, nbona, ntheta, nphia, numbnd, numang, nptra, 
	    natyp, nphb, ifpert, nbper, ngper, ndper, mbper, mgper, mdper,
	    ifbox, nmxrs, ifcap, numextra, ncopy;

	std::vector<std::string> igraph;
	std::vector<double> charge;
	std::vector<double> amass;
	std::vector<int> iac;
	std::vector<int> numex;
	std::vector<int> ico;
	std::vector<std::string> lbres;
	std::vector<int> ipres;
	std::vector<double> rk, req, tk, teq, pk, pn, phase, one_scee, one_scnb;
	std::vector<double> solty, cn1, cn2;
	std::vector<int> ibh, jbh, icbh, ib, jb, icb;
	std::vector<int> ith, jth, kth, icth, it, jt, kt, ict;
	std::vector<int> iph, jph, kph, lph, icph, ip, jp, kp, lp, icp;
	std::vector<int> inb;
	std::vector<double> asol, bsol, hbcut;
	std::vector<std::string> isymbl, itree;
	std::vector<int> join, irotat;
	std::vector<double> rborn, fs;
	int iptres, nspm, nspsol;
	std::vector<int> nsp;
	double oldbeta;
	std::vector<double> box;


	LoadAmberPrmtop(std::string filename);
	void load_prmtop(std::ifstream& fi);
	int get_num_nwt();
	void make_atomVector(std::vector<Atom>& atomVector);
};
#endif
