#include <iostream>
#include <cstdlib>
#include <iomanip>
#include <fstream>
#include "PSF.hpp"
#include "PDB.hpp"
#include "Index.hpp"
#include "ComputeRMSD.hpp"
#include "NAMDBin.hpp"
#include "ReadDCD.hpp"

using namespace std;

int main(int argc, char** argv)
{
	if (argc < 6)
	{
		cout <<
		"\nN_TRJ_RMSD\n"
		"\nTrajectory alignment\n"
		"usage: ./a.out psf[natom] dcd ref[pdb/coor] ind [savecrd? yes/no]\n\n";
		return 0;
	}

	cout << "REMARK ";
	for (int i = 0; i < argc; i++)
		cout << argv[i] << " ";
	cout << '\n';

	int natom = 0;
	vector<Atom> atomVector;

	string spsf(argv[1]);
	string spsfex = spsf.substr(spsf.find_last_of(".") + 1);
	if (spsfex == "psf")
	{
		PSF PSFFile(argv[1]);
		atomVector = PSFFile.atomVector;
	}
	else
	{
		natom = atoi(argv[1]);
		atomVector.resize(natom);
	}

	ReadDCD DCD(&atomVector);
	if (!DCD.open_dcd_read(argv[2]))
		return 1;

	vector<Atom> refAtomVector;
	refAtomVector.resize(atomVector.size());

	string sref(argv[3]);
	string srefex = sref.substr(sref.find_last_of(".") + 1);

	if (srefex == "pdb")
	{
		PDB PDBFile(argv[3]);
		if (!PDBFile.LoadCoords(refAtomVector))
			return 1;
	}
	else if (srefex == "coor")
	{
		NAMDBin BIN(argv[3], "coor");
		if (!BIN.read_fi(refAtomVector))
			return 1;
	}
	else
	{
		cerr << "\nerror: unknown file type for ref\n\n";
		return 0;
	}

	cout
	<< "REMARK Reference structure : " << argv[3] << '\n'
	<< "REMARK Target trajectory   : " << argv[2] << '\n';
	
	bool please_record = false;

	string stmp(argv[5]);
	if (stmp == "yes")
	{
		cout << "REMARK Aligned trajectory will be saved in outdcd.\n";
		DCD.open_dcd_write();
		please_record = true;
	}
	else if (stmp == "no")
	{
		cout << "REMARK Aligned trajectory will be discarded.\n";
	}
	else
	{
		err("invalid 5th argument!");
		return 0;
	}

	vector<int> alignIndex;
	
	Index IDX(argv[4]);
	if (!IDX.load_indexes(&alignIndex))
		return 0;

	if (!natom)
	{
		cout << "REMARK Alignment will be done for following " << alignIndex.size() << " atoms: ";
		IDX.show_atoms(&alignIndex, &atomVector);
	}

	ComputeRMSD SOLVER(&alignIndex, &refAtomVector);

	/* remove reference centroid */
	SOLVER.remove_ref_com();
	//PSFFile.writePDB("tmp_ref_vmd.pdb", refAtomVector, "centroid removed reference structure.");
	//cout << "REMARK 重心を除いた参照構造をtmp_ref_vmd.pdbに保存しました。\n";
	cout << scientific << setprecision(8);
	/* alinment start */
	while (DCD.read_1step())
	{
		/* set snapshot coordinates */
		SOLVER.set_tgtAtomVector(&atomVector);
		
		/* now calculate rmsd */
		SOLVER.compute_rmsd();
		
		/* output */
		cout 
		<< setw(10) << DCD._nsteps()
		<< setw(16) << SOLVER._bef_rmsd() 
		<< setw(16) << SOLVER._aft_rmsd() << '\n';
		
		/* save aligned snapshot if you want */
		if (please_record)
			DCD.write_1step();
	}
}
