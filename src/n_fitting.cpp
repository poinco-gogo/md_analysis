#include <iostream>
#include <iomanip>
#include <fstream>
#include "PSF.hpp"
#include "PDB.hpp"
#include "NAMDBin.hpp"
#include "Index.hpp"
#include "ComputeRMSD.hpp"

using namespace std;

int main(int argc, char** argv)
{
	if (argc != 5)
	{
		cout << "\nN_FITTING\n"
		"\nRMSD calculation\n"
		"usage: ./a.out psf tgt[pdb/coor] ref[pdb/coor] ind\n\n";
		return 0;
	}
	cout << "REMARK ";
	for (int i = 0; i < argc; i++)
		cout << argv[i] << ' ';
	cout << '\n';
	
	PSF PSFFile(argv[1]);

	vector<Atom> tgtAtomVector;
	vector<Atom> refAtomVector;

	tgtAtomVector.resize(PSFFile.atomVector.size());
	refAtomVector.resize(PSFFile.atomVector.size());

	string stgt(argv[2]);
	string sref(argv[3]);
	string stgtex = stgt.substr(stgt.find_last_of(".") + 1);
	string srefex = sref.substr(sref.find_last_of(".") + 1);

	if (stgtex == "pdb")
	{
		PDB TGTFile(argv[2]);

		if (!TGTFile.LoadCoords(tgtAtomVector))
			return 0;
	}
	else if (stgtex == "coor")
	{
		NAMDBin BIN(argv[2], "coor");

		BIN.read_fi(tgtAtomVector);
	}
	else
	{
		err("unknown file type");
		return 1;
	}

	if (srefex == "pdb")
	{
		PDB REFFile(argv[3]);

		if (!REFFile.LoadCoords(refAtomVector))
			return 0;
	}
	else if (srefex == "coor")
	{
		NAMDBin BIN(argv[3], "coor");

		BIN.read_fi(refAtomVector);
	}
	else
	{
		err("unknown file type");
		return 1;
	}
	
	cout << "REMARK N_FITTING\n"
	<< "REMARK Reference structure : " << argv[3] << '\n'
	<< "REMARK Target    structure : " << argv[2] << '\n';

	vector<int> alignIndex;
	Index IDX(argv[4]);
	if (!IDX.load_indexes(&alignIndex))
		return 0;
	cout << "REMARK Alignment will be done for following " << alignIndex.size() << " atoms: ";
	IDX.show_atoms(&alignIndex, &PSFFile.atomVector);
	
	
	/* prepare rmsd solver */
	ComputeRMSD SOLVER(&tgtAtomVector, &alignIndex, &refAtomVector);
	/* remove centroid and save the structure */
	SOLVER.remove_ref_com();
	PSFFile.writePDB("tmp_ref.pdb", refAtomVector, "centroid removed reference structure");
	cout << "REMARK 重心を除いた参照構造をtmp_ref.pdbに保存しました。\n";
	
	/* now calculate rmsd */
	SOLVER.compute_rmsd();
	
	/* save aligned target structure */
	PSFFile.writePDB("tmp_tgt.pdb", tgtAtomVector, "aligned target structure");
	cout << "REMARK  align済みの目標構造（重心=0）をtmp_tgt.pdbに保存しました。\n";
	SOLVER.add_ref_com(tgtAtomVector);
	PSFFile.writePDB("tmp_tgt2.pdb", tgtAtomVector, "aligned target structure");
	cout << "REMARK  align済みの目標構造（重心=ref）をtmp_tgt2.pdbに保存しました。\n";
	
	cout << scientific << setprecision(8)
	<< setw(16) << SOLVER._bef_rmsd() 
	<< setw(16) << SOLVER._aft_rmsd() << '\n';
}
