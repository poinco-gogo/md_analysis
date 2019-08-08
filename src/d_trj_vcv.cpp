#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <ctime>
#include <cstdlib>
#include <Eigen/Eigenvalues>
#include "MatrixUtil.hpp"
#include "PSF.hpp"
#include "PDB.hpp"
#include "NAMDBin.hpp"
#include "ReadDCD.hpp"
#include "Index.hpp"
using namespace std;
int main(int argc, char** argv)
{
	if (argc != 8)
	{
		cout << "\nD_TRJ_VCV\n"
		"\nVariance covariance calculation\n"
		"usage : ./a.out psf[natom] dcd ave[pdb/coor] ofilename first last ind\n\n"
		"output file is: ofilename.vcv.up\n\n";
		return 1;
	}

	cout << "REMARK ";
	for (int i = 0; i < argc; i++)
		cout << argv[i] << ' ';
	cout << '\n';

	int natom = 0;
	string spsf(argv[1]);
        string spsfex = spsf.substr(spsf.find_last_of(".") + 1);
        if (spsfex == "psf")
        {
                PSF PSFFile(argv[1]);
                natom = PSFFile.atomVector.size();
        }
        else
        {
                natom = atoi(argv[1]);
        }

	vector<Atom> refAtomVector(natom);
	
	string ofilename(argv[4]);
	
	vector<int> alignIndex;
	Index IDX(argv[7]);
	if (!IDX.load_indexes(&alignIndex))
		return 1;
	cout << "REMARK " << alignIndex.size() << " indexes loaded.\n";

	string sref(argv[3]);
	string srefex = sref.substr(sref.find_last_of(".") + 1);
	if (srefex == "pdb")
	{
		PDB PDBFile(argv[3]);
		if (!PDBFile.LoadCoords(refAtomVector))
			return 0;
	}
	else
	{
		NAMDBin JOB(argv[3], "coor");
		JOB.read_fi(refAtomVector);
	}
	
	vector<Atom> tgtAtomVector(natom);
	ReadDCD TRJFile(&tgtAtomVector);
	if (!TRJFile.open_dcd_read(argv[2]))
		return 0;

	cout << "REMARK DCD file is : " << TRJFile._ifilename() << '\n';
	cout << "REMARK vcv will be written in file : " << argv[4]
	<< ".vcv.up\n";

	// vcv calculation start...

	int dim = alignIndex.size() * 3;

	//  need a lot of memory...
	double* buff = new double[dim];

	double** vcvmat;

	MatrixUtil JOB(&vcvmat, dim);

	// set all element zero
	JOB.MatrixReset();

	int istep = 0;

	int ifirst = atoi(argv[5]);
	int ilast  = atoi(argv[6]);
	cout << "REMARK Calculation begins at " << ifirst << '\n';
	cout << "REMARK Calculation   ends at " << ilast  << '\n';

	while (TRJFile.read_1step())
	{
		if (TRJFile._nsteps() < ifirst) continue;

		Eigen::Vector3d vtmp(0., 0., 0.);
		for (int i = 0; i < dim / 3; i++)
		{
			int j = alignIndex[i] - 1;
			vtmp = tgtAtomVector[j].position 
				- refAtomVector[j].position;
			buff[i * 3    ] = vtmp.x();
			buff[i * 3 + 1] = vtmp.y();
			buff[i * 3 + 2] = vtmp.z();
		}
		for (int i = 0; i < dim; i++)
		{
			double dtmp = buff[i];
			for (int j = i; j < dim; j++)
			{
				vcvmat[i][j] += dtmp * buff[j];
			}
		}

		if (TRJFile._nsteps() == ilast) break;
	}
	istep = TRJFile._nsteps();

	cout << "REMARK " << istep << " coords loaded.\n";

	string header = "variance covariance matrix\n";
	header += "REMARK ave coord: " + sref
		+ "\nREMARK trajectory file: " + argv[2];
	
	double dstep = static_cast<double>(istep);
	
	// take average	
	JOB.MatrixDivUP(dstep);

	// write to file
	JOB.writeMatrix(ofilename + ".vcv.up", header);
	
	Eigen::MatrixXd m(dim, dim);
	for (int i = 0; i < dim; i++)
	{
		for (int j = i; j < dim; j++)
		{
			m(i, j) = vcvmat[i][j];
			m(j, i) = m(i, j);
		}
	}

        //  destruct objects
	delete[] buff;
	JOB.MatrixDispose();

	// EIGEN CALCULATION
	Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(m);

	double eigsum = 0.;
	for (int i = 0; i < dim; i++)
		eigsum += es.eigenvalues()[i];

	ofstream fo;
	fo.open((ofilename + ".vcv.eval").c_str());

	fo << "REMARK ID EIGENVALUE CONTRIBUTION(%) ACCUMULATION\n";
	fo << setw(8) << dim << " !DIMENSION\n\n";

	double contr = 0;
	double accum = 0;
	for (int i = dim - 1; i >= 0; i--)
        {
                fo << setw(8) << i + 1 <<
                scientific << setprecision(8) << setw(16)
                << es.eigenvalues()[i] << setw(16);
                contr =  es.eigenvalues()[i] * 100. / eigsum;
                fo << setw(16) << contr;
                accum += contr;
                fo << setw(16) << accum << endl;
        }

	fo.close();

	fo.open((ofilename + ".vcv.evec").c_str());

	fo << scientific << setprecision(8) << es.eigenvectors();
	fo.close();
}
