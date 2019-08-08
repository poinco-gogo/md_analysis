#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include "common.hpp"
#include "MatrixUtil.hpp"
#include "Index.hpp"
#include "PSF.hpp"
#include "ForceUtil.hpp"
using namespace std;
int main(int argc, char** argv)
{
	if (argc != 5)
	{
		cout <<
		"\nGET LINEAR RESPONSE DISP VECTORS\n"
		"usage: ./a.out psf vcv.mr4.up mr4.ind mr4force\n\n";
		return 1;
	}
	
	cout << "REMARK ";
	for (int i = 0; i < argc; i++)
		cout << argv[i] << " ";
	cout << "\nREMARK GET LRT DISPLACEMENT VECTORS\n";

	vector<Atom> atomVector;
	PSF PSFFile(argv[1], &atomVector);

	// Load matrix file
	double** MATRIX;
	string filename(argv[2]);
	MatrixUtil JOB(filename, &MATRIX);
	JOB.LoadMatrix();
	
	// Load site force
	vector<Eigen::Vector3d> vfrc;
	ForceUtil ff;
	if (ff.LoadForce(argv[4], vfrc)) return 0;

	cout << "REMARK Number of force vectors = " << vfrc.size() << '\n';

	if (vfrc.size() != JOB._dimension() / 3)
	{
		cerr << "error: vfrc.size != JOB.dimension / 3\n";
		return 0;
	}

	Index idx(argv[3]);
	vector<int> IDVector;
	if (!idx.load_indexes(&IDVector)) return 0;

	cout <<
	"REMARK Number of input indexes = " << IDVector.size() << '\n';

	// consistency check
	if (JOB._dimension() != IDVector.size() * 3)
	{
		cerr << "error: matrix size != index size * 3\n";
		return 0;
	}

	cout << "REMARK OUTPUT LRT DISPLACEMENT VECTOR DX DY DZ vabs(DR)\n";

	ForceUtil FRC; // For output.
	
	for (int i = 0; i < IDVector.size(); i++)
	{
		int ii = IDVector[i] - 1;
		
		Eigen::Vector3d dr(0., 0., 0.);
		
		for (int j = 0; j < IDVector.size(); j++)
		{
			dr.x() += MATRIX[i * 3    ][j * 3    ] * vfrc[j].x()
			        + MATRIX[i * 3    ][j * 3 + 1] * vfrc[j].y()
			        + MATRIX[i * 3    ][j * 3 + 2] * vfrc[j].z();
			dr.y() += MATRIX[i * 3 + 1][j * 3    ] * vfrc[j].x()
			        + MATRIX[i * 3 + 1][j * 3 + 1] * vfrc[j].y()
			        + MATRIX[i * 3 + 1][j * 3 + 2] * vfrc[j].z();
			dr.z() += MATRIX[i * 3 + 2][j * 3    ] * vfrc[j].x()
			        + MATRIX[i * 3 + 2][j * 3 + 1] * vfrc[j].y()
			        + MATRIX[i * 3 + 2][j * 3 + 2] * vfrc[j].z();	
		}

		dr = dr / (BOLTZMAN * 298.0);

		FRC.writeForce("stdout", atomVector[ii], dr);
	}
}
