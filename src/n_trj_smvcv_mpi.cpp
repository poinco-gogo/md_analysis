#include <boost/mpi.hpp>
#include <boost/mpi/communicator.hpp>
#include <boost/mpi/environment.hpp>
#include <iostream>
#include <cstdlib>
#include <iomanip>
#include <sstream>
#include "ComputeHanning.hpp"
#include "ReadDCD.hpp"
#include "MatrixUtil.hpp"
using namespace std;

int main (int argc, char** argv)
{
	boost::mpi::environment env(argc, argv);
        boost::mpi::communicator world;	

	if (!world.rank() && argc < 7)
	{
		cout << "\nusage: ./a.out natom trj nsmpl_per_win nwin ifirst ilast\n\n";
		env.abort(1);
		return 1;
	}

	if (!world.rank())
	{
		cout << "REMARK ";
		for (int i = 0; i < argc; i++)
			cout << argv[i] << ' ';
		cout << '\n';
	}
	
	int natom         = atoi(argv[1]);
	int nsmpl_per_win = atoi(argv[3]);
	int nwin          = atoi(argv[4]);
	int first         = atoi(argv[5]);
	int last          = atoi(argv[6]);

	int to_do = last - first + 1;

	if (!world.rank() && to_do < world.size())
	{
		cout << "too small size of to_do\n";
		env.abort(1);
		return 1;
	}

	int iifirst = first + (to_do / world.size()) * world.rank();
	int ilast  = first + (to_do / world.size()) * (world.rank() + 1) - 1;

        if (world.rank() == world.size() - 1 && to_do % world.size())
                ilast += to_do % world.size();
	
	double* xcod = new double [natom * nsmpl_per_win * nwin];
	double* ycod = new double [natom * nsmpl_per_win * nwin];
	double* zcod = new double [natom * nsmpl_per_win * nwin];

	vector<Atom> atomVector(natom);
	ReadDCD DCD(&atomVector);
	if (!DCD.open_dcd_read(argv[2]))
		return 1;

	double** vcvmat;
	MatrixUtil JOBVCV(&vcvmat, natom * 3);
	
	for (int ifirst = iifirst; ifirst <= ilast; ifirst++){
	
	DCD.read_rewind();

	JOBVCV.MatrixReset();

	vector<double> hann_wt(nsmpl_per_win, 0.);
	ComputeHanning JOB(&hann_wt);

	for (int i = 0; i < natom * nsmpl_per_win * nwin; i++)
	{
		xcod[i] = 0;
		ycod[i] = 0;
		zcod[i] = 0;
	}

	int icnt = 0;
	while (DCD.read_1step())
	{
		if (DCD._nsteps() < ifirst)
			continue;
		
		for (int i = 0; i < natom; i++)
		{
			int bs = icnt * natom + i;
			xcod[bs] = atomVector[i].position.x();
			ycod[bs] = atomVector[i].position.y();
			zcod[bs] = atomVector[i].position.z();
		}

		++icnt;
		
		
		if (DCD._nsteps() == ifirst + nsmpl_per_win * nwin - 1)
			break;
	}

	cout << "REMARK " << icnt << " coordinates loaded.\n";


	cout << setprecision(8) << fixed;

	for (int i = 0; i < nwin; i++)
	{
		double* axcod = new double [natom];
		double* aycod = new double [natom];
		double* azcod = new double [natom];
		for (int j = 0; j < natom; j++)
		{
			axcod[j] = 0;
			aycod[j] = 0;
			azcod[j] = 0;
		}

		for (int j = 0; j < nsmpl_per_win; j++)
		{
			for (int k = 0; k < natom; k++)
			{
				int bs = i * nsmpl_per_win * natom 
						+ j * natom + k;
				axcod[k] += hann_wt[j] * xcod[bs];
				aycod[k] += hann_wt[j] * ycod[bs];
				azcod[k] += hann_wt[j] * zcod[bs];
			}
		}
		for (int j = 0; j < natom; j++)
		{
			axcod[j] *= JOB._indenom();
			aycod[j] *= JOB._indenom();
			azcod[j] *= JOB._indenom();
		}

		double* buff = new double[natom * 3];
		for (int j = 0; j < nsmpl_per_win; j++)
		{
			for (int k = 0; k < natom; k++)
			{
				int bs = i * nsmpl_per_win * natom 
						+ j * natom + k;
				buff[k * 3    ] = xcod[bs] - axcod[k];
				buff[k * 3 + 1] = ycod[bs] - aycod[k];
				buff[k * 3 + 2] = zcod[bs] - azcod[k];
			}

			for (int k = 0; k < natom * 3; k++)
			{
				double dtmp = buff[k];
				for (int l = k; l < natom * 3; l++)
				{
					vcvmat[k][l] += dtmp * buff[l];
				}
			}
		}
		delete[] buff;
		delete[] axcod, aycod, azcod;
	}

	double dstep = static_cast<double>(nwin * nsmpl_per_win);
	JOBVCV.MatrixDivUP(dstep);
	ostringstream os;
	os << "tmp" << ifirst << ".vcv";
	JOBVCV.writeMatrix(os.str(), "tmp data");


	}

	JOBVCV.MatrixDispose();

	delete[] xcod, ycod, zcod;
}
