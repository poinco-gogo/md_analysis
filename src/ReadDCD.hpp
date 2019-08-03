#ifndef ___CLASS_READDCD
#define ___CLASS_READDCD

#include <vector>
#include <string>
#include "Atom.hpp"
#include "dcdplugin.h"

class ReadDCD
{
	private:

	int natom, nsets, nsteps;

	std::vector<Atom>* ptr_atomVector;

	dcdhandle* dcd;

	molfile_timestep_t ts;

	public:

	ReadDCD(std::vector<Atom>* ptr_atomVector)
	{
		this->ptr_atomVector = ptr_atomVector;

		this->nsteps         = 0;
	}
	~ReadDCD()
	{
		std::cout << "REMARK ReadDCD destructor called\n";
		close_file_read(dcd);
		free(ts.coords);
	}

	int _nsteps() { return nsteps; }

	bool open_dcd_read(std::string filename)
	{
		char ctmp;

		dcd = ::open_dcd_read(filename.c_str(), &ctmp, &natom, &nsets);

		if (natom != ptr_atomVector->size())
		{
			err("inconsistent atom number between psf and dcd!");
			return false;
		}

		ts.coords = (float *)malloc(dcd->natoms * 3 * sizeof(float));
	}

	bool read_1step()
	{
		if (!read_next_timestep(dcd, natom, &ts))
		{
			int icnt = 0;
			for (auto& at: *ptr_atomVector)
			{
				at.position.x() = dcd->x[icnt];
				at.position.y() = dcd->y[icnt];
				at.position.z() = dcd->z[icnt];
				++icnt;
			}

			++nsteps;

			return true;
		}
		else
			return false;
	}
};
#endif
