#ifndef ___CLASS_READDCD
#define ___CLASS_READDCD

#include <vector>
#include <string>
#include "Atom.hpp"
#include "uiuc/dcdplugin.h"

class ReadDCD
{
	private:

	int natom, nsets;

	std::vector<Atom>* ptr_atomVector;

	dcdhandle* dcd;

	molfile_timestep_t ts;

	public:

	ReadDCD(std::vector<Atom>* ptr_atomVector)
	{
		this->ptr_atomVector = ptr_atomVector;
	}
	~ReadDCD()
	{
		std::cout << "REMARK ReadDCD destructor called\n";
		close_file_read(dcd);
		free(ts.coords);
	}

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
			for (auto& at: *ptr_atomVector)
			{
				at.position.x() = dcd->x[i];
				at.position.y() = dcd->y[i];
				at.position.z() = dcd->z[i];
			}

			return true;
		}
		else
			return false;
	}
};
#endif
