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

	bool open_fi, open_fo;

	std::vector<Atom>* ptr_atomVector;

	dcdhandle* dcd_in;
	dcdhandle* dcd_out;

	molfile_timestep_t ts_in, ts_out;

	public:

	ReadDCD(std::vector<Atom>* ptr_atomVector)
	{
		this->ptr_atomVector = ptr_atomVector;

		this->nsteps         = 0;

		this->open_fi        = false;
		this->open_fo        = false;
	}
	~ReadDCD()
	{
		std::cout << "REMARK ReadDCD destructor called\n";

		if (open_fi)
		{
			close_file_read(dcd_in);
			free(ts_in.coords);
		}

		if (open_fo)
		{
			close_file_read(dcd_out);
			free(ts_out.coords);
		}
	}

	int _nsteps() { return nsteps; }

	void read_rewind()
	{
		dcd_rewind(dcd_in);

		nsteps = 0;
	}

	bool open_dcd_read(std::string filename)
	{
		char ctmp;

		dcd_in
		= ::open_dcd_read(filename.c_str(), &ctmp, &natom, &nsets);

		if (natom != ptr_atomVector->size())
		{
			err("inconsistent atom number between psf and dcd!");
			return false;
		}

		ts_in.coords
			= (float *)malloc(dcd_in->natoms * 3 * sizeof(float));

		open_fi = true;

		return true;
	}

	bool read_1step()
	{
		if (!read_next_timestep(dcd_in, natom, &ts_in))
		{
			int icnt = 0;
			for (auto& at: *ptr_atomVector)
			{
				at.position.x() = dcd_in->x[icnt];
				at.position.y() = dcd_in->y[icnt];
				at.position.z() = dcd_in->z[icnt];
				++icnt;
			}

			++nsteps;

			return true;
		}
		else
			return false;
	}

	void open_dcd_write()
	{
		char ctmp;
		int with_unitcell = 0;

		dcd_out = ::open_dcd_write(
				"outdcd",
				&ctmp,
				ptr_atomVector->size(),
				with_unitcell );

		ts_out.coords
			= (float *)malloc(dcd_out->natoms * 3 * sizeof(float));

		open_fo = true;
	}

	void write_1step()
	{
		for (int i = 0; i < ptr_atomVector->size(); i++)
		{
			Atom& at = ptr_atomVector->at(i);

			ts_out.coords[3 * i    ] = at.position.x();
			ts_out.coords[3 * i + 1] = at.position.y();
			ts_out.coords[3 * i + 2] = at.position.z();
		}

		::write_timestep(dcd_out, &ts_out);
	}
};
#endif
