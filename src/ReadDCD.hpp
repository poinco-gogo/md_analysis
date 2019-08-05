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

	std::string ifilename, ofilename;

	double box_size_x, box_size_y, box_size_z;

	std::vector<Atom>* ptr_atomVector;

	dcdhandle* dcd_in;
	dcdhandle* dcd_out;

	molfile_timestep_t ts_in, ts_out;

	public:

	ReadDCD(std::vector<Atom>* ptr_atomVector)
	{
		this->ptr_atomVector = ptr_atomVector;

		this->nsteps         = 0;

		this->box_size_x     = 0;
		this->box_size_y     = 0;
		this->box_size_z     = 0;

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
		if (nsteps)
		{
			dcd_rewind(dcd_in);

			nsteps = 0;
		}
	}

	bool open_dcd_read(std::string filename)
	{
		this->ifilename = filename;

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
				box_size_x = ts_in.A;
				box_size_y = ts_in.B;
				box_size_z = ts_in.C;
			}

			++nsteps;

			return true;
		}
		else
			return false;
	}

	std::string _ifilename() { return ifilename; }

	void open_dcd_write(std::string filename)
	{
		this->ofilename = filename;

		char ctmp;
		int with_unitcell = 0;

		dcd_out = ::open_dcd_write(
				ofilename.c_str(),
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

			ts_out.A = box_size_x;
			ts_out.B = box_size_y;
			ts_out.C = box_size_z;
		}

		::write_timestep(dcd_out, &ts_out);
	}

	std::string _ofilename() { return ofilename; }

	void set_box(const double x, const double y, const double z)
	{
		box_size_x = x;
		box_size_y = y;
		box_size_z = z;
	}

	double _boxx() { return box_size_x; }
	double _boxy() { return box_size_y; }
	double _boxz() { return box_size_z; }
};
#endif
