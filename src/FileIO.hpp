#ifndef ___CLASS_FILEIO
#define ___CLASS_FILEIO

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include "common.hpp"

class FileIO
{
	private:
	std::string filename;

	public:
	FileIO(std::string filename)
	{
		this->filename = filename;
	}

	inline bool load_data(int ncol, std::vector<double>* ptr_dataVector)
	{
		std::ifstream fi(filename.c_str());

		if (!fi)
		{
			std::cerr
			<< "\nerror: Could not open " << filename << "\n\n";
			return false;
		}

		std::string s;
		while (getline(fi, s))
		{
			if (s.empty() || is_comment(s))
				continue;

			std::istringstream is(s);
			double val;
			for (int i = 0; i < ncol; i++)
				is >> val;

			ptr_dataVector->push_back( val );
		}

		return true;
	}
};
#endif

