#ifndef ___CLASS_INDEX
#define ___CLASS_INDEX

#include <string>
#include <vector>
#include "Atom.hpp"

class Index
{
	private:
	std::string filename;

	public:
	Index(std::string filename);

	int load_indexes(std::vector<int>* ptr_indexVector);

	void show_atoms(std::vector<int>* ptr_indexVector, std::vector<Atom>* ptr_atomVector);

};
#endif
