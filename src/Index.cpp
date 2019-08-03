#include <iostream>
#include <fstream>
#include "Index.hpp"
using namespace std;

Index::Index(string filename)
{
	this->filename           = filename;
}

int Index::load_indexes(vector<int>* ptr_indexVector)
{
	ifstream fi(filename.c_str());
	if (!fi)
	{
		cerr << "\nerror: file " << filename << " not exists.\n\n";
		return 0;
	}

	int itmp;
	while (fi >> itmp)
		ptr_indexVector -> push_back(itmp);

	return ptr_indexVector -> size();
}

void Index::show_atoms(vector<int>* ptr_indexVector, vector<Atom>* ptr_atomVector)
{
	for (auto& i: *ptr_indexVector)
	{
		Atom& at = ptr_atomVector->at( i - 1 );

		cout
			<< at.PSFIndex
			<< at.PDBAtomName
			<< ' ';
	}

	cout << '\n';
}
