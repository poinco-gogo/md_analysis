#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <ctime>
#include "MatrixUtil.hpp"
using namespace std;

MatrixUtil::MatrixUtil(double*** matrix, int dimension)
{
	this -> matrix = matrix;
	this -> dimension = dimension;
	memory_allocated = false;
	MemoryAllocation();
}

MatrixUtil::MatrixUtil(string filename, double*** matrix)
{
	this -> matrix = matrix;
	this -> filename = filename;
	memory_allocated = false;
}

MatrixUtil::~MatrixUtil()
{
	MatrixDispose();
}

void MatrixUtil::LoadMatrix()
{
	ifstream fi(filename.c_str());

	if (!fi)
	{
		cerr << "error: file " << filename << "doesn't exist.\n";
		return;
	}

	string s;
	while (getline(fi, s))
	{
		if (s.find("DIMENSION", 0) != string::npos)
		{
			istringstream is(s);
			is >> dimension;

			MemoryAllocation();
			
			LoadFileData(&fi);
		}
	}
	cout << "REMARK LOADED MATRIX DIMENSION = " << dimension << endl;
}

void MatrixUtil::MemoryAllocation()
{
	*matrix = new double*[dimension];
	
	for (int i = 0; i < dimension; i++)
		(*matrix)[i] = new double[dimension];
	
	memory_allocated = true;
}

void MatrixUtil::MatrixReset()
{
	if (memory_allocated)
		for (int i = 0; i < dimension; i++)
			for (int j = 0; j < dimension; j++)
				(*matrix)[i][j] = 0.;
	else
		cout << "error: in MatrixUtil::MatrixReset() memory not allocated!\n";
}

void MatrixUtil::LoadFileData(ifstream* fi)
{
	if (!memory_allocated)
	{
		cerr << "error: in MatrixUtil::LoadFileData() memory not allocated!\n";
		return;
	}

	string s;
	int i = 0;
	int j = 0;
	double a;
	
	while (getline(*fi, s))
	{
		istringstream is(s);
		while(is >> a)
		{
			(*matrix)[i][j] = a;
			(*matrix)[j][i] = (*matrix)[i][j];
			j == dimension - 1 ? j = ++i : ++j;
		}
	}
}

void MatrixUtil::writeMatrix(string filename, string header)
{
	time_t now = time(NULL);
	struct tm* local = localtime(&now);

	ofstream fo(filename.c_str());
	fo << "REMARK " << asctime(local);
	fo << "REMARK " << header << '\n';
	fo << '\n';

	fo << setw(8) << dimension << " !DIMENSION\n";

	// write to filename
	int icnt = 0;
	for (int i = 0; i < dimension; i++)
	{
		for (int j = i; j < dimension; j++)
		{
			fo << scientific << setprecision(8)
			<< setw(16) << (*matrix)[i][j];
			if (icnt++ % 3 == 2)
				fo << '\n';
		}
	}
}

void MatrixUtil::MatrixSetUnit()
{
	if (!memory_allocated)
	{
		cerr << "error: in MatrixUtil::MatrixSetUnit() not allocated!\n";
		return;
	}

	MatrixReset();
	
	for (int i = 0; i < dimension; i++)
		(*matrix)[i][i] = 1.;
}

void MatrixUtil::MatrixDivUP(double d)
{
	if (!memory_allocated)
	{
		cerr << "error: in MatrixUtil::MatrixDivUP() not allocated!\n";
		return;
	}

	for (int i = 0; i < dimension; i++)
		for (int j = i; j < dimension; j++)
			(*matrix)[i][j] /= d;
}


void MatrixUtil::MatrixDispose()
{
	cout << "REMARK MatrixUtil::MatrixDispose() CALLED.\n";
	if (memory_allocated)
	{
		for (int i = 0; i < dimension; i++)
			delete[] (*matrix)[i];
		delete[] *matrix;
		memory_allocated = false;
		cout << "REMARK MEMORY DISALLOCATED.\n";
	}
	else
		cout << "REMARK ...BUT MEMORY NOT ALLOCATED.\n";
}
