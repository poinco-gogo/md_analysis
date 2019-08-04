#ifndef ___CLASS_MATRIXUTIL
#define ___CLASS_MATRIXUTIL
class MatrixUtil
{
	private:
	int dimension;
	double*** matrix;
	std::string filename;
	bool memory_allocated;
	
	public:
	// constructor
	MatrixUtil(double*** matrix, int dimension);
	MatrixUtil(std::string filename, double*** matrix);
	// destructor
	~MatrixUtil();
	// member
	std::string _filename() { return filename; }
	void LoadMatrix();
	void MatrixDispose();
	void MatrixReset();
	void MatrixSetUnit();
	void MatrixDivUP(double d);
	void writeMatrix(std::string filename, std::string header);

	inline int _dimension() { return dimension; }
	inline double operator[](int i) { return (*matrix)[i][i]; }

	private:
	// member
	void MemoryAllocation();
	void LoadFileData(std::ifstream* fi);
};
#endif
