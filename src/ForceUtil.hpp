#ifndef ___CLASS_FORCEUTIL
#define ___CLASS_FORCEUTIL
#include <string>
#include <Eigen/Core>
#include "Atom.hpp"
class ForceUtil
{
	private:
	std::string filename;
	
	public:

	ForceUtil();
	std::string _filename() { return filename; }
	bool LoadForce(std::string filename, std::vector<Eigen::Vector3d>& v);
	bool ReadLine(std::string s, Eigen::Vector3d& vtmp);
	void writeForce(std::string filename, Atom& atom, Eigen::Vector3d& f);
	void normalize(std::vector<Eigen::Vector3d>& v);

	std::vector<int> atomIndex;
};
#endif
