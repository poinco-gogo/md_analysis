#ifndef ___CLASS_NAMDBIN
#define ___CLASS_NAMDBIN

#include <string>
#include <vector>
#include "Atom.hpp"

class NAMDBin
{
	public:
	NAMDBin(std::string filename, std::string type);
	bool read_fi (std::vector<Atom>& atomVector);
	bool write_fo(std::vector<Atom>& atomVector);

	private:
	std::string filename;
	std::string type;
	void read_coor(std::ifstream& fi, std::vector<Atom>& atomVector);
	void read_vel (std::ifstream& fi, std::vector<Atom>& atomVector);
	void write_coor(std::ofstream& fo, std::vector<Atom>& atomVector);
	void write_vel (std::ofstream& fo, std::vector<Atom>& atomVector);
};
#endif
