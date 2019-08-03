#ifndef ___CLASS_COMMON
#define ___CLASS_COMMON

#include <Eigen/Core>
#include <cmath>

const Eigen::Vector3d V3ZERO = Eigen::Vector3d::Zero();
const Eigen::Matrix3d M3ZERO = Eigen::Matrix3d::Zero();
const Eigen::Matrix3d M3IDEN = Eigen::Matrix3d::Identity();

const double N_PI     = 3.14159265358979323846; // from NAMD
const double COULOMB  = 332.0636;    // from NAMD.
const double BOLTZMAN = 0.001987191; // from NAMD, in kcal/mol/K.
const double PS2ASU   = 20.455;

// derived constants
const double INVBOLTZMAN = 1. / BOLTZMAN;
const double N_TWOPI = 2. * N_PI;
const double N_SQRTPI= sqrt(N_PI);
const double DEG2RAD = N_PI / 180.;
const double RAD2DEG = 180. / N_PI;
const double BAR2CAL = 6.022141 / 4.184 * 1e-5;
const double CAL2BAR = 1. / BAR2CAL;

inline void err(std::string error_message)
{
	std::cerr << "\nerror: " << error_message << "\n\n";
}

inline void die(std::string error_message)
{
	std::cerr << '\n' << error_message << "\n\n";
	exit(1);
}
#endif
