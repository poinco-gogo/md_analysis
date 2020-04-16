#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <vector>
#include <cstdlib>
#include <cmath>
using namespace std;
int main (int argc, char ** argv)
{
	if (argc < 5)
	{
		cout << "\nCALC DIFFUSION COEFF USING MSD DATA.\n"
			"\nusage: ./a.out msd nskips ps_per_MDstep dimension('xy' or 'xyz')\n\n";
		return 1;
	}

	vector<double> _x, _y;
	ifstream fi(argv[1]);
	int nskips = atoi(argv[2]);
	string s;
	int icnt = 0;
	while (getline(fi, s))
	{
		if (s.empty() || s.find("REMARK", 0) != string::npos)
			continue;

		if (icnt++ < nskips) continue;

		istringstream is(s);
		double d1, d2;
		is >> d1 >> d2;

		_x.push_back(d1);
		_y.push_back(d2);
	}
	cout << "REMARK " << _x.size() << " data loaded.\n";

	string sdim(argv[4]);
	if (sdim == "xy")
		cout << "REMARK Calculate 2D diffusion coeff.\n";
	else if (sdim == "xyz")
		cout << "REMARK Calculate 3D diffusion coeff.\n";
	else
	{
		cerr << "\nerror: unknown argument \"" + sdim + "\"\n\n";
		return 1;
	}

	cout << "REMARK ================\n";

	double _x_ave = 0;
	double _y_ave = 0;
	for (int i = 0; i < _x.size(); i++)
	{
		_x_ave += _x[i];
		_y_ave += _y[i];
	}
	_x_ave /= _x.size();
	_y_ave /= _y.size();

	double _sxy = 0;
	double _sx  = 0;
	for (int i = 0; i < _x.size(); i++)
	{
		double delx = _x[i] - _x_ave;
		double dely = _y[i] - _y_ave;
		_sxy += delx * dely;
		_sx  += delx * delx;
	}

	double slope   = _sxy / _sx;
	double section = _y_ave - slope * _x_ave;

	cout << setprecision(8) << scientific;
	cout << "REMARK fitting to y = ax + b\n";
	cout << "REMARK a = " << setw(16) << slope << '\n';
	cout << "REMARK b = " << setw(16) << section << '\n';

	vector<double> _f;
	double _f_ave = 0;;
	for (int i = 0; i < _x.size(); i++)
	{
		double ftmp = slope * _x[i] + section;

		_f.push_back(ftmp);

		_f_ave += ftmp;
	}
	_f_ave /= _f.size();

	// correlation between predicted values and observed values
	double _syf = 0;
	double _syy = 0;
	double _sff = 0;
	for (int i = 0; i < _f.size(); i++)
	{
		double dely = _y[i] - _y_ave;
		double delf = _f[i] - _f_ave;

		_syf += dely * delf;
		_syy += dely * dely;
		_sff += delf * delf;
	}
	double r = _syf / sqrt( _syy * _sff );
	cout << "REMARK Correlation coeff. r = " << setw(16) << r << '\n';

	double step2ps = atof(argv[3]);
	double ps2sec  = 1e-12;
	double ang2cm  = 1e-8;
	double step2sec = step2ps * ps2sec;

	// slope [ang * ang / step]
	double factor = sdim == "xy" ? 1. / 4. : 1. / 6.;
	double coeff = slope * ang2cm * ang2cm / step2sec * factor;

	cout << "REMARK Diffusion coeff.   D = " << 
		setw(16) << coeff << " [cm2/s]\n";
}
