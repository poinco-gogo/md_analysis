#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include "ComputeSpline.hpp"
using namespace std;

ComputeSpline::ComputeSpline(vector<double>& x, vector<double>& y)
{
	for (auto& dtmp: x) this->x.push_back(dtmp);
	for (auto& dtmp: y) this->y.push_back(dtmp);
}

ComputeSpline::ComputeSpline(string filename)
{
	this -> filename = filename;
}

ComputeSpline::ComputeSpline(int numdata, double x_origin, double stepSize, double* data_value)
{
	x.resize(numdata);
	y.resize(numdata);
	for (int i = 0; i < numdata; i++)
	{
		x[i] = x_origin + i * stepSize;
		y[i] = data_value[i];
	}
}

bool ComputeSpline::load_data()
{	
	ifstream fi(filename.c_str());
	
	if (!fi)
	{
		cerr << "error: file " << filename << " not exists.\n"; 
		return false;
	}

	string s;
	while (getline(fi, s))
	{
		if (s.empty() || s.find("REMARK", 0) != string::npos)
			continue;

		double tmpx, tmpy;
		istringstream is(s);
		is >> tmpx >> tmpy;
		if (!is)
		{
			cerr << "error: bad data format!\n";
			return false;
		}
		x.push_back(tmpx);
		y.push_back(tmpy);
	}

	cout << "REMARK " << x.size() << " data read.\n";

	return true;
}

double ComputeSpline::calc_cubic_spline(double s)
{
	int n = x.size();
	int i = 0;

	// s が所属している区間はどこか
	// s <= x[i+1] になる最大のiを探す
	// データ点外の点なら最後の補間多項式を使う
	while (i < n - 2 && s > x[i + 1] )
		i++;
	
	// 区間 [x_i, x_(i+1)] の補間多項式を使う
	s -= x[i];

	return y[i] + (b[i] + (c[i] + d[i] * s) * s) * s;
}

double ComputeSpline::calc_gradient(double s)
{
	int n = x.size();
	int i = 0;
	while (i < n - 2 && s > x[i + 1])
		i++;

	s -= x[i];

	return b[i] + (2. * c[i] + 3. * d[i] * s) * s;

}

void ComputeSpline::calc_1d_coefficients()
{
	// データ点の数
	const int n = x.size();
	
	// データ点同士の間隔の数は n-1
	vector<double> h(n - 1, 0.);
	b.resize(n - 1, 0.);
	c.resize(n - 1, 0.);
	d.resize(n - 1, 0.);
	
	for (int i = 0; i < n - 1; i++)
	{
		h[i] = x[i + 1] - x[i];
		d[i] = (y[i + 1] - y[i]) / h[i];
	}
	
	// c[0] と b[0] は使わない。ノートを見るべし。
	c[1] = 2. * (h[0] + h[1]);
	b[1] = d[1] - d[0];
	
	for (int i = 2; i < n - 1; i++)
	{
		c[i] = 2. * (h[i - 1] + h[i]) - h[i - 1] * h[i - 1] / c[i - 1];
		b[i] = d[i] - d[i - 1] - b[i - 1] * h[i - 1] / c[i - 1];
	}

	vector<double> ws(n, 0.);
	
	// natural cubic spline.
	ws[0] = 0.;
	ws[n - 1] = 0.;

	// 後退代入
	ws[n - 2] = b[n - 2] / c[n - 2];
	
	for (int i = n - 3; i >= 1 ; i--)
		ws[i] = (b[i] - h[i] * ws[i + 1]) / c[i];

	// 以上でスプライン係数σは求まった。
	// 補間多項式が
	// p_i(x) = y_i + b_i(x-x_i) + c_i(x-x_i)**2 + d_i(x-x_i)**3
	// の形で使えるように、係数 b c d を計算しておく

	for (int i = 0; i < n - 1; i++)
	{
		b[i] = d[i] - h[i] * (ws[i + 1] + 2. * ws[i]);
		c[i] = 3. * ws[i];
		d[i] = (ws[i + 1] - ws[i]) / h[i];
	}
}

void ComputeSpline::calc_2d_coefficients(double* f, double* dfd1, double* dfd2, double* dfd1d2, double stepSize,  double (*w)[4])
{
 	const int Ainv[16][16] = {
	  1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	  0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0,
	 -3, 0, 0, 3, 0, 0, 0, 0,-2, 0, 0,-1, 0, 0, 0, 0,
	  2, 0, 0,-2, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0,
	  0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0,
	  0, 0, 0, 0,-3, 0, 0, 3, 0, 0, 0, 0,-2, 0, 0,-1,
	  0, 0, 0, 0, 2, 0, 0,-2, 0, 0, 0, 0, 1, 0, 0, 1,
	 -3, 3, 0, 0,-2,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	  0, 0, 0, 0, 0, 0, 0, 0,-3, 3, 0, 0,-2,-1, 0, 0,
	  9,-9, 9,-9, 6, 3,-3,-6, 6,-6,-3, 3, 4, 2, 1, 2,
	 -6, 6,-6, 6,-4,-2, 2, 4,-3, 3, 3,-3,-2,-1,-1,-2,
	  2,-2, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	  0, 0, 0, 0, 0, 0, 0, 0, 2,-2, 0, 0, 1, 1, 0, 0,
	 -6, 6,-6, 6,-3,-3, 3, 3,-4, 4, 2,-2,-2,-2,-1,-1,
	  4,-4, 4,-4, 2, 2,-2,-2, 2,-2,-2, 2, 1, 1, 1, 1,
	};

	double tmp[16];
	double sqstepSize = stepSize * stepSize;
	for (int i = 0; i < 4; i++)
	{
		tmp[i]      = f[i];
		tmp[i + 4]  = dfd1[i] * stepSize;
		tmp[i + 8]  = dfd2[i] * stepSize;
		tmp[i + 12] = dfd1d2[i] * sqstepSize;
	}

	double alpha[16];

	for (int i = 0; i < 16; i++)
	{
		double dtmp = 0.;
		for (int k = 0; k < 16; k++)
		{
			dtmp += Ainv[i][k] * tmp[k];
		}
		alpha[i] = dtmp;
	}

	int icnt = 0;
	for (int i = 0; i < 4; i++)
	for (int j = 0; j < 4; j++)
		w[i][j] = alpha[icnt++];
}
