#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <vector>
#include <ctime>
using namespace std;
int main(int argc, char** argv)
{
	if (argc < 2)
	{
		cout << "\nusage: ./a.out mat1 mat2 ... denom ofilename\n\n";
		return 1;
	}

	cout << "REMARK ";
	for (int i = 0; i < argc; i++)
		cout << argv[i] << " " ;
	cout << '\n';

	istringstream iss(argv[argc - 2]);
	double denom = 0.;
	iss >> denom;
	if (denom <= 0.)
	{
		cerr << "error: ./getmatsum.x : denom cannot be " 
			<< argv[argc - 2] << '\n';
		return 0;
	}

	int ndim = 0;

	vector<double> vele;
	
	for (int i = 0; i < argc - 3; i++)
	{
		ifstream fi(argv[i + 1]);
		if (!fi)
		{
			cerr << "error: ./getmatsum.x : file " << 
				argv[i + 1] << " not exsits\n";
			return 0;
		}
		istringstream is;
		string s;
		double d1 = 0.;
		double d2 = 0.;
		double d3 = 0.;
		int nline = 0;
		while (getline(fi, s))
		{
			if (s.empty() || s.find("REMARK", 0) != string::npos)
				continue;
			is.clear();
			if (s.find("DIMENSION", 0) != string::npos)
			{
				is.str(s);
				is >> ndim;
				cout << "REMARK matrix " << argv[i + 1]
					<< " has " << ndim << " dimension.\n";
				if (vele.size() == 0)
				{
					vele.resize(ndim * (ndim + 1) / 2, 0.);
					cout << "REMARK vector resized size : "
						<< vele.size() << endl;
				}
				continue;
			}
			
			nline += 1;

			is.str(s);
			is >> d1 >> d2 >> d3;
			vele.at(3 * nline - 3) += d1;
			vele.at(3 * nline - 2) += d2;
			vele.at(3 * nline - 1) += d3;
		
		}// end while getline.
		
	}// end read all files.
	
	cout << "REMARK finally divided by " << denom << '\n';
	denom = 1. / denom;
	for (int i = 0; i < vele.size(); i++)
		vele[i] *= denom;
	
	cout << "REMARK output to file " << argv[argc - 1] << '\n';
	ofstream fo(argv[argc - 1]);
	time_t now = time(NULL);
	struct tm* local = localtime(&now);
	fo << "REMARK " << asctime(local);
	fo << "REMARK sum of following matrices: \n";
	for (int i = 0; i < argc - 3; i++)
       		fo << "REMARK " << argv[i + 1] << '\n';
	fo << "REMARK the sum is divided by command line argument = "
	 << 1. / denom << '\n';	

	fo << '\n';
	fo << setw(8) << ndim << " !DIMENSION\n";
	for (int i = 0; i < vele.size() / 3; i++)
	{
		fo << scientific << setprecision(8) 
			<< setw(16) << vele[3 * i] 
			<< setw(16) << vele[3 * i + 1] 
			<< setw(16) << vele[3 * i + 2] << '\n';
	}
}
