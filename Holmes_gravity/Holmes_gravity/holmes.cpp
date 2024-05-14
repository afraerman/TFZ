#define _USE_MATH_DEFINES

double EARTH_RADIUS = 6378136.3; // m
double GM = 398600441500000.0; // m^3 /sec^2
double ALPHA = 1.0 / 298.3;

const int M = 60;

int** Order2 = new int*[M + 1];
double* Cnm = new double[(M + 2)*(M + 1) / 2];
double* Snm = new double[(M + 2)*(M + 1) / 2];

double** Ps = new double*[M + 1];
double** ml = new double*[2];

#include<iostream>
#include<string>
#include<vector>
#include<iomanip>
#include<fstream>
#include<algorithm>
#include<cmath>
#include<gravity.h>
#include<stdlib.h>

#define NUMBER_OF_FILES  66

using namespace std;

void create_map(string, string);
void get_files(string, vector<string>&, vector<string>&);

int main()
{
	string files_list = "C:/Users/Алексей/Documents/10_семестр/теория фигуры Земли/files_after_matlab.txt";
	vector<string> egmfiles(NUMBER_OF_FILES);
	vector<string> savefiles(NUMBER_OF_FILES);
	get_files(files_list, egmfiles, savefiles);


	for (int i = 0; i < NUMBER_OF_FILES; i++)
	{
		cout << "Working on file " << i << '\n';
		cout << egmfiles[i] << '\n';
		cout << savefiles[i] << '\n';
		create_map(egmfiles[i], savefiles[i]);
	}

	delete[] Ps; delete[] Order2; delete[] ml;
	system("pause");
	return 0;
}

void create_map(string egmfilename, string savefilename)
{
	// fill gravity coefficients arrays Cnm and Snm
	gravity::get_gravity_coefficients(egmfilename);
	std::ofstream fout(savefilename);
	fout << std::setprecision(17);

	double r, ctheta, stheta, slambda, clambda;
	double* gitrf = new double[3];

	for (double theta = M_PI_2 - 1.4; theta < M_PI_2 + 1.4; theta += 0.017)
	{
		ctheta = cos(theta);
		stheta = sin(theta);
		//r = EARTH_RADIUS * (1.0 - ALPHA) / sqrt((1.0 - ALPHA)*(1.0 - ALPHA)*stheta*stheta + ctheta * ctheta);
		for (double lambda = 0.0; lambda < 2.0 * M_PI; lambda += 0.017)
		{
			// create arrays of Ps
			for (int i = 0; i < M + 1; i++) Ps[i] = new double[M + 1];

			// create array for sin(ml) and cos(ml)
			ml[0] = new double[M + 1];
			ml[1] = new double[M + 1];

			// gravity

			slambda = sin(lambda);
			clambda = cos(lambda);

			gravity::allPs(ctheta);
			gravity::allml(slambda, clambda);
			gravity::earth_gravity(EARTH_RADIUS, stheta, ctheta, slambda, clambda, M, gitrf);

			fout << theta << '\t' << lambda << '\t' << sqrt(gitrf[0] * gitrf[0] + gitrf[1] * gitrf[1] + gitrf[2] * gitrf[2]) << '\n';

			delete[] ml[0]; delete[] ml[1];
			for (int j = 0; j < M + 1; j++)
			{
				delete[] Ps[j];
			}
		}
	}

	fout.close();

	for (int i = 0; i < M + 1; i++)
	{
		delete[] Order2[i];
	}
	delete[] gitrf;
}

void get_files(string files_list, vector<string>& egmfiles, vector<string>& savefiles)
{
	int i = 0;
	string egmfile, savefile;
	std::ifstream f(files_list);
	while (getline(f, egmfile))
	{
		getline(f, savefile);
		//std::cout << egmfile << '\n';
		//system("pause");
		egmfiles[i] = egmfile;
		savefiles[i] = savefile;
		i++;
	}
	f.close();
}