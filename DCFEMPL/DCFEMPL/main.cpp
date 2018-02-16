#include "eigen\Eigen\Dense"
#include <iostream>
#include <fstream>
#include "DCFEM.h"
#include "FEM.h"
#include "SOLVER.h"


using namespace Eigen;
using namespace std;

double E = 3000;
double nu = 0.16;
double h1 = 5, h2 = 5;

const int n1 = 6, n2 = 2;

int main()
{
	vd x(n1, 0), y(n2, 0);

	for (int i = 0; i < n1; i++)
		x[i] = h1 * i;
	for (int i = 0; i < n2; i++)
		y[i] = h2 * i;



	FEMPL_RECT FEMplate(n1*n2, x, y);
	FEMplate.set_mat(0, n1 - 2, 0, n2 - 2, E, nu);
	//FEMplate.add_pres(0, n1 - 2, 0, n2 - 2, 10000000);
	for (int i = 0; i < n2; i++)
	{
		//DCFEMplate.add_force(0.3, i, 1000);
		FEMplate.add_force(3, i, 100);
	}

	FEMplate.constructK();

	for (int i = 0; i < n2; i++)
	{
		FEMplate.add_bc(n2*(n1 - 1) + i, 1);
		//FEMplate.add_bc(i, 1);
	}

	for (int i = 0; i < n1; i++)
	{
		//FEMplate.add_bc(n2*i,2);
		//FEMplate.add_bc(n2*i+n2-1, 2);
	}
	//FEMplate.print_R();



	DCFEMPL DCFEMplate(n2, 0, h1*(n1 - 1), y);
	for (int i = 0; i < n2; i++)
	{
		DCFEMplate.add_force(25, i, 1);
	}
	DCFEMplate.set_mat(0, n2 - 2, E, nu);
	DCFEMplate.construct();


	SOLVER sol(DCFEMplate, FEMplate);
	sol.constructSystem();

	ofstream f("results/w_fem.txt");
	f.precision(3);
	f << fixed;
	for (int i = 0; i < n1; i++)
	{
		for (int j = 0; j < n2; j++)
		{
			f << sol.res[DCFEMplate.size_a + (i*n2 + j) * 4] << " ";
		}
		f << endl;
	}
	f.close();

	f.open("results/w_dcfem.txt");
	int count = 26;
	for (int i = 0; i < count; i++)
	{
		double x2 = i * (h1*(n1 - 1)) / (count - 1);

		f << "x = " << x2 << "  \ty = ";
		VectorXd y;
		if (i == count - 1) y = sol.calcSolution(x2, -1);
		else y = sol.calcSolution(x2, 1);
		for (int j = 0; j < n2; j++)
		{
			f << y[j * 2] << " ";
		}


		f << endl;
	}
	f.close();

	f.open("results/K_fem.txt");
	f << FEMplate.K << endl;

	f.close();

	f.open("results/K_solver.txt");
	f << sol.K << endl;
	f.close();

	complex<double> c = complex<double>(1.0, 1.0);
	cout << c << " "<<exp(c);

	return 0;
}