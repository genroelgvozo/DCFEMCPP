#include "eigen\Eigen\Dense"
#include <iostream>
#include <fstream>
#include "DCFEM.h"
#include "FEM.h"
#include "SOLVER.h"

#include "spdlog/spdlog.h"


using namespace Eigen;
using namespace std;

double E = 3000000;
double nu = 0.16;
double h1 = 0.1, h2 = 0.1;

const int n1 = 6, n2 = 6;

int main()
{
	vd x(n1, 0), y(n2, 0);

	for (int i = 0; i < n1; i++)
		x[i] = h1 * i;
	for (int i = 0; i < n2; i++)
		y[i] = h2 * i;



	FEMPL_RECT FEMplate(n1*n2, x, y);
	FEMplate.set_mat(0, n1 - 2, 0, n2 - 2, E, nu);
	FEMplate.add_pres(0, n1 - 2, 0, n2 - 2, 10000000);
	FEMplate.constructK();

	for (int i = 0; i < n2; i++)
	{
		FEMplate.add_bc(n2*(n1-1)+i, 2);
		//FEMplate.add_bc(i, 1);
	}

	for (int i = 0; i < n1; i++)
	{
		//FEMplate.add_bc(n2*i,2);
		//FEMplate.add_bc(n2*i+n2-1, 2);
	}
	//FEMplate.print_R();



	DCFEMPL DCFEMplate(n2, 0, 0.5, y);
	//DCFEMplate.add_force(0.25, n2 / 2 + 1, 10000);
	DCFEMplate.set_mat(0, n2 - 2, E, nu);
	DCFEMplate.construct();


	SOLVER sol(DCFEMplate, FEMplate);
	sol.constructSystem();

	ofstream f("results/w.txt");
	for (int i = 0; i < n1; i++)
	{
		for (int j = 0; j < n2; j++)
		{
			f << sol.res[DCFEMplate.size_a + (i*n2 + j) * 4] << " ";
		}
		f << endl;
	}
	f.close();


	return 0;
}