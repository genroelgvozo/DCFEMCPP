#pragma once
#include "FEM_base.h"
#include <math.h>
#include <iostream>


using vi = vector<int>;
using vd = vector<double>;
using vvi = vector<vi>;
using vvd = vector<vd>;

using namespace std;
using namespace Eigen;



class FEMPL_RECT : public FEM_base {
public:

	MatrixXd B;

	double Ni(int i, double x, int s1);
	double Nij(int i, int j, double x, double y, int s1, int s2);


	MatrixXd get_Kij(int i, int j, double a, double b, const Matrix3d& D);
	MatrixXd get_Klocal1(int i, int j);
	MatrixXd get_Klocal(int i, int j);
	VectorXd get_Rlocal(int i, int j);
	MatrixXd get_Bi(int i, double x, double y, double a, double b);

	FEMPL_RECT(int num_points, const vd& x, const vd& y) : FEM_base(num_points, 4),
		xmesh(x), ymesh(y), EX(vvd(x.size() - 1, vd(y.size() - 1, 0))), PRXY(vvd(x.size() - 1, vd(y.size() - 1, 0))), press(vvd(x.size() - 1, vd(y.size() - 1, 0))) {

		x_size = xmesh.size();
		y_size = ymesh.size();
		B.resize(8 * y_size, 4 * num_points);
	};

	void calc_elem_moments(VectorXd u, int i, int j, double x, double y, VectorXd& M);

	void calc_elem_result(VectorXd u, int i, int j, int num, VectorXd& M, VectorXd& V);

	void calc_node_result(VectorXd u, int i, int j, VectorXd& M, VectorXd& V);

	void set_mat(int x1, int x2, int y1, int y2, double _EX, double _PRXY)
	{
		if (x1 < 0 || x2 >= x_size - 1)
			return;
		if (y1 < 0 || y2 >= y_size - 1)
			return;
		for (int i = x1; i <= x2; ++i)
			for (int j = y1; j <= y2; ++j)
			{
				EX[i][j] = _EX;
				PRXY[i][j] = _PRXY;
			}
	}

	void add_pres(int x1, int x2, int y1, int y2, double _press)
	{
		if (x1 < 0 || x2 >= x_size - 1)
			return;
		if (y1 < 0 || y2 >= y_size - 1)
			return;
		for (int i = x1; i <= x2; ++i)
			for (int j = y1; j <= y2; ++j)
				press[i][j] = _press;
	}

	void constructB();

	void constructK();

	void add_bc(int num_node, int type);

	void print_R()
	{
		cout << R << endl;
	}

private:
	vector<double> xmesh;
	vector<double> ymesh;

	vector< vector<double>> EX, PRXY;
	vector < vector<double>> press;
	size_t x_size;
	size_t y_size;

	double h = 0.1;
};