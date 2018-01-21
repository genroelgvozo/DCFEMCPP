#pragma once
#include "eigen/Eigen/Dense"
#include <vector>
#include <memory>
#include <fstream>


using namespace Eigen;
using namespace std;

class DCFEM
{
public:


	DCFEM(int _num_points, int _num_dofs) :num_points(_num_points), num_dofs(_num_dofs), size_a(2 * _num_points*num_dofs), mat_a(MatrixXd(size_a, size_a)) {
		mat_a.fill(0);
	};

	int num_points;
	int num_dofs;
	int size_a;

	MatrixXd mat_a;
};

struct load {
	int type;
	double x2;
	int num_node;
	double value;
};

class DCFEMPL : public DCFEM
{

	Matrix4d get_K0(int i);
	Matrix4d get_K2(int i);
	Matrix4d get_K4(int i);

public:
	double Ni(int i, double x, int s1, double h);
	MatrixXd B1, B2;

	void constructB();
	void construct();
	DCFEMPL(int num_points, double _x21, double _x22, const vector<double>& x) : DCFEM(num_points, 4), xmesh(x), EX(vector<double>(num_points - 1, 0)), PRXY(vector<double>(num_points - 1, 0)) {
		x21 = _x21;
		x22 = _x22;
		x_size = xmesh.size();
	};




	void set_mat(int x1, int x2, double _EX, double _PRXY)
	{
		if (x1 < 0 || x2 >= x_size - 1)
			return;
		for (int i = x1; i <= x2; ++i)
		{
			EX[i] = _EX;
			PRXY[i] = _PRXY;
		}
	}

	void add_force(double x2, int num, double f)
	{
		load l;
		l.num_node = num;
		l.type = 1;
		l.value = f;
		l.x2 = x2;
		loads.push_back(l);
	}


	vector<double> xmesh;
	int x_size;
	double x21, x22;

	vector< double> EX, PRXY;
	double dh = 0.1;
	vector<load> loads;

	vector<pair<double, VectorXd >> R;


};