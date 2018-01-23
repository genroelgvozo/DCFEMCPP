#pragma once
#include "DCFEM.h"
#include "FEM.h"
#include "eigen/Eigen/Dense"
#include <memory>
#include "eigen\Eigen\EigenValues"
#include "spdlog/spdlog.h"
#include <fstream>


using namespace Eigen;
using namespace std;

class SOLVER
{
	DCFEMPL& model1;
	FEMPL_RECT& model2;

	int size_dcfem;
	int size_fem;
	int num_bound_points;

	int num_dcfem;
	int num_fem_x;
	int num_fem_y;

	int l;

	int size_k;


public:
	MatrixXd K;
	VectorXd R;
	VectorXd res;
	EigenSolver<MatrixXd>* eigen_sol;

	MatrixXcd T, T1, P1, P2, A1, A2, A2_2, A2_3, A2_4, A2_5, A2_6, P1_2, P2_2;
	VectorXcd e, e2;

	shared_ptr<spdlog::logger> logger;

	double eps_zero = 1e-4;

	void calcFundFuncMatrix(double x, int l1, MatrixXcd& res);
	void dotFundFuncMatrix(double x, int l1, const VectorXcd& C, VectorXcd& res);

	VectorXd calcConvolFundFunc(double x, int l);

	void constructSystem();

	void compute_components();

	SOLVER(DCFEMPL& dcfem, FEMPL_RECT& fem);
};