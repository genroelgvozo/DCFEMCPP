#pragma once
#include "DCFEM.h"
#include "FEM.h"
#include "eigen/Eigen/Dense"
#include <memory>
#include "eigen\Eigen\EigenValues"
#include "spdlog/spdlog.h"


using namespace Eigen;
using namespace std;

class SOLVER
{
	DCFEMPL& model1;
	FEMPL_RECT& model2;

	int size_dcfem;

	int l;

public:
	EigenSolver<MatrixXd>* eigen_sol;

	MatrixXcd T, T1, P2, A2, A2_2, A2_3;
	VectorXcd e;

	shared_ptr<spdlog::logger> logger;

	double eps_zero = 1e-4;

	void calcFundFuncMatrix(double x, MatrixXcd& res);
	void dotFundFuncMatrix(double x, const VectorXcd& C, VectorXcd& res);

	void constructSystem();

	SOLVER(DCFEMPL& dcfem, FEMPL_RECT& fem) : model1(dcfem), model2(fem), size_dcfem(dcfem.size_a), eigen_sol(new EigenSolver<MatrixXd>(size_dcfem)) {
		logger = spdlog::basic_logger_st("SOLVER logger", "logs/solver.txt");
		logger->info("");
		logger->info("------------------------");
		logger->info("-- SOLVER constructor --");
		logger->info("------------------------");
		logger->info("Computing eigenvalues and eigenvectors...");

		l = size_dcfem - 6;

		eigen_sol->compute(model1.mat_a);
		e = eigen_sol->eigenvalues();
		T = eigen_sol->eigenvectors();
		T.resize(size_dcfem, l);
		logger->info("T and e was computed.");

		model1.mat_a.transposeInPlace();

		logger->info("Computing A_t eigenvectors...");
		eigen_sol->compute(model1.mat_a);
		T1 = eigen_sol->eigenvectors();
		T1.resize(l, size_dcfem);
		logger->info("T1 was computed.");

		delete eigen_sol;

		model1.mat_a.transposeInPlace();


		for (int i = 0; i < l; i++)
		{
			auto c = T.col(i).dot(T1.row(i));
			c = 1.0 / c;
			T1.row(i) *= c;
		}

		P2 = T * T1;
		P2 = -P2;
		for (int i = 0; i < size_dcfem; i++)
			P2(i, i) = 1.0 + P2(i, i);
		A2 = P2 * model1.mat_a;
		A2_2 = A2 * A2;
		A2_3 = A2_2 * A2;
		logger->info("----------------------------");
		logger->info("-- END SOLVER constructor --");
		logger->info("----------------------------");
		logger->info("");
	}
};