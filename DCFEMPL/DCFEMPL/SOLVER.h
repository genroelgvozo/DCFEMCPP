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

	MatrixXcd T, T1, P2, A2, A2_2, A2_3, A2_4, A2_5, A2_6, P2_2;
	VectorXcd e, e2;

	shared_ptr<spdlog::logger> logger;

	double eps_zero = 1e-4;

	void calcFundFuncMatrix(double x, int l1, MatrixXcd& res);
	void dotFundFuncMatrix(double x, int l1, const VectorXcd& C, VectorXcd& res);

	VectorXd calcConvolFundFunc(double x, int l);

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
		e2 = eigen_sol->eigenvalues();

		T = eigen_sol->eigenvectors();
		T.resize(size_dcfem, l);
		logger->info("T and e was computed.");

		model1.mat_a.transposeInPlace();

		logger->info("Computing A_t eigenvectors...");
		eigen_sol->compute(model1.mat_a);
		e = eigen_sol->eigenvalues();


		T1 = eigen_sol->eigenvectors();
		T1.resize(l, size_dcfem);
		logger->info("T1 was computed.");

		delete eigen_sol;

		model1.mat_a.transposeInPlace();
		VectorXcd tmp;
		for (int i = 0; i < l; i++)
		{
			for (int j = 0; j < l; j++)
				if (abs(e2(j) - e(i)) < 0.001)
				{
					tmp = T1.row(i);
					T1.row(i) = T1.row(j);
					T1.row(j) = tmp;
					swap(e2(j), e2(i));
					break;
				}
		}

		ofstream f("results/eigen_values.txt");
		f.precision(3);
		f << fixed;

		for (int i = 0; i < e.size(); i++)
		{
			f << e(i) << " " << e2(i) << endl;
		}
		f.close();

		for (int i = 0; i < l; i++)
		{
			T.col(i).normalize();
		}

		for (int i = 0; i < l; i++)
		{
			auto c = T.col(i).dot(T1.row(i));
			c = 1.0 / c;
			T1.row(i) *= c;
		}


		T1 = (T1*T).inverse() *  T1;

		P2 = T * T1;


		P2 = -P2;
		for (int i = 0; i < size_dcfem; i++)
			P2(i, i) = 1.0 + P2(i, i);


		f.open("results/P2.txt");
		f << P2 << endl;
		f.close();

		f.open("results/P2xP2.txt");
		P2_2 = P2 * P2;
		f << P2_2 << endl;
		f.close();


		A2 = P2 * model1.mat_a;
		A2_2 = A2 * A2;
		A2_3 = A2_2 * A2;
		A2_4 = A2_3 * A2;
		A2_5 = A2_4 * A2;
		A2_6 = A2_5 * A2;

		f.open("results/A2.txt");
		f << A2 << endl;
		f.close();
		f.open("results/A2_2.txt");
		f << A2_2 << endl;
		f.close();
		f.open("results/A2_3.txt");
		f << A2_3 << endl;
		f.close();
		f.open("results/A2_4.txt");
		f << A2_4 << endl;
		f.close();
		f.open("results/A2_5.txt");
		f << A2_5 << endl;
		f.close();
		f.open("results/A2_6.txt");
		f << A2_6 << endl;
		f.close();

		MatrixXcd E, E2;
		calcFundFuncMatrix(0, 1, E);
		calcFundFuncMatrix(0, -1, E2);
		f.open("results/E-E2.txt");
		f << E - E2 << endl;
		f.close();

		size_fem = model2.get_size();
		num_bound_points = model1.num_points;
		size_k = size_dcfem + size_fem;
		K.resize(size_k, size_k);
		R.resize(size_k);

		logger->info("----------------------------");
		logger->info("-- END SOLVER constructor --");
		logger->info("----------------------------");
		logger->info("");
	}
};