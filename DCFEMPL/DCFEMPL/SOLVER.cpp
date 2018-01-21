#include "SOLVER.h"

void SOLVER::calcFundFuncMatrix(double x, int l1, MatrixXcd & res)
{
	double l2;
	if (l1 > 0)
		l2 = 0.00001;
	else
		l2 = -0.00001;
	res = T;
	for (int i = 0; i < l; i++)
	{
		if (e(i).real() * (x + l2) < 0)
			res.col(i) *= ((x + l2) > 0 ? 1.0 : -1.0)*exp(e(i)*x);
		else
			res.col(i).fill(0);
	}

	res = res * T1;

	double xi = (x + l2) > 0 ? 1 : 0;

	res += xi * (P2 + A2 * x + A2_2 * x*x / 2.0 + A2_3 *(x*x*x / 6.0) + A2_4*(x*x*x*x / 24.0) + A2_5*(x*x*x*x*x / 120.0));
}

void SOLVER::dotFundFuncMatrix(double x, int l1, const VectorXcd & C, VectorXcd & res)
{
	double l2;
	if (l1 > 0)
		l2 = 0.00001;
	else
		l2 = -0.00001;
	res = T1 * C;
	for (int i = 0; i < l; i++)
	{
		if (e(i).real() * (x + l2) < 0)
			res(i) *= ((x + l2) > 0 ? 1.0 : -1.0)*exp(e(i)*x);
		else
			res(i) = 0;
	}

	res = T * res;

	double xi = (x + l2) > 0 ? 1 : 0;
	res += xi * (P2 + A2 * x + A2_2 * x*x / 2.0 + A2_3 * x*x*x / 6.0 + A2_4*x*x*x*x / 24.0 + A2_5*x*x*x*x*x / 120.0) * C;
}

VectorXd SOLVER::calcConvolFundFunc(double x, int l)
{
	VectorXd S(size_dcfem);
	VectorXcd S1(size_dcfem);
	S.fill(0);

	for (int i = 0; i < model1.R.size(); i++)
	{
		S1.fill(0);
		dotFundFuncMatrix(x - model1.R[i].first, l, model1.R[i].second, S1);
		S += S1.real();
	}
	return S;
}

void SOLVER::constructSystem()
{
	model2.constructB();
	model1.constructB();

	K.fill(0);
	R.fill(0);
	K.block(size_dcfem + 4 * num_bound_points, size_dcfem, size_fem - 4 * num_bound_points, size_fem) = model2.K.block(4 * num_bound_points, 0, size_fem - 4 * num_bound_points, size_fem) / 1000000;
	R.segment(size_dcfem + 4 * num_bound_points, size_fem - 4 * num_bound_points) = model2.R.segment(4 * num_bound_points, size_fem - 4 * num_bound_points) / 1000000;
	MatrixXcd E(size_dcfem, size_dcfem), E2(size_dcfem, size_dcfem);

	calcFundFuncMatrix(0.0, 1, E);
	calcFundFuncMatrix(model1.x21 - model1.x22, 1, E2);
	K.block(0, 0, size_dcfem / 2, size_dcfem) = model1.B1 * (E - E2).real();
	R.segment(0, size_dcfem / 2) = -model1.B1 * calcConvolFundFunc(model1.x21, 1);

	calcFundFuncMatrix(model1.x22 - model1.x21, -1, E);
	calcFundFuncMatrix(0.0, -1, E2);
	K.block(size_dcfem / 2, 0, 8 * num_bound_points, size_dcfem) = model1.B2 * (E - E2).real();
	K.block(size_dcfem / 2, size_dcfem, size_dcfem, size_fem) = -model2.B;
	R.segment(size_dcfem / 2, size_dcfem) = -model1.B2 * calcConvolFundFunc(model1.x22, -1);

	res = K.colPivHouseholderQr().solve(R);
}

void SOLVER::compute_components()
{
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

	for (int i = 0; i < l; i++)
	{
		T.col(i).normalize();
	}


	T1 = (T1*T).inverse() *  T1;

	P2 = T * T1;


	P2 = -P2;
	for (int i = 0; i < size_dcfem; i++)
		P2(i, i) = 1.0 + P2(i, i);

	P2_2 = P2 * P2;
	A2 = P2 * model1.mat_a;
	A2_2 = A2 * A2;
	A2_3 = A2_2 * A2;
	A2_4 = A2_3 * A2;
	A2_5 = A2_4 * A2;
	A2_6 = A2_5 * A2;
	A2_6 = A2_6 * A2;
	A2_6 = A2_6 * A2;
	A2_6 = A2_6 * A2;
	A2_6 = A2_6 * A2;
}


SOLVER::SOLVER(DCFEMPL& dcfem, FEMPL_RECT& fem) : model1(dcfem), model2(fem), size_dcfem(dcfem.size_a), eigen_sol(new EigenSolver<MatrixXd>(size_dcfem)) {

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

	compute_components();

	ofstream f("results/eigen_values.txt");
	f.precision(3);
	f << fixed;

	for (int i = 0; i < e.size(); i++)
	{
		f << e(i) << " " << e2(i) << endl;
	}
	f.close();

	f.open("results/P2.txt");
	f << P2 << endl;
	f.close();

	f.open("results/P2xP2.txt");

	f << P2_2 << endl;
	f.close();

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