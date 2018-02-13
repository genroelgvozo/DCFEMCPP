#include "SOLVER.h"

void SOLVER::calcFundFuncMatrix(double x, int l1, MatrixXcd & res)
{
	double l2;
	if (l1 > 0)
		l2 = 0.001;
	else
		l2 = -0.001;
	MatrixXcd temp(size_dcfem, l);
	temp = T;
	for (int i = 0; i < l; i++)
	{
		if ((e(i).real() < 0.00001 && (x + l2) > 0) ||
			(e(i).real() > 0 && (x + l2) < 0))
		{
			temp.col(i) *= ((x + l2) > 0 ? 1.0 : -1.0)*exp(e(i)*x);
		}
		else
			temp.col(i).fill(0);
	}

	res = temp * T1;

	if (l < size_dcfem)
		if (x + l2 > 0) {
			res += (P2 + A2 * x + A2_2 * x*x / 2.0 + A2_3 *(x*x*x / 6.0)) / 2;
		}
		else
		{
			res -= (P2 + A2 * x + A2_2 * x*x / 2.0 + A2_3 *(x*x*x / 6.0)) / 2;
		}


}

void SOLVER::dotFundFuncMatrix(double x, int l1, const VectorXd & C1, VectorXcd & res)
{
	double l2;
	if (l1 > 0)
		l2 = 0.001;
	else
		l2 = -0.001;
	VectorXcd temp(l);
	temp = T1*C1;
	for (int i = 0; i < l; i++)
	{
		if ((e(i).real() < 0.000001 && (x + l2) > 0) || (e(i).real() > 0.000001 && (x + l2) < 0))
			temp(i) *= ((x + l2) > 0 ? 1.0 : -1.0)*exp(e(i)*x);
		else
			temp(i) = 0;
	}

	res = T * temp;
	if (l < size_dcfem)
		if (x + l2 > 0) {
			res += (P2 + A2 * x + A2_2 * x*x / 2.0 + A2_3 *(x*x*x / 6.0))*C1 / 2;
		}
		else
		{
			res -= (P2 + A2 * x + A2_2 * x*x / 2.0 + A2_3 *(x*x*x / 6.0))*C1 / 2;
		}
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

VectorXd SOLVER::calcSolution(double x, int l)
{
	VectorXd y(size_dcfem);
	VectorXcd temp(size_dcfem);

	MatrixXcd E(size_dcfem, size_dcfem), E2(size_dcfem, size_dcfem);

	//dotFundFuncMatrix(x - model1.x21, l, C, temp);
	calcFundFuncMatrix(x - model1.x21, l, E);

	//y = temp.real();
	//dotFundFuncMatrix(x - model1.x22, l, C, temp);
	calcFundFuncMatrix(x - model1.x22, l, E2);
	temp = (E - E2)*C;
	y = temp.real();
	y += calcConvolFundFunc(x, l);
	return y;
}

void SOLVER::constructSystem()
{
	model2.constructB();
	model1.constructB();

	K.fill(0);
	R.fill(0);
	MatrixXcd E(size_dcfem, size_dcfem), E2(size_dcfem, size_dcfem);
	VectorXd rr(size_dcfem / 2);
	rr.fill(0);
	for (int i = 0; i < model1.num_points; i++)
		rr[i * 4] = 1;
	calcFundFuncMatrix(0.0, 1, E);
	calcFundFuncMatrix(model1.x21 - model1.x22, 1, E2);
	E = E - E2;
	K.block(0, 0, size_dcfem / 2, size_dcfem) = model1.B1 * E.real();
	R.segment(0, size_dcfem / 2) = rr - model1.B1 * calcConvolFundFunc(model1.x21, 1);

	E.fill(0);
	E2.fill(0);
	calcFundFuncMatrix(model1.x22 - model1.x21, -1, E);
	calcFundFuncMatrix(0.0, -1, E2);
	E = E - E2;
	K.block(size_dcfem / 2, 0, size_dcfem, size_dcfem) = model1.B2 * E.real();
	K.block(size_dcfem / 2, size_dcfem, size_dcfem, size_fem) = -model2.B;
	R.segment(size_dcfem / 2, size_dcfem) = -model1.B2 * calcConvolFundFunc(model1.x22, -1);


	ofstream f("results/R_solver.txt");
	f.precision(3);
	f << fixed;
	f << R;
	f.close();

	K.block(size_dcfem + size_dcfem / 2, size_dcfem, size_fem - 4 * num_bound_points, size_fem) =
		model2.K.block(4 * num_bound_points, 0, size_fem - 4 * num_bound_points, size_fem);
	R.segment(size_dcfem + size_dcfem / 2, size_fem - 4 * num_bound_points) = model2.R.segment(4 * num_bound_points, size_fem - 4 * num_bound_points);

	res = K.colPivHouseholderQr().solve(R);

	//C = res.segment(0, size_dcfem);
	C = K.block(0, 0, size_dcfem, size_dcfem).colPivHouseholderQr().solve(R.segment(0, size_dcfem));
}

void SOLVER::compute_components()
{
	VectorXcd tmp;

	for (int i = l; i > 0; i--)
		for (int j = 0; j < i - 1; j++)
		{
			if (abs(e(j).real()) < abs(e(j + 1).real()))
			{
				tmp = T.col(j);
				T.col(j) = T.col(j + 1);
				T.col(j + 1) = tmp;
				swap(e(j), e(j + 1));
			}
		}
	int first_zero = 0;
	for (int i = 0; i < l; i++)
		if (abs(e(i).real()) < 0.001)
		{
			first_zero = i;
			break;
		}
	//for (int i = first_zero; i < l; i++)
		//e(i) = complex<double>(0.0, e(i).imag());

	for (int i = first_zero; i > 0; i--)
		for (int j = 0; j < i - 1; j++)
		{
			if (e(j).real() < e(j + 1).real())
			{
				tmp = T.col(j);
				T.col(j) = T.col(j + 1);
				T.col(j + 1) = tmp;
				swap(e(j), e(j + 1));
			}
		}

	for (int i = l; i > first_zero; i--)
	{
		for (int j = first_zero; j < i - 1; j++)
			if (e(j).imag() < e(j + 1).imag())
			{
				tmp = T.col(j);
				T.col(j) = T.col(j + 1);
				T.col(j + 1) = tmp;
				swap(e(j), e(j + 1));
			}
	}
	for (int i = 0; i < l; i++)
	{
		for (int j = i; j < l; j++)
			if (abs(e2(j).real() - e(i).real()) < 0.001 && abs(e2(j).imag() - e(i).imag()) < 0.001)
			{
				tmp = T1.row(i);
				T1.row(i) = T1.row(j);
				T1.row(j) = tmp;
				swap(e2(j), e2(i));
				break;
			}
	}


	T1 = (T1*T).inverse() *  T1;

	P1 = T * T1;
	P1_2 = P1*P1;

	A1 = P1*model1.mat_a;

	P2 = -P1;
	for (int i = 0; i < size_dcfem; i++)
		P2(i, i) = 1.0 + P2(i, i);

	P2_2 = P2 * P2;
	A2 = P2 * model1.mat_a;


	//A2 = T;
	//for (int i = 0; i < l; i++)
	//	A2.col(i) *= e(i);
	//A2 = A2*T1;
	//A2 = model1.mat_a - A2;
	A2_2 = A2 * A2;
	A2_3 = A2_2 * A2;
	A2_4 = A2_3 * A2;
	A2_5 = A2_4 * A2;
	A2_6 = A2_5 * A2;


}


SOLVER::SOLVER(DCFEMPL& dcfem, FEMPL_RECT& fem) : model1(dcfem), model2(fem), size_dcfem(dcfem.size_a), eigen_sol(EigenSolver<MatrixXd>(size_dcfem)) {


	//l = size_dcfem - 6;
	l = size_dcfem ;

	eigen_sol.compute(model1.mat_a);
	e = eigen_sol.eigenvalues();

	T = eigen_sol.eigenvectors();
	T.resize(size_dcfem, l);


	model1.mat_a.transposeInPlace();

	eigen_sol.compute(model1.mat_a);
	e2 = eigen_sol.eigenvalues();


	T1 = eigen_sol.eigenvectors();
	T1.resize(l, size_dcfem);


	model1.mat_a.transposeInPlace();

	compute_components();

	ofstream f("results/eigen_values.txt");
	f.precision(4);
	f << fixed;

	for (int i = 0; i < e.size(); i++)
	{
		f << e(i) << "\t" << e2(i) << endl;
	}
	f.close();

	f.precision(3);
	f.open("results/T.txt");
	f << T << endl;
	f.close();

	f.open("results/T1.txt");
	f << T1 << endl;
	f.close();

	f.open("results/P2.txt");
	f << P2 << endl;
	f.close();

	f.open("results/P2xP2.txt");

	f << P2_2 << endl;
	f.close();

	f.open("results/P1.txt");
	f << P1 << endl;
	f.close();

	f.open("results/P1xP1.txt");

	f << P1_2 << endl;
	f.close();

	f.open("results/A1.txt");
	f << A1 << endl;
	f.close();

	f.open("results/AP1.txt");
	f << model1.mat_a * P1 << endl;
	f.close();
	f.open("results/A1A2.txt");
	f << (A1*A2).eval() << endl;
	f.close();
	f.open("results/A2A1.txt");
	f << (A2*A1).eval() << endl;
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

	MatrixXcd E(size_dcfem, size_dcfem), E2(size_dcfem, size_dcfem);
	calcFundFuncMatrix(0, 1, E);
	calcFundFuncMatrix(0, -1, E2);
	f.open("results/E-E2.txt");
	f << E - E2 << endl;
	f.close();
	f.open("results/P1+P2.txt");
	f << P1 + P2 << endl;
	f.close();


	f.open("results/mat_a2.txt");
	f.precision(5);
	f << model1.mat_a << endl;
	f.close();

	size_fem = model2.get_size();
	num_bound_points = model1.num_points;
	size_k = size_dcfem + size_fem;
	K.resize(size_k, size_k);
	R.resize(size_k);

}