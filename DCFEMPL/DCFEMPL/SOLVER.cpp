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
