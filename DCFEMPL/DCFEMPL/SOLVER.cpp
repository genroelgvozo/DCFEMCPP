#include "SOLVER.h"

void SOLVER::calcFundFuncMatrix(double x, MatrixXcd & res)
{
	res = T;
	for(int i=0;i<l;i++)
	{
		if (e(i).real() * x < 0)
			res.col(i) *= (x > 0 ? 1.0 : -1.0)*exp(e(i)*x);
		else
			res.col(i).fill(0);
	}

	res = res * T1;

	double xi = x > 0 ? 0.5 : -0.5;

	res += xi * (P2 + A2*x + A2_2 * x*x / 2.0 + A2_3 * x*x*x / 6.0);
}

void SOLVER::dotFundFuncMatrix(double x, const VectorXcd & C, VectorXcd & res)
{
	res = T1 * C;
	for (int i = 0; i < l; i++)
	{
		if (e(i).real() * x < 0)
			res(i) *= (x > 0 ? 1.0 : -1.0)*exp(e(i)*x);
		else
			res(i) = 0;
	}

	res = T * res;

	double xi = x > 0 ? 0.5 : -0.5;
	res += xi * (P2 + A2 * x + A2_2 * x*x / 2.0 + A2_3 * x*x*x / 6.0) * C;
}

void SOLVER::constructSystem()
{
	model2.constructB();

}
