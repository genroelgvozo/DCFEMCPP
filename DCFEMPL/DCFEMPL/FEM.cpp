#include "FEM.h"
#include <iostream>
#include <fstream>

double FEMPL_RECT::Ni(int i, double x, int s1)
{
	switch (s1)
	{
	case 0:
		switch (i)
		{
		case 1:
			return 1 - 3 * x*x + 2 * x*x*x;
		case 2:
			return x - 2 * x*x + x * x*x;
		case 3:
			return 3 * x*x - 2 * x*x*x;
		case 4:
			return x * x*x - x * x;
		default:
			break;
		}
		break;
	case 1:
		switch (i)
		{
		case 1:
			return -6 * x + 6 * x*x;
		case 2:
			return 1 - 4 * x + 3 * x*x;
		case 3:
			return 6 * x - 6 * x*x;
		case 4:
			return 3 * x*x - 2 * x;
		default:
			break;
		}
		break;
	case 2:
		switch (i)
		{
		case 1:
			return 12 * x - 6;
		case 2:
			return 6 * x - 4;
		case 3:
			return 6 - 12 * x;
		case 4:
			return 6 * x - 2;
		default:
			break;
		}
		break;
	case 3:
		switch (i)
		{
		case 1:
			return 12;
		case 2:
			return 6;
		case 3:
			return -12;
		case 4:
			return 6;
		default:
			break;
		}
	default:
		break;
	}
	return 0.0;
}

double FEMPL_RECT::Nij(int i, int j, double x, double y, int s1, int s2)
{
	return Ni(i, x, s1)*Ni(j, y, s2);
}

MatrixXd FEMPL_RECT::get_Kij(int i, int j, double a, double b, const Matrix3d & D)
{
	double x = 1.0 / sqrt(3.0);
	double weight = 1;
	vector<pair<double, double>> points(4);
	points[0] = make_pair(0.5 - x / 2, 0.5 - x / 2);
	points[1] = make_pair(0.5 + x / 2, 0.5 - x / 2);
	points[2] = make_pair(0.5 - x / 2, 0.5 + x / 2);
	points[3] = make_pair(0.5 + x / 2, 0.5 + x / 2);


	MatrixXd Bi(3, 4), Bj(3, 4);
	MatrixXd Kij(4, 4);
	Kij.fill(0);
	for (int num = 0; num < points.size(); num++)
	{
		Bi = get_Bi(i, points[num].first, points[num].second, a, b);
		Bj = get_Bi(j, points[num].first, points[num].second, a, b);
		Kij += Bi.transpose() * D * Bj * weight;
	}
	Kij = Kij * a*b;
	return Kij;
}

MatrixXd FEMPL_RECT::get_Klocal1(int i, int j)
{
	double a, b;
	a = xmesh[i + 1] - xmesh[i];
	b = ymesh[j + 1] - ymesh[j];
	double E, nu;
	E = EX[i][j];
	nu = PRXY[i][j];
	Matrix3d D;
	D(0, 0) = E * h*h*h / (12 * (1 - nu * nu));
	D(1, 1) = D(0, 0);
	D(0, 1) = D(1, 0) = D(0, 0) * nu;
	D(2, 2) = D(0, 0) * (1 - nu) / 2;
	MatrixXd Kl(16, 16);

	for (int i = 0; i < 4; i++)
		for (int j = 0; j < 4; j++)
		{
			Kl.block(i * 4, j * 4, 4, 4) = get_Kij(i + 1, j + 1, a, b, D);
		}
	return Kl;
}

MatrixXd FEMPL_RECT::get_Klocal(int i, int j)
{
	double a, b;
	a = xmesh[i + 1] - xmesh[i];
	b = ymesh[j + 1] - ymesh[j];
	double E, nu;
	E = EX[i][j];
	nu = PRXY[i][j];
	double Dx, Dy, D1, Dxy;
	Dx = E * h*h*h / (12 * (1 - nu * nu));
	Dy = Dx;
	D1 = Dx * nu;
	Dxy = Dx * (1 - nu) / 2;
	MatrixXd Kl(16, 16);
	Kl.fill(0);
	Kl(0, 0) = Kl(4, 4) = Kl(8, 8) = Kl(12, 12) = 156 * b * Dx / (35 * a*a*a) + 72 * D1 / (25 * a*b) + 156 * a*Dy / (35 * b*b*b) + 144 * Dxy / (25 * a*b);

	Kl(0, 1) = Kl(12, 13) = 78 * a*Dy / (35 * b*b*b) + 36 * D1 / (25 * a*b) + 22 * b*Dx / (35 * a*a*a) + 12 * Dxy / (25 * a*b);
	Kl(4, 5) = Kl(8, 9) = -Kl(0, 1);

	Kl(0, 2) = Kl(4, 6) = Kl(12, 14) = 22 * a*Dy / (35 * b*b*b) + 36 * D1 / (25 * a*b) + 78 * b*Dx / (35 * a*a*a) + 12 * Dxy / (25 * a*b);
	Kl(8, 10) = -Kl(0, 2);

	Kl(0, 3) = Kl(8, 11) = 11 * b*Dx / (35 * a*a*a) + 11 * D1 / (50 * a*b) + 11 * a*Dy / (35 * b*b*b) + Dxy / (25 * a*b);
	Kl(4, 7) = Kl(12, 15) = -Kl(0, 3);

	Kl(0, 4) = Kl(8, 12) = -156 * b*Dx / (36 * a*a*a) - 72 * D1 / (25 * a*b) + 54 * a*Dy / (35 * b*b*b) - 144 * Dxy / (25 * a*b);

	Kl(0, 5) = Kl(9, 12) = -22 * b*Dx / (35 * a*a*a) - 36 * D1 / (25 * a*b) + 27 * a*Dy / (35 * b*b*b) - 12 * Dxy / (25 * a*b);
	Kl(1, 4) = Kl(8, 13) = -Kl(0, 5);

	Kl(0, 6) = Kl(2, 4) = 78 * b*Dx / (35 * a*a*a) + 6 * D1 / (25 * a*b) - 13 * a*Dy / (25 * b*b*b) + 12 * Dxy / (25 * a*b);
	Kl(8, 14) = Kl(10, 12) = -Kl(0, 6);

	Kl(0, 7) = Kl(1, 6) = Kl(8, 15) = Kl(9, 14) = 11 * b*Dx / (35 * a*a*a) + 3 * D1 / (25 * a*b) - 13 * a*Dy / (70 * b*b*b) + Dxy / (25 * a*b);
	Kl(2, 5) = Kl(3, 4) = Kl(11, 12) = Kl(10, 13) = -Kl(0, 7);

	Kl(0, 8) = Kl(4, 12) = -54 * b*Dx / (35 * a*a*a) + 72 * D1 / (25 * a*b) - 54 * a*Dy / (35 * b*b*b) + 144 * Dxy / (25 * a*b);

	Kl(0, 9) = Kl(5, 12) = 13 * b*Dx / (35 * a*a*a) - 6 * D1 / (25 * a*b) + 27 * a*Dy / (35 * b*b*b) - 12 * Dxy / (25 * a*b);
	Kl(1, 8) = Kl(4, 13) = -Kl(0, 9);

	Kl(0, 10) = Kl(4, 14) = 27 * b*Dx / (35 * a*a*a) - 6 * D1 / (25 * a*b) + 13 * a*Dy / (35 * b*b*b) - 12 * Dxy / (25 * a*b);
	Kl(2, 8) = Kl(6, 12) = -Kl(0, 10);

	Kl(0, 11) = Kl(3, 8) = Kl(5, 14) = Kl(6, 13) = -13 * b*Dx / (70 * a*a*a) + D1 / (50 * a*b) - 13 * a*Dy / (70 * b*b*b) + Dxy / (25 * a*b);
	Kl(1, 10) = Kl(2, 9) = Kl(4, 15) = Kl(7, 12) = -Kl(0, 11);

	Kl(0, 12) = Kl(4, 8) = 54 * b*Dx / (35 * a*a*a) - 72 * D1 / (25 * a*b) - 156 * a*Dy / (35 * b*b*b) - 144 * Dxy / (25 * a*b);

	Kl(0, 13) = Kl(1, 12) = -13 * b*Dx / (35 * a*a*a) + 6 * D1 / (25 * a*b) + 78 * a*D1 / (35 * b*b*b) + 12 * Dxy / (25 * a*b);
	Kl(4, 9) = Kl(5, 8) = -Kl(0, 13);

	Kl(0, 14) = Kl(4, 10) = 27 * b*Dx / (35 * a*a*a) - 36 * D1 / (25 * a*b) - 22 * a*Dy / (35 * b*b*b) - 12 * Dxy / (25 * a*b);
	Kl(2, 12) = Kl(6, 8) = -Kl(0, 14);

	Kl(0, 15) = Kl(1, 14) = Kl(6, 9) = Kl(7, 8) = -13 * b*Dx / (70 * a*a*a) + 3 * D1 / (25 * a*b) + 11 * a*Dy / (35 * b*b*b) + Dxy / (25 * a*b);
	Kl(2, 13) = Kl(3, 12) = Kl(4, 11) = Kl(5, 10) = -Kl(0, 15);

	Kl(1, 1) = Kl(5, 5) = Kl(9, 9) = Kl(13, 13) = 52 * a*Dy / (35 * b*b*b) + 8 * D1 / (25 * a*b) + 4 * b*Dx / (35 * a*a*a) + 16 * Dxy / (25 * a*b);

	Kl(1, 2) = Kl(9, 10) = 11 * b*Dx / (35 * a*a*a) + 61 * D1 / (50 * a*b) + 11 * a*Dy / (35 * b*b*b) + Dxy / (25 * a*b);
	Kl(5, 6) = Kl(13, 14) = -Kl(1, 2);

	Kl(1, 3) = Kl(5, 7) = 22 * a*Dy / (105 * b*b*b) + 4 * D1 / (25 * a*b) + 2 * b*Dx / (35 * a*a*a) + 4 * Dxy / (75 * a*b);
	Kl(9, 11) = Kl(13, 15) = -Kl(1, 3);

	Kl(1, 5) = Kl(9, 13) = -4 * b*Dx / (35 * a*a*a) - 8 * D1 / (25 * a*b) + 18 * a*Dy / (35 * b*b*b) - 16 * Dxy / (25 * a*b);

	Kl(1, 7) = Kl(3, 5) = 2 * b*Dx / (35 * a*a*a) + 2 * D1 / (75 * a*b) - 13 * a*Dy / (105 * b*b*b) + 4 * Dxy / (75 * a*b);
	Kl(9, 15) = Kl(11, 13) = -Kl(1, 7);

	Kl(1, 9) = Kl(5, 13) = 9 * a*Dy / (35 * b*b*b) + 2 * D1 / (25 * a*b) + 3 * b*Dx / (35 * a*a*a) + 4 * Dxy / (25 * a*b);

	Kl(1, 11) = Kl(5, 15) = -13 * a*Dy / (210 * b*b*b) - D1 / (150 * a*b) - 3 * b*Dx / (70 * a*a*a) - Dxy / (75 * a*b);
	Kl(3, 9) = Kl(7, 13) = -Kl(1, 11);

	Kl(1, 13) = Kl(5, 9) = -3 * b*Dx / (35 * a*a*a) - 6 * D1 / (75 * a*b) + 26 * a*Dy / (35 * b*b*b) - 12 * Dxy / (75 * a*b);

	Kl(1, 15) = Kl(5, 11) = -3 * b*Dx / (70 * a*a*a) - D1 / (25 * a*b) + 11 * a*Dy / (105 * b*b*b) - Dxy / (75 * a*b);
	Kl(3, 13) = Kl(7, 9) = -Kl(1, 15);

	Kl(2, 2) = Kl(6, 6) = Kl(10, 10) = Kl(14, 14) = 4 * a*Dy / (35 * b*b*b) + 8 * D1 / (25 * a*b) + 52 * b*Dx / (35 * a*a*a) + 16 * Dxy / (25 * a*b);

	Kl(2, 3) = Kl(14, 15) = 2 * a*Dy / (35 * b*b*b) + 4 * D1 / (25 * a*b) + 22 * b*Dx / (105 * a*a*a) + 4 * Dxy / (75 * a*b);
	Kl(6, 7) = Kl(10, 11) = -Kl(2, 3);

	Kl(2, 6) = Kl(10, 14) = 26 * b*Dx / (35 * a*a*a) - 2 * D1 / (25 * a*b) - 3 * a*Dy / (35 * b*b*b) - 4 * Dxy / (25 * a*b);

	Kl(2, 7) = Kl(11, 14) = 11 * b*Dx / (105 * a*a*a) - D1 / (25 * a*b) - 3 * a*Dy / (70 * b*b*b) - Dxy / (75 * a*b);
	Kl(3, 6) = Kl(10, 15) = -Kl(2, 7);

	Kl(2, 10) = Kl(6, 14) = 3 * a*Dy / (35 * b*b*b) + 2 * D1 / (25 * a*b) + 9 * b*Dx / (35 * a*a*a) + 12 * Dxy / (75 * a*b);

	Kl(2, 11) = Kl(7, 14) = -3 * a*Dy / (70 * b*b*b) - D1 / (150 * a*b) - 13 * b*Dx / (210 * a*a*a) - Dx / (75 * a*b);
	Kl(3, 10) = Kl(6, 15) = -Kl(2, 11);

	Kl(2, 14) = Kl(6, 10) = 18 * b*Dx / (35 * a*a*a) - 8 * D1 / (25 * a*b) - 4 * a*Dy / (35 * b*b*b) - 16 * Dxy / (25 * a*b);

	Kl(2, 15) = Kl(3, 14) = -13 * b*Dx / (105 * a*a*a) + 2 * D1 / (75 * a*b) + 2 * a*Dy / (35 * b*b*b) + 4 * Dxy / (75 * a*b);
	Kl(6, 11) = Kl(7, 10) = -Kl(2, 15);

	Kl(3, 3) = Kl(7, 7) = Kl(11, 11) = Kl(15, 15) = 4 * b*Dx / (105 * a*a*a) + 8 * D1 / (225 * a*b) + 4 * a*Dy / (105 * b*b*b) + 16 * Dxy / (225 * a*b);

	Kl(3, 7) = Kl(11, 15) = 2 * b*Dx / (105 * a*a*a) - 2 * D1 / (225 * a*b) - a * Dy / (35 * b*b*b) - 4 * Dxy / (225 * a*b);

	Kl(3, 11) = Kl(7, 15) = -b * Dx / (70 * a*a*a) + D1 / (450 * a*b) - a * Dy / (70 * b*b*b) + Dxy / (225 * a*b);

	Kl(3, 15) = Kl(7, 11) = -b * Dx / (35 * a*a*a) - 2 * D1 / (225 * a*b) + 2 * a*Dy / (105 * b*b*b) - 4 * Dxy / (225 * a*b);

	for (int i = 0; i < 16; i++)
		for (int j = 0; j < 16; j++)
			if (i > j) Kl(i, j) = Kl(j, i);

	//std::cout << "Make Kl " << i << "," << j << endl;
	return Kl;
}

VectorXd FEMPL_RECT::get_Rlocal(int i, int j)
{
	double a, b;
	double p = press[i][j];
	a = xmesh[i + 1] - xmesh[i];
	b = ymesh[j + 1] - ymesh[j];
	VectorXd Rl(16);
	Rl << 1, 1.0 / 6.0, 1.0 / 6.0, 1.0 / 36.0, 1, 1.0 / 6.0, -1.0 / 6.0, -1.0 / 36.0, 1, -1.0 / 6.0, -1.0 / 6.0, 1.0 / 36.0, 1, -1.0 / 6.0, 1.0 / 6.0, -1.0 / 36.0;
	Rl = Rl * p*a*b / 4;
	return Rl;
}

MatrixXd FEMPL_RECT::get_Bi(int i, double x, double y, double a, double b)
{

	MatrixXd B(3, 4);
	switch (i)
	{
	case 1:
		B.row(0) << Nij(1, 1, x, y, 2, 0), Nij(1, 2, x, y, 2, 0), Nij(2, 1, x, y, 2, 0), Nij(2, 2, x, y, 2, 0);
		B.row(1) << Nij(1, 1, x, y, 0, 2), Nij(1, 2, x, y, 0, 2), Nij(2, 1, x, y, 0, 2), Nij(2, 2, x, y, 0, 2);
		B.row(2) << Nij(1, 1, x, y, 1, 1), Nij(1, 2, x, y, 1, 1), Nij(2, 1, x, y, 1, 1), Nij(2, 2, x, y, 1, 1);
		break;
	case 2:
		B.row(0) << Nij(3, 1, x, y, 2, 0), Nij(3, 2, x, y, 2, 0), Nij(4, 1, x, y, 2, 0), Nij(4, 2, x, y, 2, 0);
		B.row(1) << Nij(3, 1, x, y, 0, 2), Nij(3, 2, x, y, 0, 2), Nij(4, 1, x, y, 0, 2), Nij(4, 2, x, y, 0, 2);
		B.row(2) << Nij(3, 1, x, y, 1, 1), Nij(3, 2, x, y, 1, 1), Nij(4, 1, x, y, 1, 1), Nij(4, 2, x, y, 1, 1);
		break;
	case 4:
		B.row(0) << Nij(1, 3, x, y, 2, 0), Nij(1, 4, x, y, 2, 0), Nij(2, 3, x, y, 2, 0), Nij(2, 4, x, y, 2, 0);
		B.row(1) << Nij(1, 3, x, y, 0, 2), Nij(1, 4, x, y, 0, 2), Nij(2, 3, x, y, 0, 2), Nij(2, 4, x, y, 0, 2);
		B.row(2) << Nij(1, 3, x, y, 1, 1), Nij(1, 4, x, y, 1, 1), Nij(2, 3, x, y, 1, 1), Nij(2, 4, x, y, 1, 1);
		break;
	case 3:
		B.row(0) << Nij(3, 3, x, y, 2, 0), Nij(3, 4, x, y, 2, 0), Nij(4, 3, x, y, 2, 0), Nij(4, 4, x, y, 2, 0);
		B.row(1) << Nij(3, 3, x, y, 0, 2), Nij(3, 4, x, y, 0, 2), Nij(4, 3, x, y, 0, 2), Nij(4, 4, x, y, 0, 2);
		B.row(2) << Nij(3, 3, x, y, 1, 1), Nij(3, 4, x, y, 1, 1), Nij(4, 3, x, y, 1, 1), Nij(4, 4, x, y, 1, 1);
		break;
	}
	B.row(0) = B.row(0) / a / a;
	B.row(1) = B.row(0) / b / b;
	B.row(2) = 2 * B.row(0) / a / b;
	B = -B;
	return B;
}



void FEMPL_RECT::calc_elem_moments(VectorXd u, int i, int j, double x, double y, VectorXd& M)
{
	double a, b;
	a = xmesh[i + 1] - xmesh[i];
	b = ymesh[j + 1] - ymesh[j];
	double E, nu;
	E = EX[i][j];
	nu = PRXY[i][j];
	Matrix3d D;
	D(0, 0) = E * h*h*h / (12 * (1 - nu * nu));
	D(1, 1) = D(0, 0);
	D(0, 1) = D(1, 0) = D(0, 0) * nu;
	D(2, 2) = D(0, 0) * (1 - nu) / 2;
	M.fill(0);
	MatrixXd B(3, 16);
	for (int i1 = 0; i1 < 4; i1++)
		B.block(0, i1 * 4, 3, 4) = get_Bi(i1, x, y, a, b);
	M = D * B * u;
}

void FEMPL_RECT::calc_elem_result(VectorXd u, int i, int j, int num, VectorXd& M, VectorXd& V)
{
	double x, y;
	switch (num)
	{
	case 1:
		x = 0;
		y = 0;
		break;
	case 2:
		x = 1;
		y = 0;
		break;
	case 3:
		x = 1;
		y = 1;
		break;
	case 4:
		x = 0;
		y = 1;
		break;
	default:
		break;
	}

	double a, b;
	a = xmesh[i + 1] - xmesh[i];
	b = ymesh[j + 1] - ymesh[j];
	double E, nu;
	E = EX[i][j];
	nu = PRXY[i][j];

	MatrixXd D(4,5);
	D(0, 0) = E * h*h*h / (12 * (1 - nu * nu));
	D(1, 1) = D(0, 0);
	D(0, 1) = D(1, 0) = D(0, 0) * nu;
	D(2, 2) = D(0, 0) * (1 - nu) / 2;
	D(3, 3) = D(0, 0);
	D(3, 4) = D(0, 0)*nu;
	MatrixXi Ind(2, 16);
	Ind.row(0) << 1, 1, 2, 2, 3, 3, 4, 4, 3, 3, 4, 4, 1, 1, 2, 2;
	Ind.row(1) << 1, 2, 1, 2, 1, 2, 1, 2, 3, 4, 3, 4, 3, 4, 3, 4;

	VectorXd u_elem(16);
	Vector4i num_n;
	num_n << i * y_size + j, (i + 1)*y_size + j, (i + 1)*y_size + j + 1, i*y_size + j + 1;
	for (int i = 0; i < 4; i++)
	{
		u_elem.segment(i * 4, 4) = u.segment(num_n(i) * 4, 4);
	}

	MatrixXd Bm(5, 16), Bv(3, 16);

	for (int i1 = 0; i1 < 16; i1++)
	{
		Bm(0, i1) = Nij(Ind(0, i1), Ind(1, i1), x, y, 2, 0) / a / a;
		Bm(1, i1) = Nij(Ind(0, i1), Ind(1, i1), x, y, 0, 2) / b / b;
		Bm(2, i1) = 2 * Nij(Ind(0, i1), Ind(1, i1), x, y, 1, 1) / a / b;
		Bm(3, i1) = Nij(Ind(0, i1), Ind(1, i1), x, y, 2, 1) / a / a / b;
		Bm(4, i1) = Nij(Ind(0, i1), Ind(1, i1), x, y, 0, 3) / b / b / b;
		Bv(0, i1) = Nij(Ind(0, i1), Ind(1, i1), x, y, 3, 0) / a / a / a + (2 - nu)*Nij(Ind(0, i1), Ind(1, i1), x, y, 1, 2) / b / b / a;
		Bv(1, i1) = Nij((2 - nu)*Ind(0, i1), Ind(1, i1), x, y, 2, 1) / a / a / b + Nij(Ind(0, i1), Ind(1, i1), x, y, 0, 3) / b / b / b;
		Bv(2, i1) = Nij(Ind(0, i1), Ind(1, i1), x, y, 3, 1) / a / a / b + (2 - nu)*Nij(Ind(0, i1), Ind(1, i1), x, y, 1, 3) / b / b / b;
	}
	Bv = -D(0, 0) * Bv;
	V = Bv * u_elem;

	M = D* Bm * u_elem;
}

void FEMPL_RECT::calc_node_result(VectorXd u, int i, int j, VectorXd& M, VectorXd& V)
{
	M.fill(0);
	V.fill(0);
	VectorXd M_temp(4);
	VectorXd V_temp(3);

	int count = 0;
	if (i > 0)
	{
		if (j > 0)
		{
			calc_elem_result(u, i - 1, j - 1, 3, M_temp, V_temp);
			M += M_temp;
			V += V_temp;
			count++;
		}
		if (j < y_size - 1)
		{
			calc_elem_result(u, i - 1, j, 2, M_temp, V_temp);
			M += M_temp;
			V += V_temp;
			count++;
		}
	}
	if (i < x_size - 1)
	{
		if (j > 0)
		{
			calc_elem_result(u, i, j - 1, 4, M_temp, V_temp);
			M += M_temp;
			V += V_temp;
			count++;
		}
		if (j < y_size - 1)
		{
			calc_elem_result(u, i, j, 1, M_temp, V_temp);
			M += M_temp;
			V += V_temp;
			count++;
		}
	}
	M = M / count;
	V = V / count;
}

void FEMPL_RECT::constructB()
{
	ofstream f("results/B_FEM.txt");
	B.fill(0);

	for (int i = 0; i < y_size; i++)
	{
		double a, b;
		a = xmesh[1] - xmesh[0];
		b = ymesh[i + 1] - ymesh[i];
		B(8 * i, 4 * i) = 1;
		B(8 * i + 3, 4 * i + 3) = 1;
		B(8 * i + 1, 4 * i + 1) = 1.0/a;
		B(8 * i + 2, 4 * i + 2) = 1.0/b;
	}

	VectorXd u(size_k);
	VectorXd M(4);
	VectorXd V(3);

	for (int i = 0; i < y_size; i++)
	{
		for (int j = 0; j < size_k; ++j) {
			u.fill(0);
			u(j) = 1;
			calc_node_result(u, 0, i, M, V);
			B(8 * i + 4, j) = M(0);
			B(8 * i + 5, j) = V(0);
			B(8 * i + 6, j) = M(3);
			B(8 * i + 7, j) = V(2);
		}
	}
	f << B;
}

void FEMPL_RECT::constructK()
{
	MatrixXd Kl(16, 16);
	VectorXd Rl(16);
	for (int i = 0; i < x_size - 1; i++)
		for (int j = 0; j < y_size - 1; j++)
		{
			Kl = get_Klocal1(i, j);
			Rl = get_Rlocal(i, j);

			Vector4i num;
			num << i * y_size + j, (i + 1)*y_size + j, (i + 1)*y_size + j + 1, i*y_size + j + 1;
			for (int num1 = 0; num1 < 4; ++num1)
				for (int num2 = 0; num2 < 4; ++num2)
				{
					K.block(num[num1] * node_dofs, num[num2] * node_dofs, node_dofs, node_dofs) = K.block(num[num1] * node_dofs, num[num2] * node_dofs, node_dofs, node_dofs) + Kl.block(num1 * node_dofs, num2 * node_dofs, node_dofs, node_dofs);
				}
			for (int num1 = 0; num1 < 4; ++num1)
				R.segment(num[num1] * node_dofs, node_dofs) = R.segment(num[num1] * node_dofs, node_dofs) + Rl.segment(num1 * node_dofs, node_dofs);
		}
}

void FEMPL_RECT::add_bc(int num_node, int type)
{
	switch (type)
	{
	case 1:
		for (int i = 0; i < 4; i++)
		{
			K.row(num_node * node_dofs + i).fill(0);
			K.col(num_node * node_dofs + i).fill(0);
			K(num_node * node_dofs + i, num_node * node_dofs + i) = 1;
			R(num_node*node_dofs + i) = 0;
		}
		break;
	case 2:
		for (int i = 3; i < 4; i++)
		{
			K.row(num_node * node_dofs + i).fill(0);
			K.col(num_node * node_dofs + i).fill(0);
			K(num_node * node_dofs + i, num_node * node_dofs + i) = 1;
			R(num_node*node_dofs + i) = 0;
		}
		break;
	default:
		break;
	}
}
