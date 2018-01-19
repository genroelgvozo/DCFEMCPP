#include "DCFEM.h"


double DCFEMPL::Ni(int i, double x, int s1, double h)
{
	switch (s1)
	{
	case 0:
		switch (i)
		{
		case 1:
			return 1 - 3 * x*x + 2 * x*x*x;
		case 2:
			return h * (x - 2 * x*x + x * x*x);
		case 3:
			return h * (3 * x*x - 2 * x*x*x);
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
			return h * (1 - 4 * x + 3 * x*x);
		case 3:
			return h * (6 * x - 6 * x*x);
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
			return h * (6 * x - 4);
		case 3:
			return h * (6 - 12 * x);
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
			return 6 * h;
		case 3:
			return -12 * h;
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

Matrix4d DCFEMPL::get_K0(int i)
{
	double h = xmesh[i + 1] - xmesh[i];
	double E = EX[i];
	double nu = PRXY[i];
	double D = E * dh*dh*dh / (12 * (1 - nu * nu));
	Matrix4d K0, H;
	H.fill(0);
	H(0, 0) = 1;
	H(1, 1) = h;
	H(2, 2) = 1;
	H(3, 3) = h;
	K0 << 6, 3, -6, 3, 3, 2, -3, 1, -6, -3, 6, -3, 3, 1, -3, 2;

	K0 = H * K0*H;
	K0 = 2 * D*K0 / h / h / h;

	return K0;

}

Matrix4d DCFEMPL::get_K4(int i)
{
	double h = xmesh[i + 1] - xmesh[i];
	double E = EX[i];
	double nu = PRXY[i];
	double D = E * dh*dh*dh / (12 * (1 - nu * nu));
	Matrix4d K4, H;
	H.fill(0);
	H(0, 0) = 1;
	H(1, 1) = h;
	H(2, 2) = 1;
	H(3, 3) = h;
	K4 << 156, 22, 54, -13,
		22, 4, 13, -3,
		54, 13, 156, -22,
		-13, -3, -22, 4;
	K4 = H * K4*H;
	K4 = D * h*K4 / 420;
	return K4;
}

void DCFEMPL::constructB()
{
	B1.resize(num_points*num_dofs, size_a);
	B1.fill(0);

	for (int i = 0; i < num_points; i++)
	{
		B1(i*num_dofs, i*num_dofs) = 1;
		B1(i*num_dofs + 1, i*num_dofs + 1) = 1;
		B1(i*num_dofs + 2, i*num_dofs + 4) = 1;
		B1(i*num_dofs + 3, i*num_dofs + 5) = 1;
	}

	B2.resize(size_a, size_a);
	B2.fill(0);
	for (int i = 0; i < num_points; i++)
	{
		B1(i * 2 * num_dofs, i * 2 * num_dofs) = 1;				// displacement
		B1(i * 2 * num_dofs + 1, i * 2 * num_dofs + 4) = 1;		// theta_x = d_1(y) = z1
		B1(i * 2 * num_dofs + 2, i * 2 * num_dofs + 1) = 1;		// theta_y = d_2(y) = y2
		B1(i * 2 * num_dofs + 3, i * 2 * num_dofs + 5) = 1;		// wxy = z2
		double h = xmesh[i + 1] - xmesh[i];
		double E = EX[i];
		double nu = PRXY[i];
		double D = E * dh*dh*dh / (12 * (1 - nu * nu));
	}
}

void DCFEMPL::construct()
{
	int bsize = size_a / 4;
	for (int i = 0; i < 3 * bsize; ++i)
		mat_a(i, i + bsize) = 1;
	MatrixXd K4(bsize, bsize), K0(bsize, bsize), K2(bsize, bsize);

	Matrix4d K4l, K0l, K2l;

	for (int i = 0; i < xmesh.size() - 1; ++i)
	{
		K4l = get_K4(i);
		K0l = get_K0(i);
		K2l = get_K2(i);

		K4.block(i * 2, i * 2, 4, 4) += K4l;
		K0.block(i * 2, i * 2, 4, 4) += K0l;
		K2.block(i * 2, i * 2, 4, 4) += K2l;

	}
	ofstream f("logs/K4.txt");
	f << K4 << endl;
	f.close();
	K4 = K4.inverse();
	f.open("logs/K4_inv.txt");
	f << K4 << endl;
	f.close();
	mat_a.block(3 * bsize, 0, bsize, bsize) = K4 * K0;
	mat_a.block(3 * bsize, 2 * bsize, bsize, bsize) = K4 * K2;
	for (int i = 0; i < loads.size(); ++i)
	{
		VectorXd R1(size_a);
		R1.fill(0);
		if (loads[i].type == 1)
			R1[size_a / 2 + (loads[i].num_node - 1)*num_dofs + 1] = loads[i].value;
		R1 = K4 * R1;
		R.push_back({ loads[i].x2,R1 });
	}

}

Matrix4d DCFEMPL::get_K2(int i)
{
	double h = xmesh[i + 1] - xmesh[i];
	double E = EX[i];
	double nu = PRXY[i];
	double D = E * dh*dh*dh / (12 * (1 - nu * nu));
	Matrix4d K2_3, K2_2, H;

	H.fill(0);
	H(0, 0) = 1;
	H(1, 1) = h;
	H(2, 2) = 1;
	H(3, 3) = h;
	K2_3 << -36, -33, 36, -3,
		-3, -4, 3, 1,
		36, 3, -36, 33,
		-3, 1, 3, -4;
	K2_3 = H * K2_3*H;
	K2_3 = -D * nu*K2_3 / 30 / h;
	K2_3 = K2_3.eval() + K2_3.transpose().eval();
	K2_2 << 36, 3, -36, 3,
		3, 4, -3, -1,
		-36, -3, 36, -3,
		3, -1, -3, 4;
	K2_2 = H * K2_2*H;
	K2_2 = D * (1 - nu)*K2_2 / 15 / h;

	return K2_2 + K2_3;
}
