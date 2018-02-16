#include "DCFEM.h"
#include <fstream>

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
			return (3 * x*x - 2 * x*x*x);
		case 4:
			return h*(x * x*x - x * x);
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
			return (6 * x - 6 * x*x);
		case 4:
			return h*(3 * x*x - 2 * x);
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
			return (6 - 12 * x);
		case 4:
			return h*(6 * x - 2);
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
			return -12;
		case 4:
			return 6 * h;
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
	Matrix4d K0, H, K0_2;
	H.fill(0);
	H(0, 0) = 1;
	H(1, 1) = h;
	H(2, 2) = 1;
	H(3, 3) = h;
	K0 << 6, 3, -6, 3, 
		3, 2, -3, 1, 
		-6, -3, 6, -3, 
		3, 1, -3, 2;

	K0 = H * K0*H;
	K0 = -2 * D*K0 / h / h / h;

	K0_2 << 156, 22, 54, -13,
		22, 4, 13, -3,
		54, 13, 156, -22,
		-13, -3, -22, 4;
	K0_2 = H* K0_2 * H;
	K0_2 = D / 2 * K0_2 * h / 420;
	//K0 = K0 +  K0_2;
	return K0;

}

Matrix4d DCFEMPL::get_K4(int i)
{
	double h = xmesh[i + 1] - xmesh[i];
	double E = EX[i];
	double nu = PRXY[i];
	double D = E * dh*dh*dh / (12 * (1 - nu * nu));
	Matrix4d K4, H;
	K4.fill(0);
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
	ofstream f("results/B1_DCFEM.txt");
	B1.resize(num_points * 4, num_points * 8);
	B1.fill(0);

	for (int i = 0; i < num_points; i++)
	{
		B1(i * 4, 2 * i) = 1;
		B1(i * 4 + 1, 2 * i + 1) = 1;
		B1(i * 4 + 2, 2 * num_points + 2 * i) = 1;
		B1(i * 4 + 3, 2 * num_points + 2 * i + 1) = 1;
	}

	f << B1;
	f.close();
	f.open("results/B2_DCFEM.txt");

	B2.resize(size_a, size_a);
	B2.fill(0);
	for (int i = 0; i < num_points; i++)
	{
		B2(i * 4 + 0, 2 * i) = 1;						// displacement
		B2(i * 4 + 1, 2 * i + 1) = 1;						// theta_x = d_1(y) = z1
		B2(i * 4 + 2, 2 * num_points + 2 * i) = 1;	// theta_y = d_2(y) = y2
		B2(i * 4 + 3, 2 * num_points + 2 * i + 1) = 1;		// wxy = z2
		int factor = 0;
		double h, E, nu, D;
		/*
		if (i > 0) {
			factor++;
			h = xmesh[i] - xmesh[i - 1];
			E = EX[i - 1];
			nu = PRXY[i - 1];
			D = E * dh*dh*dh / (12 * (1 - nu * nu));

			// M2
			B2(i * 8 + 4, 2 * i - 2) += -D * nu*Ni(1, 1, 2, h) / h / h;
			B2(i * 8 + 4, 2 * i - 1) += -D * nu*Ni(2, 1, 2, h) / h / h;
			B2(i * 8 + 4, 2 * i) += -D * nu*Ni(3, 1, 2, h) / h / h;
			B2(i * 8 + 4, 2 * i + 1) += -D * nu*Ni(4, 1, 2, h) / h / h;
			B2(i * 8 + 4, 4 * num_points + 2 * i) += -D;

			// V2
			B2(i * 8 + 5, 2 * num_points + 2 * i - 2) += -D * (2 - nu)*Ni(1, 1, 2, h) / h / h;
			B2(i * 8 + 5, 2 * num_points + 2 * i - 1) += -D * (2 - nu)*Ni(2, 1, 2, h) / h / h;
			B2(i * 8 + 5, 2 * num_points + 2 * i) += -D * (2 - nu)*Ni(3, 1, 2, h) / h / h;
			B2(i * 8 + 5, 2 * num_points + 2 * i + 1) += -D * (2 - nu)*Ni(4, 1, 2, h) / h / h;
			B2(i * 8 + 5, 6 * num_points + 2 * i) += -D;

			// d1_M2
			B2(i * 8 + 6, 2 * i - 2) += -D * nu*Ni(1, 1, 3, h) / h / h / h;
			B2(i * 8 + 6, 2 * i - 1) += -D * nu*Ni(2, 1, 3, h) / h / h / h;
			B2(i * 8 + 6, 2 * i) += -D * nu*Ni(3, 1, 3, h) / h / h / h;
			B2(i * 8 + 6, 2 * i + 1) += -D * nu*Ni(4, 1, 3, h) / h / h / h;
			B2(i * 8 + 6, 4 * num_points + 2 * i + 1) += -D;

			// d1_V2
			B2(i * 8 + 7, 2 * num_points + 2 * i - 2) += -D * (2 - nu)*Ni(1, 1, 3, h) / h / h / h;
			B2(i * 8 + 7, 2 * num_points + 2 * i - 1) += -D * (2 - nu)*Ni(2, 1, 3, h) / h / h / h;
			B2(i * 8 + 7, 2 * num_points + 2 * i) += -D * (2 - nu)*Ni(3, 1, 3, h) / h / h / h;
			B2(i * 8 + 7, 2 * num_points + 2 * i + 1) += -D * (2 - nu)*Ni(4, 1, 3, h) / h / h / h;
			B2(i * 8 + 7, 6 * num_points + 2 * i + 1) += -D;

		}
		if (i < num_points)
		{
			factor++;
			h = xmesh[i + 1] - xmesh[i];
			E = EX[i];
			nu = PRXY[i];
			D = E * dh*dh*dh / (12 * (1 - nu * nu));

			// M2
			B2(i * 8 + 4, 2 * i) += -D * nu*Ni(1, 0, 2, h) / h / h;
			B2(i * 8 + 4, 2 * i + 1) += -D * nu*Ni(2, 0, 2, h) / h / h;
			B2(i * 8 + 4, 2 * i + 2) += -D * nu*Ni(3, 0, 2, h) / h / h;
			B2(i * 8 + 4, 2 * i + 3) += -D * nu*Ni(4, 0, 2, h) / h / h;
			B2(i * 8 + 4, 4 * num_points + 2 * i) += -D;

			// V2
			B2(i * 8 + 5, 2 * num_points + 2 * i) += -D * (2 - nu)*Ni(1, 0, 2, h) / h / h;
			B2(i * 8 + 5, 2 * num_points + 2 * i + 1) += -D * (2 - nu)*Ni(2, 0, 2, h) / h / h;
			B2(i * 8 + 5, 2 * num_points + 2 * i + 2) += -D * (2 - nu)*Ni(3, 0, 2, h) / h / h;
			B2(i * 8 + 5, 2 * num_points + 2 * i + 3) += -D * (2 - nu)*Ni(4, 0, 2, h) / h / h;
			B2(i * 8 + 5, 6 * num_points + 2 * i) += -D;

			// d1_M2
			B2(i * 8 + 6, 2 * i) += -D * nu*Ni(1, 0, 3, h) / h / h / h;
			B2(i * 8 + 6, 2 * i + 1) += -D * nu*Ni(2, 0, 3, h) / h / h / h;
			B2(i * 8 + 6, 2 * i + 2) += -D * nu*Ni(3, 0, 3, h) / h / h / h;
			B2(i * 8 + 6, 2 * i + 3) += -D * nu*Ni(4, 0, 3, h) / h / h / h;
			B2(i * 8 + 6, 4 * num_points + 2 * i + 1) += -D;

			// d1_V2
			B2(i * 8 + 7, 2 * num_points + 2 * i) += -D * (2 - nu)*Ni(1, 0, 3, h) / h / h / h;
			B2(i * 8 + 7, 2 * num_points + 2 * i + 1) += -D * (2 - nu)*Ni(2, 0, 3, h) / h / h / h;
			B2(i * 8 + 7, 2 * num_points + 2 * i + 2) += -D * (2 - nu)*Ni(3, 0, 3, h) / h / h / h;
			B2(i * 8 + 7, 2 * num_points + 2 * i + 3) += -D * (2 - nu)*Ni(4, 0, 3, h) / h / h / h;
			B2(i * 8 + 7, 6 * num_points + 2 * i + 1) += -D;
		}
		*/
		//B2.block(i * 8 + 4, 0, 4, size_a) /= factor;
	}


	f << B2;
}

void DCFEMPL::construct()
{
	int bsize = size_a / 4;
	for (int i = 0; i < 3 * bsize; ++i)
		mat_a(i, i + bsize) = 1;
	MatrixXd K4(bsize, bsize), K0(bsize, bsize), K2(bsize, bsize);
	K4.fill(0);
	K0.fill(0);
	K2.fill(0);
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
	ofstream f("results/K4.txt");
	f.precision(3);
	f << fixed;
	f << K4 << endl;
	f.close();
	K4 = K4.inverse();
	f.open("results/K4_inv.txt");
	f << K4 << endl;
	f.close();
	f.open("results/K2.txt");
	f << K2 << endl;
	f.close();
	mat_a.block(3 * bsize, 0, bsize, bsize) = K4 * K0;
	mat_a.block(3 * bsize, 2 * bsize, bsize, bsize) = K4 * K2;
	f.open("results/mat_a.txt");
	f << mat_a << endl;
	f.close();

	for (int i = 0; i < loads.size(); ++i)
	{
		VectorXd R1(size_a);
		R1.fill(0);
		if (loads[i].type == 1)
			R1[3 * bsize + (loads[i].num_node - 1) * 2] = loads[i].value;
		R1.segment(3 * bsize, bsize) = -K4 * R1.segment(3 * bsize, bsize);
		R.push_back({ loads[i].x2,R1 });
	}

}

Matrix4d DCFEMPL::get_K2(int i)
{
	double h = xmesh[i + 1] - xmesh[i];
	double E = EX[i];
	double nu = PRXY[i];
	double D = E * dh*dh*dh / (12 * (1 - nu * nu));
	Matrix4d K2, K2_3, K2_2, H;

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
	K2 = K2_3;
	K2 += K2_3.transpose().eval();
	K2_2 << 36, 3, -36, 3,
		     3, 4, -3, -1,
		  -36, -3, 36, -3,
		    3, -1, -3, 4;
	K2_2 = H * K2_2*H;
	K2_2 = -D * (1 - nu)*K2_2 / 15 / h;
	K2 = K2 + K2_2;
	return K2;
}
