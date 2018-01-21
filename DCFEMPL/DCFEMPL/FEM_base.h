#pragma once
#include <memory>
#include <vector>
#include "eigen\Eigen\Dense"

using namespace std;
using namespace Eigen;

class FEM_base {


protected:
	unsigned int node_dofs;
	int size_k;

public:
	MatrixXd K;
	VectorXd R;

	FEM_base(int num_nodes, int _node_dofs) : node_dofs(_node_dofs), size_k(num_nodes*_node_dofs), K(MatrixXd(size_k,size_k)), R(VectorXd(size_k)) {
		K.fill(0);
		R.fill(0);
	}


	int get_size() {
		return size_k;
	}

	VectorXd solve() {
		return K.colPivHouseholderQr().solve(R);
	}
};