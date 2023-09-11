#include "FMMnDTree.hpp"

class AFMMnD {
public:
	int ndim;
	int n; // n=pow(N,1/ndim)
	int nParticlesInLeafAlong1D; // assuming the particles are located at tensor product chebyNodes
	double L;
	int TOL_POW;
	int Qchoice;
	int nLevels;
	std::vector<AptsnD> particles;
	userkernel* mykernel;
	FMM3DTree<userkernel>* A;

	AFMMnD (int ndim, int N, int nParticlesInLeafAlong1D, double L, int TOL_POW, int Qchoice, int nLevels, std::vector<AptsnD> points) {
		this->ndim = ndim;
		this->n = n;
		this->nParticlesInLeafAlong1D = nParticlesInLeafAlong1D;
		this->L = L;
		this->TOL_POW = TOL_POW;
		this->Qchoice = Qchoice;
		this->nLevels = nLevels;

		mykernel		=	new userkernel(particles, Qchoice);
		// std::cout << "AFMM N: " << N << std::endl;
		// std::cout << "nLevels: " << nLevels << std::endl;
		A	=	new FMM3DTree<userkernel>(ndim, mykernel, N, nLevels, nParticlesInLeafAlong1D, L, TOL_POW);
		// std::cout << "points.size(): " << points.size() << std::endl;
		A->setNodes(points);
		A->createTree();
		A->assign_Center_Location();
		A->assign_Tree_Interactions();
		A->assignChargeLocations();
		A->assignNonLeafChargeLocations();
	}

	void assembleAFMMnD() {
		A->getNodes();
		// std::cout << "getNodes done" << std::endl;
		A->assemble_M2L();
	}

	Eigen::VectorXd matVecProduct(Eigen::VectorXd x) {
		Eigen::VectorXd Ax;
		A->assignCharges(x);
		// std::cout << "x.size(): " << x.size() << std::endl;
		// std::cout << "assignCharges done" << std::endl;
		A->evaluate_M2M();
		// std::cout << "evaluate_M2M done" << std::endl;
		A->evaluate_M2L();
		// std::cout << "evaluate_M2L done" << std::endl;
		A->evaluate_L2L();
		// std::cout << "evaluate_L2L done" << std::endl;
		A->evaluate_NearField();
		// std::cout << "evaluate_NearField done" << std::endl;
		A->collectPotential(Ax);
		// std::cout << "Ax.size(): " << Ax.size() << std::endl;
		// std::cout << "collectPotential done" << std::endl;
		A->reorder(Ax);
		return Ax;
	}

	~AFMMnD() {
		delete A;
		delete mykernel;
	};
};

std::vector<AptsnD> give_random_Nodes(int N, int ndim) {
	std::vector<AptsnD> points;
	for (size_t i = 0; i < N; i++) {
		AptsnD temp1;
		for (size_t j = 0; j < ndim; j++) {
			temp1.push_back((double(rand())/double(RAND_MAX)-0.5)*2);
		}
		points.push_back(temp1);
	}
	return points;
}
