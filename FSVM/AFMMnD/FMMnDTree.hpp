#ifndef _FMMnDTreeRAMeff2_HPP__
#define _FMMnDTreeRAMeff2_HPP__

#include "kernel.hpp"
#include "ACA.hpp"
#include <boost/dynamic_bitset.hpp>

class FMMnDBox {
public:
	int boxNumber;
	int parentNumber;
	std::vector<int> childrenNumbers;
	std::vector<int> interactionList;
	std::vector<int> neighborNumbers;
	FMMnDBox (int ndim) {
		boxNumber		=	-1;
		parentNumber	=	-1;
		childrenNumbers.resize(pow(2,ndim));
	}

	AptsnD center;
	Eigen::VectorXd charges, potential;
  std::vector<int> chargeLocations;
	Eigen::VectorXd outgoing_charges;//equivalent densities {f_{k}^{B,o}}
	Eigen::VectorXd incoming_charges;//equivalent densities {f_{k}^{B,i}}
	Eigen::VectorXd incoming_potential;//check potentials {u_{k}^{B,i}}
	std::vector<int> incoming_chargePoints;//equivalent points {y_{k}^{B,i}}
	std::vector<int> incoming_checkPoints;//check points {x_{k}^{B,i}}
	std::map<int, Eigen::MatrixXd> M2L;
	Eigen::MatrixXd L2P;					//	Transfer from multipoles of 4 children to multipoles of parent.
};

template <typename kerneltype>
class FMM3DTree {
public:
	kerneltype* K;
	int nLevels;			//	Number of levels in the tree.
	int N;					//	Number of particles.
	// int ndimRootN;					//	Number of particles.
	double L;				//	Semi-length of the simulation box.
	double smallestBoxSize;	//	This is L/2.0^(nLevels).

	std::vector<int> nBoxesPerLevel;			//	Number of boxes at each level in the tree.
	std::vector<double> boxRadius;				//	Box radius at each level in the tree assuming the box at the root is [-1,1]^2
	std::vector<std::vector<FMMnDBox> > tree;	//	The tree storing all the information.

	double ACA_epsilon;
	int nParticlesInLeafAlong1D;
  int nParticlesInLeaf;
  std::vector<double> Nodes1D;
	std::vector<AptsnD> Nodes;
  std::vector<AptsnD> gridPoints; //all particles in domain
	int TOL_POW;
	double assTime, matVecTime;
	double pow2d;
	int ndim;
	double tolerance;
// public:
FMM3DTree(int ndim, kerneltype* K, int N, int nLevels, int nParticlesInLeafAlong1D, double L, int TOL_POW) {
	tolerance = pow(10,-15);
	this->ndim					=	ndim;
	this->pow2d = pow(2,ndim);
	this->K					=	K;
	// this->ndimRootN	=	ndimRootN;
	this->nLevels		=	nLevels;
		this->L					=	L;
    this->nParticlesInLeafAlong1D = nParticlesInLeafAlong1D;
    this->nParticlesInLeaf = pow(nParticlesInLeafAlong1D, pow2d);
		this->TOL_POW = TOL_POW;
    nBoxesPerLevel.push_back(1);
		boxRadius.push_back(L);
		for (int k=1; k<=nLevels; ++k) {
			nBoxesPerLevel.push_back(pow2d*nBoxesPerLevel[k-1]);
			boxRadius.push_back(0.5*boxRadius[k-1]);
		}
		this->smallestBoxSize	=	boxRadius[nLevels];
		K->a					=	smallestBoxSize;
		this->N					=	N;
		// std::cout << "N: " << N << std::endl;
		// std::cout << "ndimRootN: " << ndimRootN << std::endl;
		// std::cout << "ndim: " << ndim << std::endl;
		this->assTime = 0.0;
		this->matVecTime = 0.0;

	}

	void setNodes(std::vector<AptsnD> &points) {
		for (size_t i = 0; i < N; i++) {
			K->particles.push_back(points[i]);
		}
	}

  void set_random_Nodes() {
		for (size_t i = 0; i < N; i++) {
			AptsnD temp1;
			for (size_t j = 0; j < ndim; j++) {
				temp1.push_back((double(rand())/double(RAND_MAX)-0.5)*2);
			}
			K->particles.push_back(temp1);
		}
	}

	void createTree() {
		//	First create root and add to tree
		FMMnDBox root(ndim);
		root.boxNumber		=	0;
		root.parentNumber	=	-1;
		// #pragma omp parallel for
		for (int l=0; l<int(pow2d); ++l) {
			root.childrenNumbers[l]=l;
		}
		std::vector<FMMnDBox> rootLevel;
		rootLevel.push_back(root);
		tree.push_back(rootLevel);

		for (int j=1; j<=nLevels; ++j) {
			std::vector<FMMnDBox> level;
			for (int k=0; k<nBoxesPerLevel[j]; ++k) {
				FMMnDBox box(ndim);
				box.boxNumber		=	k;
				box.parentNumber	=	k/pow2d;
				for (int l=0; l<pow2d; ++l) {
					box.childrenNumbers[l]	=	pow2d*k+l;
				}
				level.push_back(box);
			}
			tree.push_back(level);
		}
	}

	bool admCondition(int j, int k, int l) {
		for (size_t i = 0; i < ndim; i++) {
			if (abs(tree[j][k].center[i]-tree[j][l].center[i]) > 4*boxRadius[j] - tolerance) {
				return true; // IL
			}
		}
		return false; // Neighbor
 	}

 	//	Assigns the interactions for the children of a box
	void assign_Box_Interactions(int j, int k) {
		//children of parents neighbors
		int pj = j-1;
		int pk = k/pow2d;
		for (size_t i = 0; i < tree[pj][pk].neighborNumbers.size(); i++) {
			int pn = tree[pj][pk].neighborNumbers[i];
			for (size_t i = 0; i < pow2d; i++) {
				int cpn = pn*pow2d + i;
				if (admCondition(j,k,cpn)) {
					tree[j][k].interactionList.push_back(cpn);
				}
				else {
					tree[j][k].neighborNumbers.push_back(cpn);
				}
			}
		}
		// siblings
		for (size_t i = 0; i < pow2d; i++) {
			if (pk*pow2d+i != k) { //negibors doesnt include self
				tree[j][k].neighborNumbers.push_back(pk*pow2d+i);
			}
		}
	}

	//	Assigns the interactions for the children all boxes at a given level
	void assign_Level_Interactions(int j) {
		// #pragma omp parallel for
		for (int k=0; k<nBoxesPerLevel[j]; ++k) {
			assign_Box_Interactions(j,k);
		}
	}

	//	Assigns the interactions for the children all boxes in the tree
	void assign_Tree_Interactions() {
		for (int j=1; j<=nLevels; ++j) {
			assign_Level_Interactions(j);
		}
	}

	void printIL() {
		for (size_t j = 0; j <= nLevels; j++) {
			for (size_t k = 0; k < nBoxesPerLevel[j]; k++) {
				std::cout << "j: " << j << "	k: " << k << std::endl;
				for (size_t i = 0; i < tree[j][k].interactionList.size(); i++) {
					std::cout << tree[j][k].interactionList[i] << ",";
				}
				 std::cout << std::endl;
			}
		}
	}

	void printN() {
		for (size_t j = 0; j <= nLevels; j++) {
			for (size_t k = 0; k < nBoxesPerLevel[j]; k++) {
				std::cout << "j: " << j << "	k: " << k << std::endl;
				for (size_t i = 0; i < tree[j][k].neighborNumbers.size(); i++) {
					std::cout << tree[j][k].neighborNumbers[i] << ",";
				}
				 std::cout << std::endl;
			}
		}
	}

	void assign_Center_Location() {
		int J;
		for (size_t i = 0; i < ndim; i++) {
			tree[0][0].center.push_back(0.0);
		}
		for (int j=0; j<nLevels; ++j) {
			J	=	j+1;
			double shift	=	0.5*boxRadius[j];
			// #pragma omp parallel for
			for (int k=0; k<nBoxesPerLevel[j]; ++k) {
				for (int c = 0; c < pow2d; c++) {
					boost::dynamic_bitset<> b(ndim, c);
  				for (size_t i = 0; i < ndim; i++) {// w,z,y,x
						if (!b[i]) {
							tree[J][pow2d*k+c].center.push_back(tree[j][k].center[i]-shift);
						}
						else {
							tree[J][pow2d*k+c].center.push_back(tree[j][k].center[i]+shift);
						}
					}
				}
			}
		}
	}

	void print_centers() {
		for (size_t j = 0; j <= nLevels; j++) {
			for (size_t k = 0; k < nBoxesPerLevel[j]; k++) {
				std::cout << "j: " << j << "	k: " << k << "	(";
				for (size_t i = 0; i < ndim; i++) {
					std::cout << tree[j][k].center[i] << ",";
				}
				std::cout << ")" << std::endl;
			}
		}
	}

	void assignChargeLocations() {
		for (size_t i = 0; i < N; i++) {
			tree[0][0].chargeLocations.push_back(i);
		}
		for (size_t j = 0; j < nLevels; j++) { //assign particles to its children
			for (size_t k = 0; k < nBoxesPerLevel[j]; k++) {
				int J = j+1;
				int Kp = pow2d*k;
				for (size_t i = 0; i < tree[j][k].chargeLocations.size(); i++) {
					int index = tree[j][k].chargeLocations[i];
					char b[ndim];
					for (size_t i = 0; i < ndim; i++) {
						if (K->particles[index][i] <= tree[j][k].center[i]) { //children 0,1,2,3
							b[ndim-1-i] = '0';
						}
						else {
							b[ndim-1-i] = '1';
						}
					}
					int c = std::stoi(b, 0, 2);
					tree[J][Kp+c].chargeLocations.push_back(index);
				}
			}
		}
	}

	void printcharges() {
		for (size_t i = 0; i < N; i++) {
			std::cout << "i: " << i << "	";
			for (size_t h = 0; h < ndim; h++) {
				std::cout << K->particles[i][h] << ", ";
			}
			std::cout << std::endl;
		}
	}

	void printchargeLocations() {
		for (size_t j = 0; j <= nLevels; j++) {
			for (size_t k = 0; k < nBoxesPerLevel[j]; k++) {
				std::cout << "j: " << j << "	k: " << k << "	";
				for (size_t i = 0; i < tree[j][k].chargeLocations.size(); i++) {
					std::cout << tree[j][k].chargeLocations[i] << ",";
				}
				std::cout << std::endl;
			}
		}
	}

	void assignCharges(Eigen::VectorXd &charges) {
		int start = 0;
		for (size_t j = 1; j <= nLevels; j++) {
			for (size_t k = 0; k < nBoxesPerLevel[j]; k++) {
				tree[j][k].charges = Eigen::VectorXd::Zero(tree[j][k].chargeLocations.size());
				for (size_t i = 0; i < tree[j][k].chargeLocations.size(); i++) {
					tree[j][k].charges[i] = charges[tree[j][k].chargeLocations[i]];
				}
			}
		}
	}

	void assignNonLeafChargeLocations() {
		for (int j = nLevels-1; j >= 1; j--) {
			for (size_t k = 0; k < nBoxesPerLevel[j]; k++) {
				tree[j][k].chargeLocations.clear();
				for (size_t c = 0; c < pow2d; c++) {
					tree[j][k].chargeLocations.insert(tree[j][k].chargeLocations.end(), tree[j+1][pow2d*k+c].chargeLocations.begin(), tree[j+1][pow2d*k+c].chargeLocations.end());
				}
			}
		}
	}

	void getNodes() {
		// std::cout << "nLevels: " << nLevels << std::endl;
		for (int j=nLevels; j>=2; j--) {
			// std::cout << "j: " << j << std::endl;
			getNodes_incoming_level(j);
		}
	}

	void getParticlesFromChildren_incoming(int j, int k, std::vector<int>& searchNodes) {
		if (j==nLevels) {
			searchNodes.insert(searchNodes.end(), tree[j][k].chargeLocations.begin(), tree[j][k].chargeLocations.end());
		}
		else {
			int J = j+1;
			for (int c = 0; c < pow2d; c++) {
				searchNodes.insert(searchNodes.end(), tree[J][pow2d*k+c].incoming_checkPoints.begin(), tree[J][pow2d*k+c].incoming_checkPoints.end());
			}
		}
	}

	void getNodes_incoming_box(int j, int k, int& ComputedRank) {
		std::vector<int> boxA_Nodes;
		getParticlesFromChildren_incoming(j, k, boxA_Nodes);
		std::vector<int> IL_Nodes;//indices
		for(int in=0; in<tree[j][k].interactionList.size(); ++in) {
			int kIL = tree[j][k].interactionList[in];
			std::vector<int> chargeLocations;
			getParticlesFromChildren_incoming(j, kIL, chargeLocations);
			IL_Nodes.insert(IL_Nodes.end(), chargeLocations.begin(), chargeLocations.end());
		}
		std::vector<int> row_bases, col_bases;
		Eigen::MatrixXd Ac, Ar, L, R;
		LowRank* LR		=	new LowRank(K, TOL_POW, boxA_Nodes, IL_Nodes);
		if (boxA_Nodes.size()>0 && IL_Nodes.size()>0) {
			LR->ACA_only_nodes(row_bases, col_bases, ComputedRank, Ac, Ar, L, R);
		}
		if(ComputedRank > 0) {
			for (int r = 0; r < row_bases.size(); r++) {
				tree[j][k].incoming_checkPoints.push_back(boxA_Nodes[row_bases[r]]);
			}
			for (int c = 0; c < col_bases.size(); c++) {
				tree[j][k].incoming_chargePoints.push_back(IL_Nodes[col_bases[c]]);
			}
			Eigen::MatrixXd temp = R.triangularView<Eigen::Upper>().solve<Eigen::OnTheRight>(Ac);
			tree[j][k].L2P = L.triangularView<Eigen::Lower>().solve<Eigen::OnTheRight>(temp);
		}
	}

	void getNodes_incoming_level(int j) {
		int ComputedRank;
		for (int k=0; k<nBoxesPerLevel[j]; ++k) {
			// std::cout << "j: " << j << "	k: " << k << std::endl;
			getNodes_incoming_box(j, k, ComputedRank);
		}
	}

	void assemble_M2L() {
		// #pragma omp parallel for
		for (size_t j = 2; j <= nLevels; j++) {
			// #pragma omp parallel for
			for (size_t k = 0; k < nBoxesPerLevel[j]; k++) {
				// #pragma omp parallel for
				for(int in=0; in<tree[j][k].interactionList.size(); ++in) {
					int kIL = tree[j][k].interactionList[in];
					if (tree[j][k].M2L[kIL].size() == 0)
						tree[j][k].M2L[kIL] = K->getMatrix(tree[j][k].incoming_checkPoints, tree[j][kIL].incoming_checkPoints);
					if (tree[j][kIL].M2L[k].size() == 0)
						tree[j][kIL].M2L[k] = tree[j][k].M2L[kIL].transpose();
				}
			}
		}
	}

	void evaluate_M2M() {
		for (size_t j = nLevels; j >= 2; j--) {
			// #pragma omp parallel for
			for (size_t k = 0; k < nBoxesPerLevel[j]; k++) {
				Eigen::VectorXd source_densities;
				if (j==nLevels) {
					source_densities = tree[j][k].charges;
				}
				else {
					int J = j+1;
					int Veclength = 0;
					for (int child = 0; child < pow2d; child++) {
						Veclength += tree[J][pow2d*k+child].outgoing_charges.size();
					}
					source_densities = Eigen::VectorXd::Zero(Veclength);// = tree[j][k].multipoles//source densities
					int start = 0;
					for (int child = 0; child < pow2d; child++) {
						int NumElem = tree[J][pow2d*k+child].outgoing_charges.size();
						source_densities.segment(start, NumElem) = tree[J][pow2d*k+child].outgoing_charges;
						start += NumElem;
					}
				}
				if (tree[j][k].L2P.size()>0) {
					tree[j][k].outgoing_charges = tree[j][k].L2P.transpose()*source_densities;//f^{B,o} //solve system: A\tree[j][k].outgoing_potential
				}
			}
		}
	}

	void evaluate_M2L() {
		// #pragma omp parallel for
		for (int j=2; j<=nLevels; ++j) {
			// #pragma omp parallel for
			for (int k=0; k<nBoxesPerLevel[j]; ++k) {//BoxA
				tree[j][k].incoming_potential	=	Eigen::VectorXd::Zero(tree[j][k].incoming_checkPoints.size());
				// #pragma omp parallel for
				for(int in=0; in<tree[j][k].interactionList.size(); ++in) {
					int kIL = tree[j][k].interactionList[in];
					tree[j][k].incoming_potential += tree[j][k].M2L[kIL]*tree[j][kIL].outgoing_charges;
				}
			}
		}
	}

	void evaluate_L2L() {
		for (size_t j = 2; j <= nLevels; j++) {
			// #pragma omp parallel for
			for (size_t k = 0; k < nBoxesPerLevel[j]; k++) {
				if (j != nLevels) {
					Eigen::VectorXd temp = tree[j][k].L2P*tree[j][k].incoming_potential;
					int start = 0;
					for (size_t c = 0; c < pow2d; c++) {
						tree[j+1][pow2d*k+c].incoming_potential += temp.segment(start, tree[j+1][pow2d*k+c].incoming_checkPoints.size());
						start += tree[j+1][pow2d*k+c].incoming_checkPoints.size();
					}
				}
				else {
					tree[j][k].potential = tree[j][k].L2P*tree[j][k].incoming_potential;
				}
			}
		}
	}

	void assignLeafCharges(Eigen::VectorXd &charges) {
		int start = 0;
		for (size_t k = 0; k < nBoxesPerLevel[nLevels]; k++) {
			tree[nLevels][k].charges	=	charges.segment(start, nParticlesInLeaf);
			start += nParticlesInLeaf;
		}
	}

	void assignNonLeafCharges() {
		for (int j = nLevels-1; j >= 2; j--) {
			for (size_t k = 0; k < nBoxesPerLevel[j]; k++) {
				int sum = 0;
				for (size_t c = 0; c < pow2d; c++) {
					sum += tree[j+1][pow2d*k+c].charges.size();
				}
				tree[j][k].charges	=	Eigen::VectorXd(sum);
				int start = 0;
				for (size_t c = 0; c < pow2d; c++) {
					tree[j][k].charges.segment(start, tree[j+1][pow2d*k+c].charges.size()) = tree[j+1][pow2d*k+c].charges;
					start += tree[j+1][pow2d*k+c].charges.size();
				}
			}
		}
	}

	int getAvgRank() {
		int sum = 0;
		int NumBoxes = 0;
		for (size_t j = 2; j <= nLevels; j++) {
			for (size_t k = 0; k < nBoxesPerLevel[j]; k++) {
				NumBoxes += 1;
				sum += tree[j][k].incoming_checkPoints.size();
			}
		}
		return sum/NumBoxes;
	}

	void evaluate_NearField() { // evaluating at chargeLocations
		// #pragma omp parallel for
		for (size_t k = 0; k < nBoxesPerLevel[nLevels]; k++) {
			// #pragma omp parallel for
			for (size_t n = 0; n < tree[nLevels][k].neighborNumbers.size(); n++) {
				int nn = tree[nLevels][k].neighborNumbers[n];
					Eigen::MatrixXd R = K->getMatrix(tree[nLevels][k].chargeLocations, tree[nLevels][nn].chargeLocations);
					if (nLevels == 1 && n==0) {
						tree[nLevels][k].potential = R*tree[nLevels][nn].charges;
					}
					else if (tree[nLevels][k].charges.size()>0) { // if target points are present
						if (tree[nLevels][k].potential.size()==0) { // if potential is not assigned earlier
							tree[nLevels][k].potential = R*tree[nLevels][nn].charges;
						}
						else {
							tree[nLevels][k].potential += R*tree[nLevels][nn].charges;
						}
					}
			}
			Eigen::MatrixXd R = K->getMatrix(tree[nLevels][k].chargeLocations, tree[nLevels][k].chargeLocations);
			if (tree[nLevels][k].charges.size()>0) {
				if (tree[nLevels][k].potential.size()==0) {
					tree[nLevels][k].potential = R*tree[nLevels][k].charges; //self Interaction
				}
				else {
					tree[nLevels][k].potential += R*tree[nLevels][k].charges; //self Interaction
				}
			}
		}
	}

	// void evaluate_NearField() { // evaluating at chargeLocations
	// 	// #pragma omp parallel for reduction(+:assTime,matVecTime)
	// 	double assTimeLevel = 0.0;
	// 	double matVecTimeLevel = 0.0;
	// 	// #pragma omp parallel
	// 	// {
	// 		double start, end;
	// 		double assTimeThread = 0.0;
	// 		double matVecTimeThread = 0.0;
	// 		// #pragma omp for nowait
	// 		for (size_t k = 0; k < nBoxesPerLevel[nLevels]; k++) {
	// 			for (size_t n = 0; n < tree[nLevels][k].neighborNumbers.size(); n++) {
	// 				int nn = tree[nLevels][k].neighborNumbers[n];
	// 				if(nn != -1) {
	// 					start = omp_get_wtime();
	// 					Eigen::MatrixXd R = K->getMatrix(tree[nLevels][k].chargeLocations, tree[nLevels][nn].chargeLocations);
	// 					end = omp_get_wtime();
	// 					assTimeThread += end-start;
	//
	// 					start = omp_get_wtime();
	// 					tree[nLevels][k].potential += R*tree[nLevels][nn].charges;
	// 					end = omp_get_wtime();
	// 					matVecTimeThread += end-start;
	// 				}
	// 			}
	// 			start = omp_get_wtime();
	// 			Eigen::MatrixXd R = K->getMatrix(tree[nLevels][k].chargeLocations, tree[nLevels][k].chargeLocations);
	// 			end = omp_get_wtime();
	// 			assTimeThread += end-start;
	//
	// 			start = omp_get_wtime();
	// 			tree[nLevels][k].potential += R*tree[nLevels][k].charges; //self Interaction
	// 			end = omp_get_wtime();
	// 			matVecTimeThread += end-start;
	// 		}
	// 		// #pragma omp critical
	// 		// {
	// 			if (matVecTimeLevel < matVecTimeThread) {
	// 				matVecTimeLevel = matVecTimeThread;
	// 			}
	// 			if (assTimeLevel < assTimeThread) {
	// 				assTimeLevel = assTimeThread;
	// 			}
	// 		// }
	// 	// }
	// 	assTime += assTimeLevel;
	// 	matVecTime += matVecTimeLevel;
	// }

	void collectPotential(Eigen::VectorXd &potential) {
		potential = Eigen::VectorXd::Zero(N);
		// std::cout << "N: " << N << std::endl;
		// std::cout << "potential.size(): " << potential.size() << std::endl;
		// #pragma omp parallel for
		for (size_t j = 1; j <= nLevels; j++) {
			int start = 0;
			Eigen::VectorXd potentialTemp = Eigen::VectorXd::Zero(N);
			// #pragma omp parallel for
			for (size_t k = 0; k < nBoxesPerLevel[j]; k++) { //using the fact that all the leaves have same number of particles
				potentialTemp.segment(start, tree[j][k].potential.size()) = tree[j][k].potential;
				start += tree[j][k].potential.size();
			}
			potential = potential + potentialTemp;
		}
		// std::cout << "potential.size(): " << potential.size() << std::endl;
	}

	void reorder(Eigen::VectorXd &potential) {
		Eigen::VectorXd potentialTemp = potential;
		int start = 0;
		for (size_t k = 0; k < nBoxesPerLevel[nLevels]; k++) {
			for (size_t i = 0; i < tree[nLevels][k].chargeLocations.size(); i++) {
				int index = tree[nLevels][k].chargeLocations[i];
				// std::cout << index << std::endl;
				potential(index) = potentialTemp(start);
				start++;
			}
		}
	}

	// void findMemory(double &sum) {
	// 	sum = 0;
	// 	for (size_t j = 2; j <= nLevels; j++) {
	// 		for (size_t k = 0; k < nBoxesPerLevel[j]; k++) {
	// 			for (size_t i = 0; i < tree[j][k].interactionList.size(); i++) {
	// 				int ki = tree[j][k].interactionList[i];
	// 				sum += tree[j][k].L[i].cols()*tree[j][k].chargeLocations.size();
	// 				sum += tree[j][k].L[i].cols()*tree[j][k].L[i].cols();
	// 				sum += tree[j][ki].chargeLocations.size()*tree[j][k].L[i].cols();
	// 			}
	// 		}
	// 	}
	// 	for (size_t k = 0; k < nBoxesPerLevel[nLevels]; k++) {
	// 		sum += tree[nLevels][k].chargeLocations.size()*tree[nLevels][k].chargeLocations.size(); //self
	// 		for (size_t n = 0; n < tree[nLevels][k].neighborNumbers.size(); n++) {
	// 			int nn = tree[nLevels][k].neighborNumbers[n];
	// 			// if(nn != -1) {
	// 				sum += tree[nLevels][k].chargeLocations.size()*tree[nLevels][nn].chargeLocations.size();
	// 			// }
	// 		}
	// 	}
	// }
	//
	// void findMemory2(double &sum) {
	// 	sum = 0;
	// 	for (size_t j = 2; j <= nLevels; j++) {
	// 		for (size_t k = 0; k < nBoxesPerLevel[j]; k++) {
	// 			sum += tree[j][k].L2P.size();
	// 		}
	// 	}
	// 	for (size_t j = 2; j <= nLevels; j++) {
	// 		for (size_t k = 0; k < nBoxesPerLevel[j]; k++) {
	// 			for(int in=0; in<tree[j][k].interactionList.size(); ++in) {
	// 				int kIL = tree[j][k].interactionList[in];
	// 					sum += tree[j][k].M2L[kIL].size();
	// 			}
	// 		}
	// 	}
	//
	// 	for (size_t k = 0; k < nBoxesPerLevel[nLevels]; k++) {
	// 		sum += tree[nLevels][k].chargeLocations.size()*tree[nLevels][k].chargeLocations.size(); //self
	// 		for (size_t n = 0; n < tree[nLevels][k].neighborNumbers.size(); n++) {
	// 			int nn = tree[nLevels][k].neighborNumbers[n];
	// 			// if(nn != -1) {
	// 				sum += tree[nLevels][k].chargeLocations.size()*tree[nLevels][nn].chargeLocations.size();
	// 			// }
	// 		}
	// 	}
	// }

};

#endif
