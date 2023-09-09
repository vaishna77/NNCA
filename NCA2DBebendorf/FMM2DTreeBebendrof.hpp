#ifndef _FMM2DTree_HPP__
#define _FMM2DTree_HPP__
class FMM2DBox {
public:
	int boxNumber;
	int parentNumber;
	int childrenNumbers[4];
	int neighborNumbers[8];
	int innerNumbers[16];
	int outerNumbers[24];

	FMM2DBox () {
		boxNumber		=	-1;
		parentNumber	=	-1;
		for (int l=0; l<4; ++l) {
			childrenNumbers[l]	=	-1;
		}
		for (int l=0; l<8; ++l) {
			neighborNumbers[l]	=	-1;
		}
		for (int l=0; l<16; ++l) {
			innerNumbers[l]		=	-1;
		}
		for (int l=0; l<24; ++l) {
			outerNumbers[l]		=	-1;
		}
	}

	pts2D center;
	Eigen::VectorXd charges;

	std::map<int, Eigen::MatrixXd> M2L;
	Eigen::MatrixXd M2M[4];					//	Transfer from multipoles of 4 children to multipoles of parent.
	//	The following will be stored only at the leaf nodes
	std::vector<pts2D> chebNodes;
	std::vector<pts2D> particle_loc;
  std::vector<int> chargeLocations;

	Eigen::VectorXd outgoing_charges;//equivalent densities {f_{k}^{B,o}}
	Eigen::VectorXd incoming_charges;//equivalent densities {f_{k}^{B,i}}
	Eigen::VectorXd incoming_potential;//check potentials {u_{k}^{B,i}}
	std::vector<int> incoming_chargePoints;//equivalent points {y_{k}^{B,i}}
	std::vector<int> incoming_checkPoints;//check points {x_{k}^{B,i}}
	Eigen::VectorXd potential;//used only at leaf level
	Eigen::MatrixXd L,R;
};

template <typename kerneltype>
class FMM2DTree {
public:
	kerneltype* K;
	int nLevels;			//	Number of levels in the tree.
	int N;					//	Number of particles.
	double L;				//	Semi-length of the simulation box.
	double smallestBoxSize;	//	This is L/2.0^(nLevels).

	std::vector<int> nBoxesPerLevel;			//	Number of boxes at each level in the tree.
	std::vector<double> boxRadius;				//	Box radius at each level in the tree assuming the box at the root is [-1,1]^2
	std::vector<std::vector<FMM2DBox> > tree;	//	The tree storing all the information.

	double ACA_epsilon;
	int nParticlesInLeafAlong1D;
  int nParticlesInLeaf;
  std::vector<double> Nodes1D;
	std::vector<pts2D> Nodes;
	std::vector<double> standardChebNodes1D;
	std::vector<pts2D> standardChebNodes;
  std::vector<pts2D> gridPoints; //all particles in domain
	int TOL_POW;
	int nChebNodes1D;
// public:
	FMM2DTree(kerneltype* K, int nLevels, int nChebNodes1D, int nParticlesInLeafAlong1D, double L, int TOL_POW) {
		this->K					=	K;
		this->nLevels		=	nLevels;
		this->nChebNodes1D		=	nChebNodes1D;
		this->L					=	L;
    this->nParticlesInLeafAlong1D = nParticlesInLeafAlong1D;
    this->nParticlesInLeaf = nParticlesInLeafAlong1D*nParticlesInLeafAlong1D;
		this->TOL_POW = TOL_POW;
    nBoxesPerLevel.push_back(1);
		boxRadius.push_back(L);
		for (int k=1; k<=nLevels; ++k) {
			nBoxesPerLevel.push_back(4*nBoxesPerLevel[k-1]);
			boxRadius.push_back(0.5*boxRadius[k-1]);
		}
		this->smallestBoxSize	=	boxRadius[nLevels];
		// K->a					=	0.01;
		K->a					=	smallestBoxSize;
		this->N					=	nParticlesInLeaf*nBoxesPerLevel[nLevels];
		this->ACA_epsilon = 1.0e-10;
	}

  void set_Standard_Cheb_Nodes() {
		for (int k=0; k<nChebNodes1D; ++k) {
			standardChebNodes1D.push_back(-cos((k+0.5)/nParticlesInLeafAlong1D*PI));
		}
		pts2D temp1;
		for (int j=0; j<nChebNodes1D; ++j) {
			for (int k=0; k<nChebNodes1D; ++k) {
				temp1.x	=	standardChebNodes1D[k];
				temp1.y	=	standardChebNodes1D[j];
				standardChebNodes.push_back(temp1);
			}
		}
	}

	void set_Uniform_Nodes() {
		for (int k=0; k<nParticlesInLeafAlong1D; ++k) {
			Nodes1D.push_back(-1.0+2.0*(k+1.0)/(nParticlesInLeafAlong1D+1.0));
		}
		pts2D temp1;
		for (int j=0; j<nParticlesInLeafAlong1D; ++j) {
			for (int k=0; k<nParticlesInLeafAlong1D; ++k) {
				temp1.x	=	Nodes1D[k];
				temp1.y	=	Nodes1D[j];
				Nodes.push_back(temp1);
			}
		}
	}

  void shift_Nodes(double radius, pts2D center, std::vector<pts2D> &particle_loc) {
    for (size_t i = 0; i < Nodes.size(); i++) {
      pts2D temp;
      temp.x = Nodes[i].x*radius + center.x;
      temp.y = Nodes[i].y*radius + center.y;
      particle_loc.push_back(temp);
    }
  }

	void shift_Cheb_Nodes(pts2D radius, pts2D center, std::vector<pts2D> &particle_loc) {
    for (size_t i = 0; i < nChebNodes1D*nChebNodes1D; i++) {
      pts2D temp;
      temp.x = standardChebNodes[i].x*radius.x + center.x;
      temp.y = standardChebNodes[i].y*radius.y + center.y;
      particle_loc.push_back(temp);
    }
  }

	void createTree() {
		//	First create root and add to tree
		FMM2DBox root;
		root.boxNumber		=	0;
		root.parentNumber	=	-1;
		#pragma omp parallel for
		for (int l=0; l<4; ++l) {
			root.childrenNumbers[l]	=	l;
		}
		#pragma omp parallel for
		for (int l=0; l<8; ++l) {
			root.neighborNumbers[l]	=	-1;
		}
		#pragma omp parallel for
		for (int l=0; l<16; ++l) {
			root.innerNumbers[l]	=	-1;
		}
		#pragma omp parallel for
		for (int l=0; l<24; ++l) {
			root.outerNumbers[l]	=	-1;
		}
		std::vector<FMM2DBox> rootLevel;
		rootLevel.push_back(root);
		tree.push_back(rootLevel);

		for (int j=1; j<=nLevels; ++j) {
			std::vector<FMM2DBox> level;
			for (int k=0; k<nBoxesPerLevel[j]; ++k) {
				FMM2DBox box;
				box.boxNumber		=	k;
				box.parentNumber	=	k/4;
				for (int l=0; l<4; ++l) {
					box.childrenNumbers[l]	=	4*k+l;
				}
				level.push_back(box);
			}
			tree.push_back(level);
		}
	}

	//	Assigns the interactions for child0 of a box
	void assign_Child0_Interaction(int j, int k) {
		int nL	=	j+1;
		int nC	=	4*k;
		int nN;

		//	Assign siblings
		{
			/************************/
			/*	 _____________		*/
			/*	|	   |	  |		*/
			/*	|  N5  |  N4  |		*/
			/*	|______|______|		*/
			/*	|	   |	  |		*/
			/*	|  **  |  N3  |		*/
			/*	|______|______|		*/
			/*						*/
			/************************/
			tree[nL][nC].neighborNumbers[3]	=	nC+1;
			tree[nL][nC].neighborNumbers[4]	=	nC+2;
			tree[nL][nC].neighborNumbers[5]	=	nC+3;
		}

		//	Assign children of parent's zeroth neighbor
		{
			/****************************/
			/*				   ______	*/
			/*				  |		 |	*/
			/*				  |	 **  |	*/
			/*	 _____________|______|  */
			/*	|	   |	  |			*/
			/*	|  I15 |  N0  |			*/
			/*	|______|______|			*/
			/*	|	   |	  |			*/
			/*	|  I0  |  I1  |			*/
			/*	|______|______|			*/
			/*							*/
			/****************************/
			nN	=	tree[j][k].neighborNumbers[0];
			if (nN != -1) {
				tree[nL][nC].innerNumbers[0]	=	tree[j][nN].childrenNumbers[0];
				tree[nL][nC].innerNumbers[1]	=	tree[j][nN].childrenNumbers[1];
				tree[nL][nC].neighborNumbers[0]	=	tree[j][nN].childrenNumbers[2];
				tree[nL][nC].innerNumbers[15]	=	tree[j][nN].childrenNumbers[3];
			}
		}

		//	Assign children of parent's first neighbor
		{
			/****************************/
			/*	 ______			  		*/
			/*	|	   |	  			*/
			/*	|  **  |				*/
			/*	|______|______			*/
			/*	|	   |	  |			*/
			/*	|  N1  |  N2  |			*/
			/*	|______|______|			*/
			/*	|	   |	  |			*/
			/*	|  I2  |  I3  |			*/
			/*	|______|______|			*/
			/*							*/
			/****************************/
			nN	=	tree[j][k].neighborNumbers[1];
			if (nN != -1) {
				tree[nL][nC].innerNumbers[2]	=	tree[j][nN].childrenNumbers[0];
				tree[nL][nC].innerNumbers[3]	=	tree[j][nN].childrenNumbers[1];
				tree[nL][nC].neighborNumbers[2]	=	tree[j][nN].childrenNumbers[2];
				tree[nL][nC].neighborNumbers[1]	=	tree[j][nN].childrenNumbers[3];
			}
		}

		//	Assign children of parent's second neighbor
		{
			/************************************/
			/*	 ______			  				*/
			/*	|	   |	  					*/
			/*	|  **  |						*/
			/*	|______|	   _____________	*/
			/*				  |	     |	    |	*/
			/*				  |  I5  |  O8  |	*/
			/*				  |______|______|	*/
			/*				  |	     |	    |	*/
			/*				  |  I4  |  O7  |	*/
			/*				  |______|______|	*/
			/*									*/
			/************************************/
			nN	=	tree[j][k].neighborNumbers[2];
			if (nN != -1) {
				tree[nL][nC].innerNumbers[4]	=	tree[j][nN].childrenNumbers[0];
				tree[nL][nC].outerNumbers[7]	=	tree[j][nN].childrenNumbers[1];
				tree[nL][nC].outerNumbers[8]	=	tree[j][nN].childrenNumbers[2];
				tree[nL][nC].innerNumbers[5]	=	tree[j][nN].childrenNumbers[3];
			}
		}

		//	Assign children of parent's third neighbor
		{
			/************************************/
			/*				   _____________	*/
			/*				  |	     |	    |	*/
			/*				  |  I7  |  O10 |	*/
			/*	 ______		  |______|______|	*/
			/*	|	   |	  |	     |	    |	*/
			/*	|  **  |	  |  I6  |  O9  |	*/
			/*	|______|	  |______|______|	*/
			/*									*/
			/************************************/
			nN	=	tree[j][k].neighborNumbers[3];
			if (nN!=-1) {
				tree[nL][nC].innerNumbers[6]	=	tree[j][nN].childrenNumbers[0];
				tree[nL][nC].outerNumbers[9]	=	tree[j][nN].childrenNumbers[1];
				tree[nL][nC].outerNumbers[10]	=	tree[j][nN].childrenNumbers[2];
				tree[nL][nC].innerNumbers[7]	=	tree[j][nN].childrenNumbers[3];
			}

		}

		//	Assign children of parent's fourth neighbor
		{
			/************************************/
			/*				   _____________	*/
			/*				  |	     |	    |	*/
			/*				  |  O13 |  O12 |	*/
			/*				  |______|______|	*/
			/*				  |	     |	    |	*/
			/*		    	  |  I8  |  O11 |	*/
			/*		    	  |______|______|	*/
			/*									*/
			/*									*/
			/*	 ______							*/
			/*  |      |						*/
			/*  |  **  |						*/
			/*  |______|						*/
			/*  								*/
			/************************************/
			nN	=	tree[j][k].neighborNumbers[4];
			if (nN != -1) {
				tree[nL][nC].innerNumbers[8]	=	tree[j][nN].childrenNumbers[0];
				tree[nL][nC].outerNumbers[11]	=	tree[j][nN].childrenNumbers[1];
				tree[nL][nC].outerNumbers[12]	=	tree[j][nN].childrenNumbers[2];
				tree[nL][nC].outerNumbers[13]	=	tree[j][nN].childrenNumbers[3];
			}
		}

		//	Assign children of parent's fifth neighbor
		{
			/************************/
			/*	 _____________		*/
			/*	|	   |	  |		*/
			/*	|  O15 |  O14 |		*/
			/*	|______|______|		*/
			/*	|	   |	  |		*/
			/*	|  I10 |  I9  |		*/
			/*	|______|______|		*/
			/*						*/
			/*						*/
			/*	 ______				*/
			/*  |	   |			*/
			/*	|  **  |			*/
			/*	|______|			*/
			/*  					*/
			/************************/
			nN	=	tree[j][k].neighborNumbers[5];
			if (nN != -1) {
				tree[nL][nC].innerNumbers[10]	=	tree[j][nN].childrenNumbers[0];
				tree[nL][nC].innerNumbers[9]	=	tree[j][nN].childrenNumbers[1];
				tree[nL][nC].outerNumbers[14]	=	tree[j][nN].childrenNumbers[2];
				tree[nL][nC].outerNumbers[15]	=	tree[j][nN].childrenNumbers[3];
			}
		}

		//	Assign children of parent's sixth neighbor
		{
			/****************************/
			/*	 _____________			*/
			/*	|	   |	  |			*/
			/*	|  O17 |  O16 |			*/
			/*	|______|______|			*/
			/*	|	   |	  |			*/
			/*	|  I12 |  I11 |			*/
			/*	|______|______|			*/
			/*							*/
			/*							*/
			/*				   ______	*/
			/*  			  |		 |	*/
			/*				  |	 **  |	*/
			/*				  |______|	*/
			/*  						*/
			/****************************/
			nN	=	tree[j][k].neighborNumbers[6];
			if (nN != -1) {
				tree[nL][nC].innerNumbers[12]	=	tree[j][nN].childrenNumbers[0];
				tree[nL][nC].innerNumbers[11]	=	tree[j][nN].childrenNumbers[1];
				tree[nL][nC].outerNumbers[16]	=	tree[j][nN].childrenNumbers[2];
				tree[nL][nC].outerNumbers[17]	=	tree[j][nN].childrenNumbers[3];
			}
		}

		//	Assign children of parent's seventh neighbor
		{
			/****************************/
			/*							*/
			/*	 _____________			*/
			/*	|	   |	  |			*/
			/*	|  I13 |  N6  |			*/
			/*	|______|______|______	*/
			/*  |	   |	  |		 |	*/
			/*	|  I14 |  N7  |	 **  |	*/
			/*	|______|______|______|	*/
			/*  						*/
			/****************************/
			nN	=	tree[j][k].neighborNumbers[7];
			if (nN != -1) {
				tree[nL][nC].innerNumbers[14]	=	tree[j][nN].childrenNumbers[0];
				tree[nL][nC].neighborNumbers[7]	=	tree[j][nN].childrenNumbers[1];
				tree[nL][nC].neighborNumbers[6]	=	tree[j][nN].childrenNumbers[2];
				tree[nL][nC].innerNumbers[13]	=	tree[j][nN].childrenNumbers[3];
			}
		}
	}

	//	Assigns the interactions for child1 of a box
	void assign_Child1_Interaction(int j, int k) {
		int nL	=	j+1;
		int nC	=	4*k+1;
		int nN;

		//	Assign siblings
		{
			/************************/
			/*	 _____________		*/
			/*	|	   |	  |		*/
			/*	|  N6  |  N5  |		*/
			/*	|______|______|		*/
			/*	|	   |	  |		*/
			/*	|  N7  |  **  |		*/
			/*	|______|______|		*/
			/*						*/
			/************************/
			tree[nL][nC].neighborNumbers[7]	=	nC-1;
			tree[nL][nC].neighborNumbers[5]	=	nC+1;
			tree[nL][nC].neighborNumbers[6]	=	nC+2;
		}

		//	Assign children of parent's zeroth neighbor
		{
			/************************************/
			/*				   		  ______	*/
			/*				  	     |		|	*/
			/*				         |	**  |	*/
			/*	 _____________       |______|  	*/
			/*	|	   |	  |					*/
			/*	|  O22 |  I15 |					*/
			/*	|______|______|					*/
			/*	|	   |	  |					*/
			/*	|  O23 |  I0  |					*/
			/*	|______|______|					*/
			/*									*/
			/************************************/
			nN	=	tree[j][k].neighborNumbers[0];
			if (nN != -1) {
				tree[nL][nC].outerNumbers[23]	=	tree[j][nN].childrenNumbers[0];
				tree[nL][nC].innerNumbers[0]	=	tree[j][nN].childrenNumbers[1];
				tree[nL][nC].innerNumbers[15]	=	tree[j][nN].childrenNumbers[2];
				tree[nL][nC].outerNumbers[22]	=	tree[j][nN].childrenNumbers[3];
			}
		}

		//	Assign children of parent's first neighbor
		{
			/****************************/
			/*	 		______		  	*/
			/*		   |	  |			*/
			/*	       |  **  |			*/
			/*	 ______|______|			*/
			/*	|	   |	  |			*/
			/*	|  N0  |  N1  |			*/
			/*	|______|______|			*/
			/*	|	   |	  |			*/
			/*	|  I1  |  I2  |			*/
			/*	|______|______|			*/
			/*							*/
			/****************************/
			nN	=	tree[j][k].neighborNumbers[1];
			if (nN != -1) {
				tree[nL][nC].innerNumbers[1]	=	tree[j][nN].childrenNumbers[0];
				tree[nL][nC].innerNumbers[2]	=	tree[j][nN].childrenNumbers[1];
				tree[nL][nC].neighborNumbers[1]	=	tree[j][nN].childrenNumbers[2];
				tree[nL][nC].neighborNumbers[0]	=	tree[j][nN].childrenNumbers[3];
			}
		}

		//	Assign children of parent's second neighbor
		{
			/****************************/
			/*	 ______		  			*/
			/*	|	   |				*/
			/*	|  **  |	  			*/
			/*	|______|_____________	*/
			/*		   |	  |	     |	*/
			/*		   |  N2  |  I5  |	*/
			/*		   |______|______|	*/
			/*		   |	  |	     |	*/
			/*	       |  I3  |  I4  |	*/
			/*		   |______|______|	*/
			/*							*/
			/****************************/
			nN	=	tree[j][k].neighborNumbers[2];
			if (nN != -1) {
				tree[nL][nC].innerNumbers[3]	=	tree[j][nN].childrenNumbers[0];
				tree[nL][nC].innerNumbers[4]	=	tree[j][nN].childrenNumbers[1];
				tree[nL][nC].innerNumbers[5]	=	tree[j][nN].childrenNumbers[2];
				tree[nL][nC].neighborNumbers[2]	=	tree[j][nN].childrenNumbers[3];
			}
		}

		//	Assign children of parent's third neighbor
		{
			/****************************/
			/*			_____________	*/
			/*		   |	  |	     |	*/
			/*		   |  N4  |	 I7	 |  */
			/*	 ______|______|______|	*/
			/*	|	   |	  |	     |	*/
			/*	|  **  |  N3  |  I6  |  */
			/*	|______|______|______| 	*/
			/*							*/
			/****************************/
			nN	=	tree[j][k].neighborNumbers[3];
			if (nN != -1){
				tree[nL][nC].neighborNumbers[3]	=	tree[j][nN].childrenNumbers[0];
				tree[nL][nC].innerNumbers[6]	=	tree[j][nN].childrenNumbers[1];
				tree[nL][nC].innerNumbers[7]	=	tree[j][nN].childrenNumbers[2];
				tree[nL][nC].neighborNumbers[4]	=	tree[j][nN].childrenNumbers[3];
			}
		}

		//	Assign children of parent's fourth neighbor
		{
			/****************************/
			/*			_____________	*/
			/*		   |	  |	     |	*/
			/*		   |  O14 |  O13 |	*/
			/*		   |______|______|	*/
			/*		   |	  |	     |	*/
			/*		   |  I9  |  I8  |	*/
			/*		   |______|______|	*/
			/*				  			*/
			/*				  			*/
			/*	 ______					*/
			/*	|	   |				*/
			/*  |  **  |				*/
			/*  |______|				*/
			/*  						*/
			/****************************/
			nN	=	tree[j][k].neighborNumbers[4];
			if (nN != -1) {
				tree[nL][nC].innerNumbers[9]	=	tree[j][nN].childrenNumbers[0];
				tree[nL][nC].innerNumbers[8]	=	tree[j][nN].childrenNumbers[1];
				tree[nL][nC].outerNumbers[13]	=	tree[j][nN].childrenNumbers[2];
				tree[nL][nC].outerNumbers[14]	=	tree[j][nN].childrenNumbers[3];
			}
		}

		//	Assign children of parent's fifth neighbor
		{
			/************************/
			/*	 _____________		*/
			/*	|	   |	  |		*/
			/*	|  O16 |  O15 |		*/
			/*	|______|______|		*/
			/*	|	   |	  |		*/
			/*	|  I11 |  I10 |		*/
			/*	|______|______|		*/
			/*						*/
			/*						*/
			/*		    ______		*/
			/* 		   |	  |		*/
			/*		   |  **  |		*/
			/*		   |______|		*/
			/*  					*/
			/************************/
			nN	=	tree[j][k].neighborNumbers[5];
			if (nN != -1) {
				tree[nL][nC].innerNumbers[11]	=	tree[j][nN].childrenNumbers[0];
				tree[nL][nC].innerNumbers[10]	=	tree[j][nN].childrenNumbers[1];
				tree[nL][nC].outerNumbers[15]	=	tree[j][nN].childrenNumbers[2];
				tree[nL][nC].outerNumbers[16]	=	tree[j][nN].childrenNumbers[3];
			}
		}

		//	Assign children of parent's sixth neighbor
		{
			/************************************/
			/*	 _____________					*/
			/*	|	   |	  |					*/
			/*	|  O18 |  O17 |					*/
			/*	|______|______|					*/
			/*	|	   |	  |					*/
			/*	|  O19 |  I12 |					*/
			/*	|______|______|					*/
			/*									*/
			/*									*/
			/*				   		  ______	*/
			/*  			  		 |		|	*/
			/*				  		 |	** 	|	*/
			/*				  		 |______|	*/
			/*  								*/
			/************************************/
			nN	=	tree[j][k].neighborNumbers[6];
			if (nN != -1) {
				tree[nL][nC].outerNumbers[19]	=	tree[j][nN].childrenNumbers[0];
				tree[nL][nC].innerNumbers[12]	=	tree[j][nN].childrenNumbers[1];
				tree[nL][nC].outerNumbers[17]	=	tree[j][nN].childrenNumbers[2];
				tree[nL][nC].outerNumbers[18]	=	tree[j][nN].childrenNumbers[3];
			}
		}

		//	Assign children of parent's seventh neighbor
		{
			/************************************/
			/*									*/
			/*	 _____________					*/
			/*	|	   |	  |					*/
			/*	|  O20 |  I13 |					*/
			/*	|______|______|		  ______	*/
			/*  |	   |	  |		 |		|	*/
			/*	|  O21 |  I14 |	 	 |	**  |	*/
			/*	|______|______|		 |______|	*/
			/*  								*/
			/************************************/
			nN	=	tree[j][k].neighborNumbers[7];
			if (nN != -1) {
				tree[nL][nC].outerNumbers[21]	=	tree[j][nN].childrenNumbers[0];
				tree[nL][nC].innerNumbers[14]	=	tree[j][nN].childrenNumbers[1];
				tree[nL][nC].innerNumbers[13]	=	tree[j][nN].childrenNumbers[2];
				tree[nL][nC].outerNumbers[20]	=	tree[j][nN].childrenNumbers[3];
			}
		}
	}

	//	Assigns the interactions for child2 of a box
	void assign_Child2_Interaction(int j, int k) {
		int nL	=	j+1;
		int nC	=	4*k+2;
		int nN;

		//	Assign siblings
		{
			/************************/
			/*	 _____________		*/
			/*	|	   |	  |		*/
			/*	|  N7  |  **  |		*/
			/*	|______|______|		*/
			/*	|	   |	  |		*/
			/*	|  N0  |  N1  |		*/
			/*	|______|______|		*/
			/*						*/
			/************************/
			tree[nL][nC].neighborNumbers[0]	=	nC-2;
			tree[nL][nC].neighborNumbers[1]	=	nC-1;
			tree[nL][nC].neighborNumbers[7]	=	nC+1;
		}

		//	Assign children of parent's zeroth neighbor
		{
			/************************************/
			/*				   		  ______	*/
			/*				  	     |		|	*/
			/*				         |	**  |	*/
			/*				         |______|  	*/
			/*									*/
			/*									*/
			/*	 _____________					*/
			/*	|	   |	  |					*/
			/*	|  O23 |  I0  |					*/
			/*	|______|______|					*/
			/*	|	   |	  |					*/
			/*	|  O0  |  O1  |					*/
			/*	|______|______|					*/
			/*									*/
			/************************************/
			nN	=	tree[j][k].neighborNumbers[0];
			if (nN != -1) {
				tree[nL][nC].outerNumbers[0]	=	tree[j][nN].childrenNumbers[0];
				tree[nL][nC].outerNumbers[1]	=	tree[j][nN].childrenNumbers[1];
				tree[nL][nC].innerNumbers[0]	=	tree[j][nN].childrenNumbers[2];
				tree[nL][nC].outerNumbers[23]	=	tree[j][nN].childrenNumbers[3];
			}
		}

		//	Assign children of parent's first neighbor
		{
			/****************************/
			/*	 		______		  	*/
			/*		   |	  |			*/
			/*	       |  **  |			*/
			/*	 	   |______|			*/
			/*							*/
			/*							*/
			/*	 _____________			*/
			/*	|	   |	  |			*/
			/*	|  I1  |  I2  |			*/
			/*	|______|______|			*/
			/*	|	   |	  |			*/
			/*	|  O2  |  O3  |			*/
			/*	|______|______|			*/
			/*							*/
			/****************************/
			nN	=	tree[j][k].neighborNumbers[1];
			if (nN != -1) {
				tree[nL][nC].outerNumbers[2]	=	tree[j][nN].childrenNumbers[0];
				tree[nL][nC].outerNumbers[3]	=	tree[j][nN].childrenNumbers[1];
				tree[nL][nC].innerNumbers[2]	=	tree[j][nN].childrenNumbers[2];
				tree[nL][nC].innerNumbers[1]	=	tree[j][nN].childrenNumbers[3];
			}
		}

		//	Assign children of parent's second neighbor
		{
			/****************************/
			/*	 ______		  			*/
			/*	|	   |				*/
			/*	|  **  |	  			*/
			/*	|______|				*/
			/*							*/
			/*							*/
			/*			_____________	*/
			/*		   |	  |	     |	*/
			/*		   |  I3  |  I4  |	*/
			/*		   |______|______|	*/
			/*		   |	  |	     |	*/
			/*	       |  O4  |  O5  |	*/
			/*		   |______|______|	*/
			/*							*/
			/****************************/
			nN	=	tree[j][k].neighborNumbers[2];
			if (nN != -1) {
				tree[nL][nC].outerNumbers[4]	=	tree[j][nN].childrenNumbers[0];
				tree[nL][nC].outerNumbers[5]	=	tree[j][nN].childrenNumbers[1];
				tree[nL][nC].innerNumbers[4]	=	tree[j][nN].childrenNumbers[2];
				tree[nL][nC].innerNumbers[3]	=	tree[j][nN].childrenNumbers[3];
			}
		}

		//	Assign children of parent's third neighbor
		{
			/****************************/
			/*	 ____________________	*/
			/*	|	   |	  |	     |	*/
			/*	|  **  |  N3  |	 I6	 |  */
			/*	|______|______|______|	*/
			/*		   |	  |	     |	*/
			/*		   |  N2  |  I5  |  */
			/*		   |______|______| 	*/
			/*							*/
			/****************************/
			nN	=	tree[j][k].neighborNumbers[3];
			if (nN != -1) {
				tree[nL][nC].neighborNumbers[2]	=	tree[j][nN].childrenNumbers[0];
				tree[nL][nC].innerNumbers[5]	=	tree[j][nN].childrenNumbers[1];
				tree[nL][nC].innerNumbers[6]	=	tree[j][nN].childrenNumbers[2];
				tree[nL][nC].neighborNumbers[3]	=	tree[j][nN].childrenNumbers[3];
			}
		}

		//	Assign children of parent's fourth neighbor
		{
			/****************************/
			/*			_____________	*/
			/*		   |	  |	     |	*/
			/*		   |  I9  |  I8  |	*/
			/*		   |______|______|	*/
			/*		   |	  |	     |	*/
			/*		   |  N4  |  I7  |	*/
			/*	 ______|______|______|	*/
			/*	|	   |	  			*/
			/*	|  **  |	  			*/
			/*	|______|	  			*/
			/*							*/
			/****************************/
			nN	=	tree[j][k].neighborNumbers[4];
			if (nN != -1) {
				tree[nL][nC].neighborNumbers[4]	=	tree[j][nN].childrenNumbers[0];
				tree[nL][nC].innerNumbers[7]	=	tree[j][nN].childrenNumbers[1];
				tree[nL][nC].innerNumbers[8]	=	tree[j][nN].childrenNumbers[2];
				tree[nL][nC].innerNumbers[9]	=	tree[j][nN].childrenNumbers[3];
			}
		}

		//	Assign children of parent's fifth neighbor
		{
			/************************/
			/*	 _____________		*/
			/*	|	   |	  |		*/
			/*	|  I11 |  I10 |		*/
			/*	|______|______|		*/
			/*	|	   |	  |		*/
			/*	|  N6  |  N5  |		*/
			/*	|______|______|		*/
			/*		   |	  |		*/
			/*		   |  **  |		*/
			/*		   |______|		*/
			/*  					*/
			/************************/
			nN	=	tree[j][k].neighborNumbers[5];
			if (nN != -1) {
				tree[nL][nC].neighborNumbers[6]	=	tree[j][nN].childrenNumbers[0];
				tree[nL][nC].neighborNumbers[5]	=	tree[j][nN].childrenNumbers[1];
				tree[nL][nC].innerNumbers[10]	=	tree[j][nN].childrenNumbers[2];
				tree[nL][nC].innerNumbers[11]	=	tree[j][nN].childrenNumbers[3];
			}
		}

		//	Assign children of parent's sixth neighbor
		{
			/************************************/
			/*	 _____________					*/
			/*	|	   |	  |					*/
			/*	|  O19 |  I12 |					*/
			/*	|______|______|					*/
			/*	|	   |	  |					*/
			/*	|  O20 |  I13 |					*/
			/*	|______|______|		  ______	*/
			/*  			  		 |		|	*/
			/*				  		 |	** 	|	*/
			/*				  		 |______|	*/
			/*  								*/
			/************************************/
			nN	=	tree[j][k].neighborNumbers[6];
			if (nN != -1) {
				tree[nL][nC].outerNumbers[20]	=	tree[j][nN].childrenNumbers[0];
				tree[nL][nC].innerNumbers[13]	=	tree[j][nN].childrenNumbers[1];
				tree[nL][nC].innerNumbers[12]	=	tree[j][nN].childrenNumbers[2];
				tree[nL][nC].outerNumbers[19]	=	tree[j][nN].childrenNumbers[3];
			}
		}

		//	Assign children of parent's seventh neighbor
		{
			/************************************/
			/*									*/
			/*	 _____________		  ______	*/
			/*	|	   |	  |		 |	    |	*/
			/*	|  O21 |  I14 |		 |	**	|	*/
			/*	|______|______|		 |______|	*/
			/*  |	   |	  |		 			*/
			/*	|  O22 |  I15 |	 	 			*/
			/*	|______|______|		 			*/
			/*  								*/
			/************************************/
			nN	=	tree[j][k].neighborNumbers[7];
			if (nN != -1) {
				tree[nL][nC].outerNumbers[22]	=	tree[j][nN].childrenNumbers[0];
				tree[nL][nC].innerNumbers[15]	=	tree[j][nN].childrenNumbers[1];
				tree[nL][nC].innerNumbers[14]	=	tree[j][nN].childrenNumbers[2];
				tree[nL][nC].outerNumbers[21]	=	tree[j][nN].childrenNumbers[3];
			}
		}
	}

	//	Assigns the interactions for child3 of a box
	void assign_Child3_Interaction(int j, int k) {
		int nL	=	j+1;
		int nC	=	4*k+3;
		int nN;

		//	Assign siblings
		{
			/************************/
			/*	 _____________		*/
			/*	|	   |	  |		*/
			/*	|  **  |  N3  |		*/
			/*	|______|______|		*/
			/*	|	   |	  |		*/
			/*	|  N1  |  N2  |		*/
			/*	|______|______|		*/
			/*						*/
			/************************/
			tree[nL][nC].neighborNumbers[1]	=	nC-3;
			tree[nL][nC].neighborNumbers[2]	=	nC-2;
			tree[nL][nC].neighborNumbers[3]	=	nC-1;
		}

		//	Assign children of parent's zeroth neighbor
		{
			/****************************/
			/*				   ______	*/
			/*				  |		 |	*/
			/*				  |	 **  |	*/
			/*				  |______|  */
			/*							*/
			/*							*/
			/*	 _____________			*/
			/*	|	   |	  |			*/
			/*	|  I0  |  I1  |			*/
			/*	|______|______|			*/
			/*	|	   |	  |			*/
			/*	|  O1  |  O2  |			*/
			/*	|______|______|			*/
			/*							*/
			/****************************/
			nN	=	tree[j][k].neighborNumbers[0];
			if (nN != -1) {
				tree[nL][nC].outerNumbers[1]	=	tree[j][nN].childrenNumbers[0];
				tree[nL][nC].outerNumbers[2]	=	tree[j][nN].childrenNumbers[1];
				tree[nL][nC].innerNumbers[1]	=	tree[j][nN].childrenNumbers[2];
				tree[nL][nC].innerNumbers[0]	=	tree[j][nN].childrenNumbers[3];
			}
		}

		//	Assign children of parent's first neighbor
		{
			/****************************/
			/*	 ______		  			*/
			/*	|	   |				*/
			/*	|  **  |				*/
			/*	|______|				*/
			/*							*/
			/*							*/
			/*	 _____________			*/
			/*	|	   |	  |			*/
			/*	|  I2  |  I3  |			*/
			/*	|______|______|			*/
			/*	|	   |	  |			*/
			/*	|  O3  |  O4  |			*/
			/*	|______|______|			*/
			/*							*/
			/****************************/
			nN	=	tree[j][k].neighborNumbers[1];
			if (nN != -1) {
				tree[nL][nC].outerNumbers[3]	=	tree[j][nN].childrenNumbers[0];
				tree[nL][nC].outerNumbers[4]	=	tree[j][nN].childrenNumbers[1];
				tree[nL][nC].innerNumbers[3]	=	tree[j][nN].childrenNumbers[2];
				tree[nL][nC].innerNumbers[2]	=	tree[j][nN].childrenNumbers[3];
			}
		}

		//	Assign children of parent's second neighbor
		{
			/************************************/
			/*	 ______		  					*/
			/*	|	   |						*/
			/*	|  **  |	  					*/
			/*	|______|						*/
			/*									*/
			/*									*/
			/*				   _____________	*/
			/*		   		  |	     |	    |	*/
			/*		   		  |  I4  |  O7  |	*/
			/*		   		  |______|______|	*/
			/*		   		  |	  	 |	    |	*/
			/*	       		  |  O5  |  O6  |	*/
			/*		   		  |______|______|	*/
			/*									*/
			/************************************/
			nN	=	tree[j][k].neighborNumbers[2];
			if (nN != -1) {
				tree[nL][nC].outerNumbers[5]	=	tree[j][nN].childrenNumbers[0];
				tree[nL][nC].outerNumbers[6]	=	tree[j][nN].childrenNumbers[1];
				tree[nL][nC].outerNumbers[7]	=	tree[j][nN].childrenNumbers[2];
				tree[nL][nC].innerNumbers[4]	=	tree[j][nN].childrenNumbers[3];
			}
		}

		//	Assign children of parent's third neighbor
		{
			/************************************/
			/*	 ______		   _____________	*/
			/*	|	   |	  |	     |		|	*/
			/*	|  **  |      |	 I6	 |  O9	|	*/
			/*	|______|	  |______|______|	*/
			/*		   		  |	  	 |	    |	*/
			/*		   		  |  I5  |  O8  |  	*/
			/*		   		  |______|______| 	*/
			/*									*/
			/************************************/
			nN	=	tree[j][k].neighborNumbers[3];
			if (nN != -1) {
				tree[nL][nC].innerNumbers[5]	=	tree[j][nN].childrenNumbers[0];
				tree[nL][nC].outerNumbers[8]	=	tree[j][nN].childrenNumbers[1];
				tree[nL][nC].outerNumbers[9]	=	tree[j][nN].childrenNumbers[2];
				tree[nL][nC].innerNumbers[6]	=	tree[j][nN].childrenNumbers[3];
			}
		}

		//	Assign children of parent's fourth neighbor
		{
			/************************************/
			/*				   _____________	*/
			/*		   		  |	  	 |	    |	*/
			/*		   		  |  I8  |  O11 |	*/
			/*		   		  |______|______|	*/
			/*		   		  |	  	 |	    |	*/
			/*		   		  |  I7  |  O10 |	*/
			/*	 ______	      |______|______|	*/
			/*	|	   |	  					*/
			/*	|  **  |	  					*/
			/*	|______|	  					*/
			/*									*/
			/************************************/
			nN	=	tree[j][k].neighborNumbers[4];
			if (nN != -1) {
				tree[nL][nC].innerNumbers[7]	=	tree[j][nN].childrenNumbers[0];
				tree[nL][nC].outerNumbers[10]	=	tree[j][nN].childrenNumbers[1];
				tree[nL][nC].outerNumbers[11]	=	tree[j][nN].childrenNumbers[2];
				tree[nL][nC].innerNumbers[8]	=	tree[j][nN].childrenNumbers[3];
			}
		}

		//	Assign children of parent's fifth neighbor
		{
			/************************/
			/*	 _____________		*/
			/*	|	   |	  |		*/
			/*	|  I10 |  I9  |		*/
			/*	|______|______|		*/
			/*	|	   |	  |		*/
			/*	|  N5  |  N4  |		*/
			/*	|______|______|		*/
			/*	|	   |			*/
			/*	|  **  |			*/
			/*	|______|			*/
			/*  					*/
			/************************/
			nN	=	tree[j][k].neighborNumbers[5];
			if (nN != -1) {
				tree[nL][nC].neighborNumbers[5]	=	tree[j][nN].childrenNumbers[0];
				tree[nL][nC].neighborNumbers[4]	=	tree[j][nN].childrenNumbers[1];
				tree[nL][nC].innerNumbers[9]	=	tree[j][nN].childrenNumbers[2];
				tree[nL][nC].innerNumbers[10]	=	tree[j][nN].childrenNumbers[3];
			}
		}

		//	Assign children of parent's sixth neighbor
		{
			/****************************/
			/*	 _____________			*/
			/*	|	   |	  |			*/
			/*	|  I12 |  I11 |			*/
			/*	|______|______|			*/
			/*	|	   |	  |			*/
			/*	|  I13 |  N6  |			*/
			/*	|______|______|______	*/
			/*  			  |		 |	*/
			/*				  |	 **  |	*/
			/*				  |______|	*/
			/*  						*/
			/****************************/
			nN	=	tree[j][k].neighborNumbers[6];
			if (nN != -1) {
				tree[nL][nC].innerNumbers[13]	=	tree[j][nN].childrenNumbers[0];
				tree[nL][nC].neighborNumbers[6]	=	tree[j][nN].childrenNumbers[1];
				tree[nL][nC].innerNumbers[11]	=	tree[j][nN].childrenNumbers[2];
				tree[nL][nC].innerNumbers[12]	=	tree[j][nN].childrenNumbers[3];
			}
		}

		//	Assign children of parent's seventh neighbor
		{
			/****************************/
			/*							*/
			/*	 ____________________	*/
			/*	|	   |	  |		 |	*/
			/*	|  I14 |  N7  |	 **	 |	*/
			/*	|______|______|______|	*/
			/*  |	   |	  |		 	*/
			/*	|  I15 |  N0  |	 	 	*/
			/*	|______|______|		 	*/
			/*  						*/
			/****************************/
			nN	=	tree[j][k].neighborNumbers[7];
			if (nN != -1) {
				tree[nL][nC].innerNumbers[15]	=	tree[j][nN].childrenNumbers[0];
				tree[nL][nC].neighborNumbers[0]	=	tree[j][nN].childrenNumbers[1];
				tree[nL][nC].neighborNumbers[7]	=	tree[j][nN].childrenNumbers[2];
				tree[nL][nC].innerNumbers[14]	=	tree[j][nN].childrenNumbers[3];
			}
		}
	}


	//	Assigns the interactions for the children of a box
	void assign_Box_Interactions(int j, int k) {
		assign_Child0_Interaction(j,k);
		assign_Child1_Interaction(j,k);
		assign_Child2_Interaction(j,k);
		assign_Child3_Interaction(j,k);
	}

	//	Assigns the interactions for the children all boxes at a given level
	void assign_Level_Interactions(int j) {
		#pragma omp parallel for
		for (int k=0; k<nBoxesPerLevel[j]; ++k) {
			assign_Box_Interactions(j,k);
		}
	}

	//	Assigns the interactions for the children all boxes in the tree
	void assign_Tree_Interactions() {
		for (int j=0; j<nLevels; ++j) {
			assign_Level_Interactions(j);
		}
	}

	void assign_Center_Location() {
		int J;
		tree[0][0].center.x	=	0.0;
		tree[0][0].center.y	=	0.0;
		for (int j=0; j<nLevels; ++j) {
			J	=	j+1;
			double shift	=	0.5*boxRadius[j];
			#pragma omp parallel for
			for (int k=0; k<nBoxesPerLevel[j]; ++k) {
				tree[J][4*k].center.x		=	tree[j][k].center.x-shift;
				tree[J][4*k+1].center.x	=	tree[j][k].center.x+shift;
				tree[J][4*k+2].center.x	=	tree[j][k].center.x+shift;
				tree[J][4*k+3].center.x	=	tree[j][k].center.x-shift;

				tree[J][4*k].center.y		=	tree[j][k].center.y-shift;
				tree[J][4*k+1].center.y	=	tree[j][k].center.y-shift;
				tree[J][4*k+2].center.y	=	tree[j][k].center.y+shift;
				tree[J][4*k+3].center.y	=	tree[j][k].center.y+shift;
			}
		}
	}

  void assignLeafChargeLocations() {
		for (size_t k = 0; k < nBoxesPerLevel[nLevels]; k++) {
      int startIndex = gridPoints.size();
      for (size_t i = 0; i < nParticlesInLeaf; i++) {
        tree[nLevels][k].chargeLocations.push_back(startIndex+i);
      }
      shift_Nodes(boxRadius[nLevels], tree[nLevels][k].center, gridPoints);
    }
    K->particles_X = gridPoints;//object of base class FMM_Matrix
    K->particles_Y = gridPoints;
  }

	void assignNonLeafChargeLocations() {
	  for (int j=nLevels-1; j>=2; --j) {
	    for (size_t k = 0; k < nBoxesPerLevel[j]; k++) {
	      for (int c=0; c<4; ++c) {
	        // for (size_t i = 0; i < tree[j+1][4*k+c].chargeLocations.size(); i++) {
	          tree[j][k].chargeLocations.insert(tree[j][k].chargeLocations.end(), tree[j+1][4*k+c].chargeLocations.begin(), tree[j+1][4*k+c].chargeLocations.end());
	        // }
	      }
	    }
	  }
	}

	void getNodes() {
		for (int j=2; j<=nLevels; j++) {
			getNodes_incoming_level(j);
		}
	}

	void getParticlesFromChildren_incoming(int j, int k, std::vector<int>& searchNodes) {
		if (j==nLevels) {
			searchNodes.insert(searchNodes.end(), tree[j][k].chargeLocations.begin(), tree[j][k].chargeLocations.end());
		}
		else {
			int J = j+1;
			for (int c = 0; c < 4; c++) {
				searchNodes.insert(searchNodes.end(), tree[J][4*k+c].incoming_checkPoints.begin(), tree[J][4*k+c].incoming_checkPoints.end());
			}
		}
	}

	void BebendrofPivots(std::vector<int>& inputIndices, std::vector<int>& searchNodes) {
		// if (j==nLevels) {
		// 	searchNodes = inputIndices;
		// }
		// else {
			if (nChebNodes1D*nChebNodes1D < inputIndices.size()) {
				// std::cout << "j: " << j << "	k: " << k << "	inputIndices.size(): " << inputIndices.size() << std::endl;
				geometricPivots(inputIndices, searchNodes);
			}
			else {
				searchNodes = inputIndices;
			}
		// }
	}

	void geometricPivots(std::vector<int>& inputIndices, std::vector<int>& searchNodes) {
		pts2D centroid;
		centroid.x = 0.0;
		centroid.y = 0.0;
		int NumPoints = inputIndices.size();
		for (size_t i = 0; i < NumPoints; i++) {
			pts2D temp = gridPoints[inputIndices[i]];
			centroid.x += temp.x;
			centroid.y += temp.y;
		}
		centroid.x /= NumPoints;
		centroid.y /= NumPoints;
		Eigen::MatrixXd U(2,NumPoints);
		for (size_t i = 0; i < NumPoints; i++) {
			pts2D temp = gridPoints[inputIndices[i]];
			U(0,i) = temp.x-centroid.x;
			U(1,i) = temp.y-centroid.y;
		}
		Eigen::MatrixXd cov = U*U.transpose();
		Eigen::EigenSolver<Eigen::MatrixXd> es(cov);
		// std::cout << "es.eigenvectors().real(): " << es.eigenvectors().real() << std::endl;
		Eigen::MatrixXd eigenVectors = es.eigenvectors().real();//size 2x2
		double angle = atan(eigenVectors(1,0)/eigenVectors(0,0));//x axis
		double angleRotate = -angle;
		Eigen::MatrixXd R(2,2);
		R(0,0) = cos(angleRotate);
		R(0,1) = -sin(angleRotate);
		R(1,0) = sin(angleRotate);
		R(1,1) = cos(angleRotate);
		Eigen::MatrixXd points(2,NumPoints);
		for (size_t i = 0; i < NumPoints; i++) {
			pts2D temp = gridPoints[inputIndices[i]];
			points(0,i) = temp.x;
			points(1,i) = temp.y;
		}
		// std::cout << "here" << std::endl;
		Eigen::MatrixXd rotatedNodes = R*points;
		pts2D radius;
		radius.x = rotatedNodes.row(0).maxCoeff() - rotatedNodes.row(0).minCoeff();
		radius.y = rotatedNodes.row(1).maxCoeff() - rotatedNodes.row(1).minCoeff();
		Eigen::Vector2d centroidVec(2);
		centroidVec(0) = centroid.x;
		centroidVec(1) = centroid.y;
		Eigen::Vector2d rotatedCenterVec = R*centroidVec;
		pts2D rotatedCenter;
		rotatedCenter.x = rotatedCenterVec(0);
		rotatedCenter.y = rotatedCenterVec(1);
		// std::cout << "sideX: " << radius.x << "	sideY: " << radius.y << std::endl;
		// std::cout << "rotatedCenterVec.X: " << rotatedCenterVec(0) << "	rotatedCenterVec.Y: " << rotatedCenterVec(1) << std::endl;
		std::vector<pts2D> shiftedChebNodes;
		shift_Cheb_Nodes(radius, rotatedCenter, shiftedChebNodes);
		// std::cout << "here" << std::endl;
		std::vector<int> pivotIndices;
		getPivots(shiftedChebNodes, rotatedNodes, pivotIndices);
		// std::cout << "shiftedChebNodes.size(): " << shiftedChebNodes.size() << std::endl;
		// std::cout << "pivotIndices.size(): " << pivotIndices.size() << std::endl;
		// for (size_t i = 0; i < pivotIndices.size(); i++) {
		// 	std::cout << pivotIndices[i] << std::endl;
		// }
		for (size_t i = 0; i < pivotIndices.size(); i++) {
			searchNodes.push_back(inputIndices[pivotIndices[i]]);
		}
		// std::cout << "searchNodes" << std::endl;
		// for (size_t i = 0; i < searchNodes.size(); i++) {
		// 	std::cout << searchNodes[i] << std::endl;
		// }
		// std::cout << "DONE" << std::endl;
		// exit(0);
	}

	void getPivots(std::vector<pts2D>& chebNodes, Eigen::MatrixXd& points, std::vector<int>& pivotIndices) {
		int NumPivots = chebNodes.size();
		std::set<int> remaining_ind;
		for(int l = 0; l < NumPivots; l++) {
				remaining_ind.insert(l);
		}
		// std::cout << "NumPivots: " << NumPivots << std::endl;
		for (size_t i = 0; i < NumPivots; i++) {
			// std::cout << "i: " << i << std::endl;
			int index = getClosestNode(chebNodes[i], remaining_ind, points); //find closest to chebNodes[i] and within the remaining indices
			pivotIndices.push_back(index);
			remaining_ind.erase(index);
		}
	}

	int getClosestNode(pts2D chebNode, const std::set<int>& allowed_indices, Eigen::MatrixXd& points) {
		std::set<int>::iterator it;
		double min   = 0.0;
		int index = *allowed_indices.begin();
		// std::cout << "here" << std::endl;
		// std::cout << "points.size(): " << points.rows() << ", " << points.cols() << std::endl;
		for(it = allowed_indices.begin(); it != allowed_indices.end(); it++) {
			// std::cout << "*it: " << *it << std::endl;
			pts2D point;
			point.x = points(0,*it);
			point.y = points(1,*it);
			double dist = squaredDistance(chebNode, point);
			if (min < dist) {
				index = *it;
				min = dist;
			}
		}
		return index;
	}

	double squaredDistance(pts2D p1, pts2D p2) {
		return((p1.x-p2.x)*(p1.x-p2.x) + (p1.y-p2.y)*(p1.y-p2.y));
	}

	void getNodes_incoming_box(int j, int k, int& ComputedRank) {
		std::vector<int> boxA_Nodes;
		BebendrofPivots(tree[j][k].chargeLocations, boxA_Nodes);
		// std::cout << "tree[j][k].chargeLocations.s: " << tree[j][k].chargeLocations.size() << "	boxA_Nodes.s: " << boxA_Nodes.size() << std::endl;
		// std::cout << "boxA_Nodes.size(): " << boxA_Nodes.size() << std::endl;
		// getParticlesFromChildren_incoming(j, k, boxA_Nodes);
		std::vector<int> IL_Nodes_temp, IL_Nodes;//indices
		for(int in=0; in<16; ++in) {
			if(tree[j][k].innerNumbers[in] != -1) {
				int kIL = tree[j][k].innerNumbers[in];
				std::vector<int> chargeLocations;
				// getParticlesFromChildren_incoming(j, kIL, chargeLocations);
				IL_Nodes_temp.insert(IL_Nodes_temp.end(), tree[j][kIL].chargeLocations.begin(), tree[j][kIL].chargeLocations.end());
			}
		}
		for(int on=0; on<24; ++on) {
			if(tree[j][k].outerNumbers[on] != -1) {
				int kIL = tree[j][k].outerNumbers[on];
				std::vector<int> chargeLocations;
				// getParticlesFromChildren_incoming(j, kIL, chargeLocations);
				IL_Nodes_temp.insert(IL_Nodes_temp.end(), tree[j][kIL].chargeLocations.begin(), tree[j][kIL].chargeLocations.end());
			}
		}
		if (j > 2) {
			IL_Nodes_temp.insert(IL_Nodes_temp.end(), tree[j-1][k/4].incoming_chargePoints.begin(), tree[j-1][k/4].incoming_chargePoints.end());
		}
		BebendrofPivots(IL_Nodes_temp, IL_Nodes);
		// std::cout << "IL_Nodes_temp.s: " << IL_Nodes_temp.size() << "	IL_Nodes.s: " << IL_Nodes.size() << std::endl;

		std::vector<int> row_bases, col_bases;
		Eigen::MatrixXd Ac, Ar;
		LowRank* LR		=	new LowRank(K, TOL_POW, boxA_Nodes, IL_Nodes);
		LR->ACA_only_nodes(row_bases, col_bases, ComputedRank, Ac, Ar, tree[j][k].L, tree[j][k].R);
		// LR->ACA_only_nodes(row_bases, col_bases, ComputedRank, Ac, Ar);
		if(ComputedRank > 0) {
			for (int r = 0; r < row_bases.size(); r++) {
				tree[j][k].incoming_checkPoints.push_back(boxA_Nodes[row_bases[r]]);
			}
			for (int c = 0; c < col_bases.size(); c++) {
				tree[j][k].incoming_chargePoints.push_back(IL_Nodes[col_bases[c]]);
			}
			// Eigen::MatrixXd D = K->getMatrix(tree[j][k].incoming_checkPoints, tree[j][k].incoming_chargePoints);
			// Eigen::PartialPivLU<Eigen::MatrixXd> D_T = Eigen::PartialPivLU<Eigen::MatrixXd>(D.transpose());
			// Eigen::PartialPivLU<Eigen::MatrixXd> D_T_T = Eigen::PartialPivLU<Eigen::MatrixXd>(D);
			// tree[j][k].L2P = D_T.solve(Ac.transpose()).transpose();

			// Eigen::MatrixXd temp = tree[j][k].R.triangularView<Eigen::Upper>().solve<Eigen::OnTheRight>(Ac);
			// tree[j][k].L2P = tree[j][k].L.triangularView<Eigen::Lower>().solve<Eigen::OnTheRight>(temp);

			// std::cout << "tree[j][k].L2P.size(): " << tree[j][k].L2P.rows() << ", " << tree[j][k].L2P.cols() << std::endl;
			// Eigen::MatrixXd Error_LR_D = L*R - D;
			// std::cout << "E_LR_D: " << Error_LR_D.norm()/D.norm() << std::endl;
			// Eigen::MatrixXd Error_AcAr_LR = U*V - Ac*D_T_T.solve(Ar);//tree[j][k].L2P*Ar;//K->getMatrix(boxA_Nodes, IL_Nodes);
			// std::cout << "E: " << Error_AcAr_LR.norm() << std::endl;
			//
			// std::cout << "L: " << std::endl << L << std::endl;
			// std::cout << "R: " << std::endl << R << std::endl;
			// std::cout << "LR: " << std::endl << L*R << std::endl;
			// std::cout << "D: " << std::endl << D << std::endl;
			// std::cout << "LR-D: " << std::endl << L*R-D << std::endl;
			// exit(0);

		}
	}

	void getNodes_incoming_level(int j) {
		int ComputedRank;
		for (int k=0; k<nBoxesPerLevel[j]; ++k) {
			// std::cout << "j: " << j << "	k: " << k << std::endl;
			getNodes_incoming_box(j, k, ComputedRank);
		}
	}

	// int getAvgRank() {
	// 	int sum = 0;
	// 	int NumBoxes = 0;
	// 	#pragma omp parallel for
	// 	for (size_t j = 2; j <= nLevels; j++) {
	// 		NumBoxes += pow(4, j);
	// 		#pragma omp parallel for
	// 		for (size_t k = 0; k < nBoxesPerLevel[j]; k++) {
	// 			sum += tree[j][k].incoming_checkPoints.size();
	// 		}
	// 	}
	// 	return sum/NumBoxes;
	// }

	int getAvgRank() {
		int sum = 0;
		int NumBoxes = 0;
		// #pragma omp parallel for
		for (size_t j = 2; j <= nLevels; j++) {
			// NumBoxes += nBoxesPerLevel[j];
			// #pragma omp parallel for
			for (size_t k = 0; k < nBoxesPerLevel[j]; k++) {
				NumBoxes += 1;
				sum += tree[j][k].incoming_checkPoints.size();
			}
		}
		return sum/NumBoxes;
	}

	void assemble_M2L() {
		#pragma omp parallel for
		for (size_t j = 2; j <= nLevels; j++) {
			#pragma omp parallel for
			for (size_t k = 0; k < nBoxesPerLevel[j]; k++) {
				#pragma omp parallel for
				for(int in=0; in<16; ++in) {
					if(tree[j][k].innerNumbers[in] != -1) {
						int kIL = tree[j][k].innerNumbers[in];
						if (tree[j][k].M2L[kIL].size() == 0)
							tree[j][k].M2L[kIL] = K->getMatrix(tree[j][k].incoming_checkPoints, tree[j][kIL].incoming_checkPoints);
						if (tree[j][kIL].M2L[k].size() == 0)
							tree[j][kIL].M2L[k] = tree[j][k].M2L[kIL].transpose();
					}
				}
				#pragma omp parallel for
				for(int on=0; on<24; ++on) {
					if(tree[j][k].outerNumbers[on] != -1) {
						int kIL = tree[j][k].outerNumbers[on];
						if (tree[j][k].M2L[kIL].size() == 0)
							tree[j][k].M2L[kIL] = K->getMatrix(tree[j][k].incoming_checkPoints, tree[j][kIL].incoming_checkPoints);
						if (tree[j][kIL].M2L[k].size() == 0)
							tree[j][kIL].M2L[k] = tree[j][k].M2L[kIL].transpose();
					}
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

	// void evaluate_M2M() {
	// 	for (size_t j = nLevels; j >= 2; j--) {
	// 		#pragma omp parallel for
	// 		for (size_t k = 0; k < nBoxesPerLevel[j]; k++) {
	// 			if (j==nLevels) {
	// 				Eigen::VectorXd temp = K->getMatrix(tree[j][k].incoming_chargePoints, tree[j][k].chargeLocations)*tree[j][k].charges;
	// 				Eigen::VectorXd temp2 = tree[j][k].R.transpose().triangularView<Eigen::Lower>().solve(temp);
	// 				tree[j][k].outgoing_charges = tree[j][k].L.transpose().triangularView<Eigen::Upper>().solve(temp2);
	// 			}
	// 			else {
	// 				tree[j][k].outgoing_charges = Eigen::VectorXd::Zero(tree[j][k].incoming_checkPoints.size());
	// 				for (size_t c = 0; c < 4; c++) {
	// 					Eigen::VectorXd temp = K->getMatrix(tree[j][k].incoming_chargePoints, tree[j+1][4*k+c].incoming_checkPoints)*tree[j+1][4*k+c].outgoing_charges;
	// 					Eigen::VectorXd temp2 = tree[j][k].R.transpose().triangularView<Eigen::Lower>().solve(temp);
	// 					tree[j][k].outgoing_charges += tree[j][k].L.transpose().triangularView<Eigen::Upper>().solve(temp2);
	// 				}
	// 			}
	// 		}
	// 	}
	// }

	void evaluate_M2L() {
		#pragma omp parallel for
		for (int j=2; j<=nLevels; ++j) {
			#pragma omp parallel for
			for (int k=0; k<nBoxesPerLevel[j]; ++k) {//BoxA
				tree[j][k].incoming_potential	=	Eigen::VectorXd::Zero(tree[j][k].incoming_checkPoints.size());
				#pragma omp parallel for
				for(int in=0; in<16; ++in) {
					if(tree[j][k].innerNumbers[in] != -1) {
						int kIL = tree[j][k].innerNumbers[in];
						tree[j][k].incoming_potential += tree[j][k].M2L[kIL]*tree[j][kIL].outgoing_charges;
					}
				}
				#pragma omp parallel for
				for(int on=0; on<24; ++on) {
					if(tree[j][k].outerNumbers[on] != -1) {
						int kIL = tree[j][k].outerNumbers[on];
						tree[j][k].incoming_potential += tree[j][k].M2L[kIL]*tree[j][kIL].outgoing_charges;
					}
				}
			}
		}
	}

	// void evaluate_L2L() {
	// 	for (size_t j = 2; j <= nLevels; j++) {
	// 		#pragma omp parallel for
	// 		for (size_t k = 0; k < nBoxesPerLevel[j]; k++) {
	// 			if (j != nLevels) {
	// 				for (size_t c = 0; c < 4; c++) {
	// 					Eigen::VectorXd temp2 = tree[j][k].L.triangularView<Eigen::Lower>().solve(tree[j][k].incoming_potential);
	// 					Eigen::VectorXd temp = tree[j][k].R.triangularView<Eigen::Upper>().solve(temp2);
	// 					tree[j+1][4*k+c].incoming_potential += K->getMatrix(tree[j+1][4*k+c].incoming_checkPoints, tree[j][k].incoming_chargePoints)*temp;
	// 				}
	// 			}
	// 			else {
	// 				Eigen::VectorXd temp2 = tree[j][k].L.triangularView<Eigen::Lower>().solve(tree[j][k].incoming_potential);
	// 				Eigen::VectorXd temp = tree[j][k].R.triangularView<Eigen::Upper>().solve(temp2);
	// 				tree[j][k].potential = K->getMatrix(tree[j][k].chargeLocations, tree[j][k].incoming_chargePoints)*temp;
	// 			}
	// 		}
	// 	}
	// }

	void assemble_M2M() {
		#pragma omp parallel for
		for (size_t j = 2; j <= nLevels; j++) {
			#pragma omp parallel for
			for (size_t k = 0; k < nBoxesPerLevel[j]; k++) {
				if (j==nLevels) {
					Eigen::MatrixXd temp = K->getMatrix(tree[j][k].incoming_chargePoints, tree[j][k].chargeLocations);
					Eigen::MatrixXd temp2 = tree[j][k].R.transpose().triangularView<Eigen::Lower>().solve(temp);
					tree[j][k].M2M[0] = tree[j][k].L.transpose().triangularView<Eigen::Upper>().solve(temp2);
				}
				else {
					for (size_t c = 0; c < 4; c++) {
						Eigen::MatrixXd temp = K->getMatrix(tree[j][k].incoming_chargePoints, tree[j+1][4*k+c].incoming_checkPoints);
						Eigen::MatrixXd temp2 = tree[j][k].R.transpose().triangularView<Eigen::Lower>().solve(temp);
						tree[j][k].M2M[c] = tree[j][k].L.transpose().triangularView<Eigen::Upper>().solve(temp2);
					}
				}
			}
		}
	}

	void evaluate_M2M() {
		for (size_t j = nLevels; j >= 2; j--) {
			#pragma omp parallel for
			for (size_t k = 0; k < nBoxesPerLevel[j]; k++) {
				if (j==nLevels) {
					tree[j][k].outgoing_charges = tree[j][k].M2M[0]*tree[j][k].charges;
				}
				else {
					tree[j][k].outgoing_charges = Eigen::VectorXd::Zero(tree[j][k].incoming_checkPoints.size());
					for (size_t c = 0; c < 4; c++) {
						tree[j][k].outgoing_charges += tree[j][k].M2M[c]*tree[j+1][4*k+c].outgoing_charges;
					}
				}
			}
		}
	}

	void evaluate_L2L() {
		for (size_t j = 2; j <= nLevels; j++) {
			#pragma omp parallel for
			for (size_t k = 0; k < nBoxesPerLevel[j]; k++) {
				if (j != nLevels) {
					for (size_t c = 0; c < 4; c++) {
						tree[j+1][4*k+c].incoming_potential += tree[j][k].M2M[c].transpose()*tree[j][k].incoming_potential;
					}
				}
				else {
					tree[j][k].potential = tree[j][k].M2M[0].transpose()*tree[j][k].incoming_potential;
				}
			}
		}
	}

	void evaluate_NearField() { // evaluating at chargeLocations
		#pragma omp parallel for
		for (size_t k = 0; k < nBoxesPerLevel[nLevels]; k++) {
			#pragma omp parallel for
			for (size_t n = 0; n < 8; n++) {
				int nn = tree[nLevels][k].neighborNumbers[n];
				if(nn != -1) {
					Eigen::MatrixXd R = K->getMatrix(tree[nLevels][k].chargeLocations, tree[nLevels][nn].chargeLocations);
					tree[nLevels][k].potential += R*tree[nLevels][nn].charges;
				}
			}
			Eigen::MatrixXd R = K->getMatrix(tree[nLevels][k].chargeLocations, tree[nLevels][k].chargeLocations);
			tree[nLevels][k].potential += R*tree[nLevels][k].charges; //self Interaction
		}
	}

	void collectPotential(Eigen::VectorXd &potential) {
		potential = Eigen::VectorXd::Zero(N);
		int start = 0;
		for (size_t k = 0; k < nBoxesPerLevel[nLevels]; k++) {
			potential.segment(start, tree[nLevels][k].potential.size()) = tree[nLevels][k].potential;
			start += tree[nLevels][k].potential.size();
		}
	}

	void findMemory(double &sum) {
		sum = 0.0;
		for (size_t j = 2; j <= nLevels; j++) {
			for (size_t k = 0; k < nBoxesPerLevel[j]; k++) {
				for(int in=0; in<16; ++in) {
					if(tree[j][k].innerNumbers[in] != -1) {
						int kIL = tree[j][k].innerNumbers[in];
						sum += tree[j][k].M2L[kIL].size();
					}
				}
				for(int on=0; on<24; ++on) {
					if(tree[j][k].outerNumbers[on] != -1) {
						int kIL = tree[j][k].outerNumbers[on];
						sum += tree[j][k].M2L[kIL].size();
					}
				}
			}
		}

		for (size_t j = 2; j <= nLevels; j++) {
			// #pragma omp parallel for
			for (size_t k = 0; k < nBoxesPerLevel[j]; k++) {
				if (j != nLevels) {
					for (size_t c = 0; c < 4; c++) {
						sum += tree[j][k].M2M[c].size();
					}
				}
				else {
					sum += tree[j][k].M2M[0].size();
				}
			}
		}
		for (size_t k = 0; k < nBoxesPerLevel[nLevels]; k++) {
			sum += tree[nLevels][k].chargeLocations.size()*tree[nLevels][k].chargeLocations.size(); //self
			for (size_t n = 0; n < 8; n++) {
				int nn = tree[nLevels][k].neighborNumbers[n];
				if(nn != -1) {
					sum += tree[nLevels][k].chargeLocations.size()*tree[nLevels][nn].chargeLocations.size();
				}
			}
		}
	}

	double compute_error(int nBox) { // evaluating at chargeLocations
		Eigen::VectorXd true_potential = Eigen::VectorXd::Zero(tree[nLevels][nBox].chargeLocations.size());
		for (size_t k = 0; k < nBoxesPerLevel[nLevels]; k++) {
			Eigen::MatrixXd R = K->getMatrix(tree[nLevels][nBox].chargeLocations, tree[nLevels][k].chargeLocations);
			true_potential += R*tree[nLevels][k].charges;
		}
		Eigen::VectorXd errVec = true_potential - tree[nLevels][nBox].potential;
		double error = errVec.norm()/true_potential.norm();
		return error;
	}
};

#endif
