#include <bits/stdc++.h>
#include <iostream>
#include <vector>
#include <Eigen/Dense>

#include <map>
typedef std::vector<double> AptsnD;

class kernel {
public:
  double a;
  std::vector<AptsnD> particles;

	kernel(std::vector<AptsnD>& particles) {
			this->particles = particles;
	}

	virtual double getMatrixEntry(const unsigned i, const unsigned j) {
		std::cout << "virtual getInteraction" << std::endl;
		return 0.0;
	}

	Eigen::VectorXd getRow(const int j, std::vector<int> col_indices) {
		int n_cols = col_indices.size();
		Eigen::VectorXd row(n_cols);
    // #pragma omp parallel for
    for(int k = 0; k < n_cols; k++) {
        row(k) = this->getMatrixEntry(j, col_indices[k]);
    }
    return row;
  }

  Eigen::VectorXd getCol(const int k, std::vector<int> row_indices) {
		int n_rows = row_indices.size();
    Eigen::VectorXd col(n_rows);
    // #pragma omp parallel for
    for (int j=0; j<n_rows; ++j) {
			col(j) = this->getMatrixEntry(row_indices[j], k);
    }
    return col;
  }

	Eigen::VectorXd getCol(const int k) {
		Eigen::VectorXd col(particles.size());
		// #pragma omp parallel for
		for (int j=0; j<particles.size(); ++j) {
			col(j) = this->getMatrixEntry(j, k);
		}

		return col;
	}

  Eigen::MatrixXd getMatrix(std::vector<int> row_indices, std::vector<int> col_indices) {
		int n_rows = row_indices.size();
		int n_cols = col_indices.size();
    Eigen::MatrixXd mat(n_rows, n_cols);
    // #pragma omp parallel for collapse(2)
		// #pragma omp parallel for
    for (int j=0; j < n_rows; ++j) {
        // #pragma omp parallel for
        for (int k=0; k < n_cols; ++k) {
            mat(j,k) = this->getMatrixEntry(row_indices[j], col_indices[k]);
        }
    }
    return mat;
  }

	Eigen::MatrixXd getMatrix(int row_start_index, int col_start_index, int row_end_index, int col_end_index) {
		Eigen::MatrixXd mat(row_end_index-row_start_index, col_end_index-col_start_index);
		// #pragma omp parallel for
		for (int j=row_start_index; j < row_end_index; ++j) {
				// #pragma omp parallel for
				for (int k=col_start_index; k < col_end_index; ++k) {
						mat(j,k) = this->getMatrixEntry(j, k);
				}
		}
		return mat;
	}

  ~kernel() {};
};

class userkernel: public kernel {
public:
	int Qchoice;
	double Kii;
	double h3;
	// double chargesFunction(const AptsnD r) {
	// 	double q = r.x; //user defined
	// 	return q;
	// };
	// #ifdef ONEOVERR
	// userkernel(std::vector<AptsnD>& particles_X, std::vector<AptsnD>& particles_Y): kernel(particles_X, particles_Y) {
	// };

	// #elif LOGR
	userkernel(std::vector<AptsnD> particles, int Qchoice): kernel(particles) {
		this->Qchoice = Qchoice;
	};

	double getr(AptsnD &r1, AptsnD &r2) {
		double r = 0.0;
		int Ndim = r1.size();
		for (size_t i = 0; i < Ndim; i++) {
			r += (r1[i]-r2[i])*(r1[i]-r2[i]);
		}
		r = sqrt(r);
		return r;
	}

	double Matern(const unsigned i, const unsigned j) {
		AptsnD r1 = particles[i];
		AptsnD r2 = particles[j];
		double R = getr(r1,r2);
		return exp(-R);
		// return 1;
	}

  double RBF_Gaussian(const unsigned i, const unsigned j) {
		AptsnD r1 = particles[i];
		AptsnD r2 = particles[j];
		double R = getr(r1,r2);
		return exp(-R*R);
		// return 1;
	}

	// RBF Logarithm
	double RBF_Logarithm(const unsigned i, const unsigned j) {
		AptsnD r1 = particles[i];
		AptsnD r2 = particles[j];
		double R = getr(r1,r2);
		double b = 1.0;
		return log(1.0+R/b);
	}


	double getMatrixEntry(const unsigned i, const unsigned j) {
    // if (i==j) {
    //   return 0.0;
    // }
    // else {
      if (Qchoice == 0)
  			return Matern(i,j);
      else if(Qchoice == 1)
        return RBF_Gaussian(i,j);
    }
	// }

	~userkernel() {};
};
