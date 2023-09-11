#ifndef SVM_HPP
#define SVM_HPP
// const int NDIM = 5;

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <algorithm>
#include <assert.h>
#include <cmath>
#include <Eigen/Dense>
#include "AFMMnD/AFMMnD.hpp"
// #include "HODLRnD/myHeaders.hpp"
// #include "HODLRnD/points_dt.hpp"
#include "myHeaders.hpp"
#include "points_dt.hpp"
std::ofstream testfile("matlab/test_label.txt");

using namespace std;
class svm{
    public:
    std::vector<ptsnD>* x = new std::vector<ptsnD>;
    vector<int>* ind = new std::vector<int>;
    double C, eta, limit;  // The hyper-params you need to tune
    double accuracy;  //overall accuracy
    double accuracy_c1, accuracy_c2;
    double b = 0.0;
    size_t correct_c1 = 0, correct_c2 = 0;

    // To store the information of support vectors
    std::vector<ptsnD>* xs = new std::vector<ptsnD>;
    std::vector<int>* ys = new std::vector<int>;
    std::vector<double>* alpha_s = new std::vector<double>;
    std::vector<ptsnD>* xs_in = new std::vector<ptsnD>;
    std::vector<int>* ys_in = new std::vector<int>;
    std::vector<double>* alpha_s_in = new std::vector<double>;
    /////////////////////////////////////////////

    const double eps = 1e-12;
    Eigen::VectorXd term1;
    double term2, term3, delta, beta = 1.0;
    double error;

    Eigen::VectorXd X, Y; // HODLRdD purpose
    Eigen::VectorXd y;  // This vector store the label of the datas
    Eigen::VectorXd alpha; // This is vector because we need it to test sample
    int Qchoice;
    int nLevels; //of AFMMnD
    bool flag;

    size_t Ns=0, Ns_in=0;
    size_t N;

    int iter=0;   // To keep the iteration count of GD

    // Empty constructor
    svm(){
        this->C = 0.5;
        this->eta = 0.05;
        this->limit = 0.2;
        this->Qchoice = 0;
        this->nLevels = 4;
    }
    // Another constructor
    svm(double C, double eta, double limit, int Qchoice, int nLevels){
        this->C = C;
        this->eta = eta;
        this->limit = limit;
        this->Qchoice = Qchoice;
        this->nLevels = nLevels;
    }
    ~svm(){}

    // Mate'rn kernel
    #ifdef USE_Matern
    double Ker_Function(const ptsnD& a, const ptsnD& b){
        double temp = 0.0;
        for(size_t i=0; i<NDIM; ++i){
            temp += (a.x[i] - b.x[i]) * (a.x[i] - b.x[i]);
        }
        return exp(-1.0 * sqrt(temp));
    }
    double Ker_Function(double r){
        return exp(-1.0*r);
    }
    #endif

    // Mate'rn kernel
    #ifdef USE_Gaussian
    double Ker_Function(const ptsnD& a, const ptsnD& b){
        double temp = 0.0;
        for(size_t i=0; i<NDIM; ++i){
            temp += (a.x[i] - b.x[i]) * (a.x[i] - b.x[i]);
        }
        return exp(-1.0 * temp);
    }
    double Ker_Function(double r){
        return exp(-1.0*r*r);
    }
    #endif

    double euclidean_distance(ptsnD& a, ptsnD& b){
        double dist = 0.0;
        for (size_t i = 0; i < NDIM; i++){
            dist += (a.x[i] - b.x[i]) * (a.x[i] - b.x[i]);
        }
        dist = sqrt(dist);
        return dist;
    }
    double getMatrixEntry(int i, int j)
    {
        double r;
        double dist = euclidean_distance(x->at(i), x->at(j));
        r = Ker_Function(dist);
        return r;
    }
    Eigen::MatrixXd getMatrix(std::vector<int> &sources, std::vector<int> &targets){
        int n_rows = sources.size();
        int n_cols = targets.size();
        Eigen::MatrixXd mat(n_rows, n_cols);
        for (int j = 0; j < n_rows; ++j){
            for (int k = 0; k < n_cols; ++k){
                mat(j, k) = getMatrixEntry(sources[j], targets[k]);
            }
        }
        return mat;
    }

    // Training parameters
    void set_training_params(const std::vector<ptsnD>& class1_data, const std::vector<ptsnD>& class2_data){
        y.resize(class1_data.size()+class2_data.size());

        // Set class 1 data
        for(auto& itr : class1_data){
            x->push_back(itr);
        }

        // Set labels of class 1 (label = 1)
        y.segment(0,class1_data.size()) = Eigen::VectorXd::Ones(class1_data.size());

        // Set class 2 data, label = 1
        for(auto& itr : class2_data){
            x->push_back(itr);
        }
        // Set labels class 2, label = -1
        y.segment(class1_data.size(),class2_data.size()) = -1.0*Eigen::VectorXd::Ones(class2_data.size());

        // Set Lagrange Multiplier and Parameters
        N = x->size();
				// std::cout << "set train N: " << N << std::endl;

        //std::vector<double> alpha = std::vector<double>(N, 0.0);   // This is vector because we need it to test sample

        // Set index set
        for(size_t i=0; i<N; ++i){
            ind->push_back(i);
        }
    }

    // Assemble the HODLR4D matrix
    void assemble_hodlr(){
        X.resize(NDIM);
        Y.resize(NDIM);
        for (int i = 0; i < NDIM; ++i){
            X(i) = -1;
            Y(i) = 1;
        }
    }

    void train(){
        // kernel_5d_test *ker = new kernel_5d_test();
        // ker->get_points(x);
        // HODLRdD_matrix Kmat = HODLRdD_matrix(ker, x, X, Y);
        // Kmat.Assemble_matrix_operators();
				std::vector<AptsnD> points;
				// std::cout << "N: " << N << std::endl;
				for (size_t i = 0; i < N; i++) {//N=pow(n,ndim)
					ptsnD temp = x->at(i);
					// std::cout << "i: " << i << std::endl;
					AptsnD Atemp;
					for (size_t j = 0; j < NDIM; j++) {
						Atemp.push_back(temp.x[j]);
					}
					points.push_back(Atemp);
				}

				AFMMnD* afmmnd = new AFMMnD(NDIM, N, 6, 1, 6, Qchoice, nLevels, points);
        #ifdef USE_AFMMnD
        afmmnd->assembleAFMMnD();
        #endif

        alpha = Eigen::VectorXd::Zero(N);
        do{
						// std::cout << "here" << std::endl;
						flag = false;
            error = 0.0;
            Eigen::VectorXd V = alpha.cwiseProduct(y);
            #ifdef USE_AFMMnD
						term1 = afmmnd->matVecProduct(V);//mat-vec using afmmnd
            #endif
						//////////////////direct mat-vec//////////////////////////////////
            #ifdef USE_directMatVec
						term1 = Eigen::VectorXd::Zero(V.size());
						for (size_t i = 0; i < V.size(); i++) {
							Eigen::VectorXd col = afmmnd->mykernel->getCol(i);
							term1 += V(i)*col;
						}
            #endif
						// term1 = Kmat * V;     //! MAT-VEC
						for(size_t i=0; i<N; ++i){
                // Term 1
                term1(i) *= y(i);

                // Term 2
                term2 = alpha.dot(y);
                term2 *= y(i);

                // Set delta (Partial derivative)
                delta = 1.0 - term1(i) - beta * term2;

                // Upadte alpha
                alpha(i) += eta * delta;
                if(alpha(i) < 0.0){
                    alpha(i) = 0.0;
                }
                else if(alpha(i) > C){
                    alpha(i) = C;
                }
                else if(fabs(delta) > limit){
                    flag = true;
                    error += fabs(delta) - limit;
                }
            }
						// std::cout << "here" << std::endl;
            // Update Beta
            term3 = 0.0;
            for(size_t i = 0; i < N; ++i){
                term3 += alpha(i) * y(i);
            }
            beta += term3 * term3 / 2.0;
            // if(error !=0)
            //     std::cout << "i: " << iter << "	err: " << error << std::endl;
            ++iter;

        }while(flag);
				// std::cout << "out of loop" << std::endl;

        // Information of support vectors
        for(size_t i = 0; i < N; ++i){
            if((eps < alpha(i)) && (alpha(i) < C - eps)){
                this->xs->push_back(x->at(i));
                this->ys->push_back(y(i));
                this->alpha_s->push_back(alpha(i));
                ++Ns;
            }
            else if(alpha(i) >= C - eps){
                this->xs_in->push_back(x->at(i));
                this->ys_in->push_back(y(i));
                this->alpha_s_in->push_back(alpha(i));
                ++Ns_in;
            }
        }
        if (Ns > 0) {
          // Update bias
          for(size_t i = 0; i < Ns; ++i){
              this->b += (double)this->ys->at(i);
              for(size_t j = 0; j < Ns; ++j){
                  double dist = euclidean_distance(this->xs->at(j), this->xs->at(i));
                  this->b -= this->alpha_s->at(j) * (double)this->ys->at(j) * this->Ker_Function(dist);
              }
              for(size_t j = 0; j < Ns_in; ++j){
                  double dist = euclidean_distance(this->xs_in->at(j), this->xs->at(i));
                  this->b -= this->alpha_s_in->at(j) * (double)this->ys_in->at(j) * this->Ker_Function(dist);
              }
          }
          this->b /= (double)Ns;
        }
    }

    // The decision function
    double fun(const ptsnD& x){
        double temp = 0.0;
        for(size_t i = 0; i < xs->size(); i++){
            temp += this->alpha_s->at(i) * this->ys->at(i) * this->Ker_Function(this->xs->at(i), x);
        }
        for(size_t i = 0; i < xs_in->size(); i++){
            temp += this->alpha_s_in->at(i) * this->ys_in->at(i) * this->Ker_Function(this->xs_in->at(i), x);
        }
        temp += this->b;
        return temp;
    }

    double dpendent_fun(const ptsnD& x){
        double fx = this->fun(x);
        double temp = (fx >= 0.0) ? 1.0 : -1.0;
        testfile << fx << " " << temp << std::endl;
        return temp;
    }

    // Now it's time to check the model
    void test(const std::vector<ptsnD>& class1_data, const std::vector<ptsnD>& class2_data){
        for(auto& itr : class1_data){
            if (this->dpendent_fun(itr) == 1.0){
                ++this->correct_c1;                         //! TODO MAT-VEC
            }
        }
        for(auto& itr : class2_data){
            if (this->dpendent_fun(itr) == -1.0){
                ++this->correct_c2;                         //! TODO MAT-VEC
            }
        }

        this->accuracy = (double)(this->correct_c1 + this->correct_c2) / (double)(class1_data.size() + class2_data.size());
        this->accuracy_c1 = (double)this->correct_c1 / (double)class1_data.size();
        this->accuracy_c2 = (double)this->correct_c2 / (double)class2_data.size();
    }

    // Print all details
    void print_details(){
        std::cout << "--------------------------------------------------" << std::endl;
        std::cout << "The number of iterations: " << this->iter << std::endl;
        std::cout << "The number of support vectors on margin: " << this->Ns << std::endl;
        std::cout << "The number of support vectors inside margin: " << this->Ns_in << std::endl;
        std::cout << "The accuracy of class1: " << this->accuracy_c1 << std::endl;
        std::cout << "The accuracy of class2: " << this->accuracy_c2 << std::endl;
        std::cout << "The over all accuracy of: " << this->accuracy << std::endl;
        std::cout << "Correct C1: " << this->correct_c1 << std::endl;
        std::cout << "Correct C2: " << this->correct_c2 << std::endl;
        std::cout << "--------------------------------------------------" << std::endl;

    }

};

#endif
