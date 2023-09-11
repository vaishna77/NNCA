#ifndef myHeaders_HPP
#define myHeaders_HPP

#include <iostream>
#include <chrono>
#include <ctime>
#include <set>
#include <string>
#include <vector>
#include <random>
#include <cmath>
#include <complex>
#include <iterator>
#include <math.h>
#include <iomanip>
#include <fstream>
#include <list>
#include <algorithm>
#include <numeric>
#include <stack>

#define _USE_MATH_DEFINES
#include <Eigen/Dense>
#include <Eigen/Sparse>

using std::string;
using namespace Eigen;

const double PI = 4.0 * atan(1);

const int nThreads = 4;
/////////////////////////////////////////////////////////////
// matrix parameters
#ifdef USE_DIM2
const int NDIM = 2;
#elif USE_DIM4
const int NDIM = 4;
#endif
const int Nmax = 500;
// const int numPoints = 75; // Along 1D // use 8,9,10,11 values to reproduce the results in 4D of the article
int numPoints; // Along 1D // use 8,9,10,11 values to reproduce the results in 4D of the article
// The admissibility is based on the max norm of the center
const int INTERACTION_TYPE_ALLOWED = 0; // This represents d'
// 2D
// d' = 0 -> Vertex, HODLR2D
// d' = 1 -> HODLR in 2D
// 3D
// d' = 0 -> Vertex, HODLR3D
// d' = 1 -> edge
// d' = 2 -> face
const double eps_ACA = pow(10,-7);
// const int N = pow(numPoints,NDIM);
int N;
const int SYS_SIZE = 2*Nmax;

/////////////////////////////////////////////////////////////

int mod(int a, int b){
    return ((a % b + b) % b);
}

// #ifdef USE_FLOAT
//     using dtype = float;
//     using dtype_base = float;
//     using Mat = Eigen::MatrixXf;
//     using Vec = Eigen::VectorXf;
//     #define Calc_dist(x1, y1, x2, y2) float(sqrt(pow(x1 - x2, 2) + pow(y1 - y2, 2)))
// #endif

// #ifdef USE_DOUBLE
    using dtype = double;
    using dtype_base = double;
    using Mat = Eigen::MatrixXd;
    using Vec = Eigen::VectorXd;
    #define abs_(x) ((x < 0) ? (-x) : (x))
    #define conj_(x) (x)
    #define Calc_dist(x1, y1, x2, y2) double(sqrt(pow(x1 - x2, 2) + pow(y1 - y2, 2)))
// #endif

// #ifdef USE_COMPLEX32
//     using dtype = std::complex<float>;
//     using dtype_base = float;
//     using Mat = Eigen::MatrixXcf;
//     using Vec = Eigen::VectorXcf;
//     const std::complex<float> I(0.0, 1.0);
// #endif

// #ifdef USE_COMPLEX64
//     using dtype = std::complex<double>;
//     using dtype_base = double;
//     using Mat = Eigen::MatrixXcd;
//     using Vec = Eigen::VectorXcd;
//     #define abs_(x) (std::abs((x)))
//     #define conj_(x) (std::conj((x)))
//     const std::complex<double> I(0.0, 1.0);
//     #define Calc_dist(x1, y1, x2, y2) float(sqrt(pow(x1 - x2, 2) + pow(y1 - y2, 2)))
// #endif


// Cheb nodes
Eigen::VectorXd cheb_nodes(double a, double b, int n)
{
    Eigen::VectorXd X(n);
    double l, l1, param;
    l = 0.5 * (a + b);
    l1 = 0.5 * (b - a);
    for (int k = 0; k < n; k++)
    {
        param = (double)(k + 0.5) / n;
        X(k) = l - l1 * cos(param * 3.1412);
    }
    return X;
}

// Uniform nodes
Eigen::VectorXd uniform_nodes(double a, double b, int n)
{
    Eigen::VectorXd X(n);
    double l = (b-a);
    for (int k = 0; k < n; k++)
    {
        X(k) = a + l/(double) (n-1) * k;
    }
    return X;
}

std::vector<double> random_nodes(double a, double b, int n)
{
    std::vector<double> X(n);
    Eigen::VectorXd Points = Eigen::VectorXd::Random(n);
    double l, l1;
    l = 0.5 * (a + b);
    l1 = 0.5 * (b - a);
    for (int k = 0; k < n; k++)
    {
        //param = (double)(k + 0.5) / n;
        X[k] = l - l1 * Points(k);
    }
    return X;
}

std::string inline getTimeStamp()
{
    std::string s;
    s = "\n" + string(__DATE__) + "," + string(__TIME__) + "\n";
    return s;
}

namespace Vec_ops
{
    dtype_base inline relative_error(Vec &X_ori, Vec &X_comp)
    {
        dtype_base result = 0.0;
        result = (X_ori - X_comp).norm() / X_ori.norm();
        return result;
    }
    dtype inline dot_product(const Vec &x, const Vec &y)
    {
        dtype tmp = dtype(0.0);
        for (int i = 0; i < x.size(); i++)
            tmp += x(i) * y(i);
        return tmp;
    }
}

#ifdef USE_CEREAL
#include <cereal/archives/portable_binary.hpp>
#include <cereal/types/vector.hpp>
namespace storedata
{
    Eigen::VectorXd inline load_vec(std::string fname)
    {
        Eigen::VectorXd X;
        std::vector<dtype_base> K;
        std::ifstream infile(fname, std::ios::binary);
        cereal::PortableBinaryInputArchive iarchive(infile); // Create an input archive
        iarchive(K);                                         // Read the data from the archive
        X = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(K.data(), K.size());
        return X;
    }

    void inline save_vec(std::string fname, Eigen::VectorXd X)
    {
        std::vector<dtype_base> K(X.data(), X.data() + X.size());
        std::ofstream outfile(fname, std::ios::binary);
        cereal::PortableBinaryOutputArchive oarchive(outfile); // Create an output archive
        oarchive(CEREAL_NVP(K));                               // Write the data to the archive
    }
}
#endif

template <typename T>
void inline print_vec(std::vector<T> &vec_to_print)
{
    for (auto iter = vec_to_print.begin(); iter != vec_to_print.end(); iter++)
    {
        if (iter != vec_to_print.begin())
            std::cout << ", ";
        std::cout << *iter;
    }
}

#endif
