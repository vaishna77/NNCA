#ifndef __points_dt__
#define __points_dt__

#include <cmath>
#include "myHeaders.hpp"


// Struture that holds the individual location of points
struct ptsnD
{
    double x[NDIM];
    size_t id=-1; // User provides a unique identifier if not then kernel is defined by considering two points
};
// void print_point(ptsnD& a){
//     for(int i=0;i<NDIM;i++)
//         std::cout << a.x[i] <<" " << std::endl;
// }
// namespace nd_points{
// // Computes the infinity norm of the |x-y| R^d -> R
// double max_norm_distance(ptsnD& a, ptsnD& b){
//     double dist = abs(a.x[0] - b.x[0]);
//     for(int i=0; i<NDIM; i++){
//         double t = abs(a.x[i] - b.x[i]);
//         if(t>dist)
//             dist = t;
//     }
//     return dist;
// }
//
// // Computes the Euclidean norm or 2 norm of the difference vector
// double euclidean_distance(ptsnD& a, ptsnD& b)
// {
//     double dist = 0.0;
//     for (int i = 0; i < NDIM; i++){
//         dist += (a.x[i] - b.x[i]) * (a.x[i] - b.x[i]);
//     }
//     dist = sqrt(dist);
//     return dist;
// }
// }
#endif
