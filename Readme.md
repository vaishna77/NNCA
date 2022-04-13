This library does matrix vector product for FMM-matrices, a sub-class of H2 matrices with well-separated admissibility condition.

The low-rank approximations of matrix sub-blocks is constructed using a new Nested Cross Approximation.

The charges are distributed at the leaf nodes and are distributed uniformly in square [-L,L]^2.

The matrix entries are to be defined in a function getMatrixEntry(i,j) in the kernel.hpp file.

The vector to be applied to the matrix is to be defined in VectorXd 'b' in the testFMM2D.cpp
