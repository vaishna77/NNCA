This library was created by Vaishnavi Gujjula, during the course of her Ph.D.

This library does matrix vector product for FMM-matrices, a sub-class of H2 matrices with well-separated admissibility condition.

The low-rank approximations of matrix sub-blocks is constructed using a new Nested Cross Approximation[[1](https://arxiv.org/abs/2203.14832)].

The particles are distributed at the leaf nodes and are distributed uniformly in square [-L,L]^2.

The matrix entries are to be defined in function "getMatrixEntry(i,j)" of the kernel.hpp file.

The vector to be applied to the matrix is to be defined in VectorXd "b" of the testFMM2D.cpp.
