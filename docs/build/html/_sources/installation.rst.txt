.. _installation:

*************************
Installation and Building
*************************

Downloading the Source
-----------------------

:math:`\texttt{NNCAlib}` is distributed using the git version control system, and is hosted on Github. The repository can be cloned using::

    git clone https://github.com/SAFRAN-LAB/NNCA.git --recursive

The ``--recursive`` flag is argument to ensure that all submodules are also checked out.

Dependencies
-------------

- Eigen Linear Algebra Library (get it `here <https://eigen.tuxfamily.org/index.php?title=Main_Page>`_)
- (optional) An OpenMP enabled compiler (e.g. gcc4.2 or above) is required to use shared-memory parallelism.

**NOTE**: On MacOS, the default compiler is `clang` which doesn't have OpenMP support. You will have to use g++ to make use of the speedups from OpenMP::

    user@computer HODLR3D$ brew install g++-13

Installation
-------------

Set the environment variable ``EIGEN_PATH`` to the location of your Eigen installation. This is needed by the Make file.::

    user@computer HODLR3D$ export EIGEN_PATH=path/to/eigen/

.. _Testing:

Testing
-------

Now, we need to ensure that all the functions of the libraries function as intended. For this purpose, we will be running the script ``NNCA2D/testFMM2D.cpp``.
Run ``NNCA2D/Makefile2D.mk`` to get your executable.
To check this on your computer, run the following lines::

    user@computer NNCA$ cd NNCA2D
    user@computer NNCA2D$ make -f Makefile2D.mk clean
    user@computer NNCA2D$ make -f Makefile2D.mk
    user@computer NNCA2D$ ./testFMM2D 5 6 1 5 1
