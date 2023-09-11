********
Tutorial
********


The project "NNCA" has the following sub-projects/directories:

  #. docs:: It holds the files related to this documentation
  #. NCA2DBebendorf:: It holds the implementation of the NCA developed by Bebendorf in C++ for 2D problems. Its performance is tested by evaluating a fast matrix-vector product using NCA.
  #. NCA2DZhao::  It holds the implementation of the NCA developed by Zhao in C++ for 2D problems. Its performance is tested by evaluating a fast matrix-vector product using NCA.
  #. NNCA2D::  It holds the implementation of the NNCA that is proposed in C++ for 2D problems. Its performance is tested by evaluating a fast matrix-vector product using NNCA.
  #. NNCA3D::  It holds the implementation of the NNCA that is proposed in C++ for 3D problems. Its performance is tested by evaluating a fast matrix-vector product using NNCA and by solving an integral equation solver.
  #. NNCAnD::  It holds the implementation of the NNCA that is proposed in C++ for nD problems. Its performance is tested by evaluating a fast matrix-vector product using NNCA.
  #. FSVM::  It holds the implementation of Support Vector Machine (SVM) using NNCA in C++ for nD problems.

Below is a description of the sub-projects listed above.

.. _settingParams:
Setting Parameters in make file
-------------------------------

#. Each of these sub-projects has a make file with extension ``.mk``. The following variable needs to be set by the user in the make file before running it:

NCA2DBebendorf, NCA2DZhao, NNCA2D, NNCA3D, NNCAnD, FSVM
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  - ``CC``: Set it to the path where the g++ compiler is located in your computer.

FSVM
^^^^
  - ``KERNEL``: Set it to ``-DUSE_Matern`` to use Matern kernel or ``-DUSE_Gaussian`` to use Gaussian kernel.
  - ``DIM``: Set it to the dimension of the problem. For example ``2`` for 2D, ``4``for 4D, etc.
  - ``MATVEC``: Set it to  ``-DUSE_directMatVec`` to implement normal SVM (where matrix-vector products are computed naively) or ``-DUSE_AFMMnD`` to implement fast SVM (where matrix-vector products are accelerated using NNCA).

Building and Executing
----------------------

#. Install the library by following the steps given in the :ref:`installation`
#. Each of these sub-projects has make file(s) with extension ``.mk``. Set the various parameters in the make file as described in :ref:`settingParams`.
#. Run the make file.
#. While running it make sure to give the required run-time inputs. To know details about the run-time inputs look at the main function in the file with name mentioned under ``SOURCES`` flag in the make file. Below is a description of the inputs required from the user for different sub-projects::

NCA2DBebendorf, NCA2DZhao, NNCA2D
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  a. ``nLevels`` : no. of levels to build in the hierarchical quad-tree
  b. ``nParticlesInLeafAlong1D`` : square root of the maximum no. of particles a leaf can have
  c. ``L`` : half side length of the square domain centered at origin (It is assumed the domain is a square centered at origin with half side length equal to L).
  d. ``TOL_POW`` : pow(10, -TOL_POW) is the tolerance the code uses to terminate ACA
  e. ``yesToUniformDistribution`` : set it to 1 to consider a uniform distribution of particles, 0 to consider Chebyshev distribution

NNCA3D - fast matrix-vector product (Makefile3D.mk)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  a. ``cubeRootN`` : no. of levels to build in the hierarchical quad-tree
  b. ``nParticlesInLeafAlong1D`` : square root of the maximum no. of particles a leaf can have
  c. ``L`` : half side length of the cube domain centered at origin (It is assumed the domain is a cube centered at origin with half side length equal to L).
  d. ``TOL_POW`` : pow(10, -TOL_POW) is the tolerance the code uses to terminate ACA
  e. ``Qchoice`` : It is used to choose the kernel function among multiple definitions defined in the kernel.hpp file. Look at the function ``getMatrixEntry`` in the ``kernel.hp`` file for the choices to be used for various definitions. For example, to use the Laplacian 3D kernel set ``Qchoice`` to 7.
  f. ``yesToUniformDistribution`` : set it to 1 to consider a uniform distribution of particles, 0 to consider Chebyshev distribution

NNCA3D - integral equation solver (Makefile3Dsolve.mk)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  a. ``cubeRootN`` : no. of levels to build in the hierarchical quad-tree
  b. ``nParticlesInLeafAlong1D`` : square root of the maximum no. of particles a leaf can have
  c. ``L`` : half side length of the cube domain centered at origin (It is assumed the domain is a cube centered at origin with half side length equal to L).
  d. ``TOL_POW`` : pow(10, -TOL_POW) is the tolerance the code uses to terminate ACA

NNCAnD
^^^^^^
  a. ``ndim`` : The dimension of the problem being solved
  b. ``n`` : It is used to set the number of particles (target points and charge points are assumed to be the same). The number of particles is equal to pow(n,ndim).
  c. ``nParticlesInLeafAlong1D`` : It is used to set the maximum no. of particles a leaf can have. the maximum no. of particles a leaf can have is equal to pow(nParticlesInLeafAlong1D,ndim).
  d. ``L`` : half side length of the hypercube domain centered at origin (It is assumed the domain is a hypercube centered at origin with half side length equal to L).
  e. ``TOL_POW`` : pow(10, -TOL_POW) is the tolerance the code uses to terminate ACA
  f. ``Qchoice`` : It is used to choose the kernel function among multiple definitions defined in the kernel.hpp file. Look at the function ``getMatrixEntry`` in the ``kernel.hp`` file for the choices to be used for various definitions. For example, to use the Laplacian 3D kernel set ``Qchoice`` to 7.

FSVM
^^^^
  a. ``numPoints`` : It is used to set the total number of data points to be used (sum of test and train data). pow(numPoints, dim) is the number of points, where ``dim`` is the dimension of the problem.


Example building and executing
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Navigate to the directory that holds the input .cpp file and run a make file with extension ``.mk`` in that directory to get your executable. The name of the executable is the name set in the ``EXECUTABLE`` flag of the make file. Then run the executable by giving the required run-time inputs. A sample run looks like this::

    user@computer NNCA$ cd NNCA2D
    user@computer NNCA2D$ make -f Makefile2D.mk clean
    user@computer NNCA2D$ make -f Makefile2D.mk
    user@computer NNCA2D$ ./testFMM2D 5 6 1 5 1




For the sake of this tutorial, we are going to be using the ``testHODLR3D.cpp`` file that is listed under ``examples/`` since it demonstrates the features of this library. For the most part, comments in the file demonstrate intended functionality. However, we go over the main functions that may be of interest to a user on this page.

**NOTE**: It is assumed that you have already completed the installation process of getting the dependencies.


Defining the kernel function
----------------------------

The matrix that needs to be solved for is abstracted through this derived class of ``kernel``. The main method that needs to be set for this class is ``getMatrixEntry`` which returns the entry at the :math:`i^{\mathrm{th}}` row and :math:`j^{\mathrm{th}}` column of the matrix. Provide the definition of the kernel function in the member function ``getMatrixEntry`` of class ``userkernel`` in the ``kernel.hpp`` file.


Defining location of particles/data points:
-------------------------------------------

NCA2DBebendorf, NCA2DZhao, NNCA2D
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The location of target points and charge points are to be defined in the member variable ``particles_X`` and ``particles_Y`` respectively of the ``userkernel`` object defined in the main function.

NNCA3D
^^^^^^

The location of particles (target points and charge points are assumed to be the same) are to be defined in the member variable ``particles`` of the ``userkernel`` object defined in the main function.

FSVM
^^^^

The location of data points are to be defined in the variable ``gridPoints`` of the main function.


Defining vector ``b``:
----------------------

NCA2DBebendorf, NCA2DZhao, NNCA2D, NNCA3D, NNCAnD
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The vector that is to be applied to the matrix to perform matrix-vector product is to be defined in the variable ``b``, which is a vector of ``Eigen::VectorXd`` type.

Defining vector ``x``:
----------------------

NNCA3D
^^^^^^
The RHS vector in the linear system ``Ax=b``, that is solved in the ``NNCA3D/testFMM3Dsolve.cpp`` file is to be defined in the variable ``x``, which is a vector of ``Eigen::VectorXd`` type.


Creating an Instance of ``AFMMnD``:
-----------------------------------

The main operations of the sub-project NNCAnD are carried out through the ``AFMMnD`` class defined in the file ``NNCA/NNCAnD/AFMMnD.hpp``. The parameters of its constructor are (ndim, N, nParticlesInLeafAlong1D, L, TOL_POW, Qchoice, nLevels, points), where N is the number of particles in the system, nLevels is the number of levels of the :math:`2^d` tree that is to be built, points is a vector holding the location of particles::

  AFMMnD* afmmnd = new AFMMnD(ndim, N, nParticlesInLeafAlong1D, L, TOL_POW, Qchoice, nLevels, points);


We will now proceed to demonstrate the individual methods available under this class.

``assembleAFMMnD``
^^^^^^^^^^^^^^^^^^

Assemble the matrix in NNCAnD structure; i.e. it finds the low rank representation of the low-rank matrix sub-blocks::

  afmmnd->assembleAFMMnD();

``matVecProduct``
^^^^^^^^^^^^^^^^^

Multiplies the matrix that is defined through object ``userkernel`` with the vector ``b`` and return the product::

  AFMM_Ab = afmmnd->matVecProduct(b);
