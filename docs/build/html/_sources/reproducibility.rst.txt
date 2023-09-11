***************
Reproducibility
***************

This part of the documentation helps in reproducing the results shown in the article titled **A new nested cross approximation**, authored by **Vaishnavi Gujjula, Sivaram Ambikasaran**.
The code is available as an open source library and can be found `here <https://github.com/SAFRAN-LAB/NNCA>`_.

Experiment 1: Matrix-vector product with uniform distribution of particles in 2D and its comparison with the existing NCAs
--------------------------------------------------------------------------------------------------------------------------

Figure 2
^^^^^^^^

I.  To obtain NNCA plots of Figure 2 of the article (labeled as 'Proposed' in Figure 2), which is about the various numerical benchmarks of NNCA vs :math:`\epsilon_{ACA}`, the following steps are to be followed

  #. Set the kernel function as :math:`\left(\frac{r(\log(r)-1)}{a(\log(a)-1)}\right)\chi_{r<a}+\left(\frac{\log(r)}{\log(a)}\right)\chi_{r\geq a}` in the ``getMatrixEntry`` function of the ``userkernel`` class located in the ``NNCA2D/kernel.hpp`` file.
  #. Run the file ``NNCA2D/testFMM2D.cpp`` using the ``NNCA2D/Makefile2D.mk``. After successful compilation, an executable by name ``NNCA2D/testFMM2D`` gets generated. Then input the following at run-time

    a. ``nLevels`` = 6
    b. ``nParticlesInLeafAlong1D`` = 5
    c. ``L`` = 1.0
    d. ``TOL_POW`` : varied from 3 through 10
    e. ``yesToUniformDistribution`` = 1

    For example, run the following command::

       NNCA2D/testFMM2D 6 5 1.0 3 1

II. To obtain NCA by Zhao plots of Figure 2 of the article (labeled as 'Zhao' in Figure 2), which is about the various numerical benchmarks of NCA by Zhao vs :math:`\epsilon_{ACA}`, the following steps are to be followed

  #. Set the kernel function as :math:`\left(\frac{r(\log(r)-1)}{a(\log(a)-1)}\right)\chi_{r<a}+\left(\frac{\log(r)}{\log(a)}\right)\chi_{r\geq a}` in the ``getMatrixEntry`` function of the ``userkernel`` class located in the ``NCA2DZhao/kernel.hpp`` file.
  #. Run the file ``NCA2DZhao/testFMM2D.cpp`` using the ``NCA2DZhao/Makefile2D.mk``. After successful compilation, an executable by name ``NCA2DZhao/testFMM2D`` gets generated. Then input the following at run-time

    a. ``nLevels`` = 6
    b. ``nParticlesInLeafAlong1D`` = 5
    c. ``L`` = 1.0
    d. ``TOL_POW`` : varied from 3 through 10
    e. ``yesToUniformDistribution`` = 1

    For example, run the following command::

       NCA2DZhao/testFMM2D 6 5 1.0 3 1

III. To obtain NCA by Bebendorf plots of Figure 2 of the article (labeled as 'Bebendorf' in Figure 2), which is about the various numerical benchmarks of NCA by Bebendorf vs :math:`\epsilon_{ACA}`, the following steps are to be followed

 #. Set the kernel function as :math:`\left(\frac{r(\log(r)-1)}{a(\log(a)-1)}\right)\chi_{r<a}+\left(\frac{\log(r)}{\log(a)}\right)\chi_{r\geq a}` in the ``getMatrixEntry`` function of the ``userkernel`` class located in the ``NCA2DBebendorf/kernel.hpp`` file.
 #. Run the file ``NCA2DBebendorf/testFMM2D.cpp`` using the ``NCA2DBebendorf/Makefile2D.mk``. After successful compilation, an executable by name ``NCA2DBebendorf/testFMM2D`` gets generated. Then input the following at run-time

   a. ``nLevels`` = 6
   b. ``nParticlesInLeafAlong1D`` = 5
   c. ``L`` = 1.0
   d. ``TOL_POW`` : varied from 3 through 10
   e. ``yesToUniformDistribution`` = 1

   For example, run the following command::

      NCA2DBebendorf/testFMM2D 6 5 1.0 3 1

Figure 3
^^^^^^^^

To obtain Figure 3 of the article, follow the same steps as that of Figure 2 except for the inputs at run-time. The following inputs are to be given at run-time for NCA by Zhao, NCA by Bebendorf, and NNCA.

1. ``nLevels`` : varied from 2 through 6
2. ``nParticlesInLeafAlong1D`` = 5
3. ``L`` = 1.0
4. ``TOL_POW`` = 9
5. ``yesToUniformDistribution`` = 1


Figure 4
^^^^^^^^

To obtain Figure 4, follow the same steps as that of Figure 2 except for the definition of the kernel function. Set the kernel function as :math:`\left(\frac{r}{a}\right)\chi_{r<a}+\left(\frac{a}{r}\right)\chi_{r\geq a}` in the ``getMatrixEntry`` function of the ``userkernel`` class located in the ``kernel.hpp`` file.

Figure 5
^^^^^^^^

To obtain Figure 5 of the article, follow the same steps as that of Figure 3 except for the definition of the kernel function. Set the kernel function as :math:`\left(\frac{r}{a}\right)\chi_{r<a}+\left(\frac{a}{r}\right)\chi_{r\geq a}` in the ``getMatrixEntry`` function of the ``userkernel`` class located in the ``kernel.hpp`` file.



Experiment 2: Matrix-vector product with Chebyshev distribution of particles in 2D
----------------------------------------------------------------------------------

Figure 6
^^^^^^^^

Line 18 of ``NNCA2D/testFMM2D.cpp`` file, which is ``A->set_Uniform_Nodes();``, has to be commented and uncomment Line 19 of ``NNCA2D/testFMM2D.cpp`` file, which is ``A->set_Standard_Cheb_Nodes();``. Further follow the below steps


#. Set the kernel function as :math:`\left(\frac{r(\log(r)-1)}{a(\log(a)-1)}\right)\chi_{r<a}+\left(\frac{\log(r)}{\log(a)}\right)\chi_{r\geq a}` in the ``getMatrixEntry`` function of the ``userkernel`` class located in the ``NNCA2D/kernel.hpp`` file.
#. Run the file ``NNCA2D/testFMM2D.cpp`` using the ``NNCA2D/Makefile2D.mk``. After successful compilation, an executable by name ``NNCA2D/testFMM2D`` gets generated. Then input the following at run-time

  a. ``nLevels`` = 6
  b. ``nParticlesInLeafAlong1D`` = 5
  c. ``L`` = 1.0
  d. ``TOL_POW`` : varied from 3 through 10
  e. ``yesToUniformDistribution`` = 0

  For example, run the following command::

     NNCA2D/testFMM2D 6 5 1.0 3 0

 Figure 7
 ^^^^^^^^

 To obtain Figure 7 of the article, follow the same steps as that of Figure 6 except for the inputs at run-time. The following inputs are to be given at run-time.

 1. ``nLevels`` : varied from 2 through 6
 2. ``nParticlesInLeafAlong1D`` = 5
 3. ``L`` = 1.0
 4. ``TOL_POW`` = 9
 5. ``yesToUniformDistribution`` = 0

 For example, run the following command::

    NNCA2D/testFMM2D 6 5 1.0 9 0


Experiment 3: Matrix-vector product with a uniform distribution of particles in 3D
----------------------------------------------------------------------------------

Figure 8
^^^^^^^^
To obtain NNCA plots of Figure 8 of the article, which is about the various numerical benchmarks of NNCA in 3D vs :math:`\epsilon_{ACA}`, the following steps are to be followed

  #. To set the kernel function as :math:`\left(\frac{r}{a}\right)\chi_{r<a}+\left(\frac{a}{r}\right)\chi_{r\geq a}`, choose ``Qchoice`` to be 17 at run time.
  #. Run the file ``NNCA3D/testFMM3D.cpp`` using the ``NNCA3D/Makefile3D.mk``. After successful compilation, an executable by name ``NNCA3D/testFMM3D`` gets generated. Then input the following at run-time

    1. ``nLevels`` : varied from 2 through 6
    2. ``nParticlesInLeafAlong1D`` = 5
    3. ``L`` = 1.0
    4. ``TOL_POW`` = 7
    5. ``Qchoice`` = 17
    6. ``yesToUniformDistribution`` = 1

    For example, run the following command::

       NNCA3D/testFMM3D 6 5 1.0 7 17 1

Experiment 4: Matrix-vector product with Chebyshev distribution of particles in 3D
----------------------------------------------------------------------------------

Figure 9
^^^^^^^^

Line 21 of ``NNCA3D/testFMM3D.cpp`` file, which is ``A->set_Uniform_Nodes(h);``, has to be commented and uncomment Line 22 of ``NNCA3D/testFMM3D.cpp`` file, which is ``A->set_Standard_Cheb_Nodes();``. Then follow the below steps

#. To set the kernel function as :math:`\left(\frac{r}{a}\right)\chi_{r<a}+\left(\frac{a}{r}\right)\chi_{r\geq a}`, choose ``Qchoice`` to be 17 at run time.
#. Run the file ``NNCA3D/testFMM3D.cpp`` using the ``NNCA3D/Makefile3D.mk``. After successful compilation, an executable by name ``NNCA3D/testFMM3D`` gets generated. Then input the following at run-time

  1. ``nLevels`` : varied from 2 through 6
  2. ``nParticlesInLeafAlong1D`` = 5
  3. ``L`` = 1.0
  4. ``TOL_POW`` = 7
  5. ``Qchoice`` = 17
  6. ``yesToUniformDistribution`` = 0

  For example, run the following command::

     NNCA3D/testFMM3D 6 5 1.0 7 17 0



Experiment 5: Integral equation solver in 3D
--------------------------------------------

Figure 10
^^^^^^^^^
To obtain NNCA plots of Figure 10 of the article, which is about the various numerical benchmarks of NNCA in 3D vs :math:`\epsilon_{ACA}`, the following steps are to be followed

  #. To set the kernel function as :math:`\left(\frac{r}{a}\right)\chi_{r<a}+\left(\frac{a}{r}\right)\chi_{r\geq a}`, choose ``Qchoice`` to be 17 at run time.
  #. Run the file ``NNCA3D/testFMM3Dsolve.cpp`` using the ``NNCA3D/Makefile3Dsolve.mk``. After successful compilation, an executable by name ``NNCA3D/testFMM3Dsolve`` gets generated. Then input the following at run-time

    1. ``nLevels`` : varied from 2 through 6
    2. ``nParticlesInLeafAlong1D`` = 5
    3. ``L`` = 1.0
    4. ``TOL_POW`` = 7

    For example, run the following command::

       NNCA3D/testFMM3D 6 5 1.0 7


Experiment 6: Kernel SVM in 4D
------------------------------

Table 5 and Figure 11
^^^^^^^^^^^^^^^^^^^^^
To obtain the results of Table 5 and Figure 11 of the article, which is about the performance of Fast SVM vs Normal SVM in 2D, the following steps are to be followed

#. Set the ``KERNEL``, ``DIM``, ``MATVEC`` flags in the  ``Makefile.mk`` file as below

   #. To get the results corresponding to the Gaussian kernel and Fast SVM use ``USE_Gaussian``, ``USE_DIM2``, and ``USE_AFMMnD`` flags in the make file ``Makefile.mk``.
   #. To get the results corresponding to the Gaussian kernel and Normal SVM use ``USE_Gaussian``, ``USE_DIM2``, and ``USE_directMatVec`` flags in the make file ``Makefile.mk``.
   #. To get the results corresponding to the Matern kernel and Fast SVM use ``USE_Matern``, ``USE_DIM2``, and ``USE_AFMMnD`` flags in the make file ``Makefile.mk``.
   #. To get the results corresponding to the Matern kernel and Normal SVM use ``USE_Matern``, ``USE_DIM2``, and ``USE_directMatVec`` flags in the make file ``Makefile.mk``.

#. Further run the file ``FSVM/.cpp`` using the ``FSVM/Makefile3D.mk``. After successful compilation, an executable by name ``FSVM/testFSVM`` gets generated. Then input ``numPoints`` at runtime. ``numPoints`` is the number of points in one dimension. In 2D, the number of data points (sum of train and test data points) is square of ``numPoints`` and  in 4D it is quadraple of ``numPoints``.

   To reproduce the results, run the following command::

      FSVM/testFSVM 75

This will print various benchmarks. Table 5 of the article is populated using these results. To obtain Figure 5, the following steps are to be followed

#. Run the ``FSVM/matlab/plot_decision2D.m`` file. It prompts the user to input the kernel being used. The prompt says ``Enter 0 for Matern and 1 for Gaussian `` . When this is inputed, the user could see the figures.

Table 6
^^^^^^^
To obtain the results of Table 6 of the article, which is about the performance of Fast SVM vs Normal SVM in 4D, the following steps are to be followed

#. Set the ``KERNEL`` flag as ``USE_Matern``, ``DIM`` as ``USE_DIM2`` in the ``Makefile.mk``. Set the ``MATVEC`` flag in the ``Makefile.mk`` as below
   #. To get the results corresponding Fast SVM use ``USE_AFMMnD``.
   #. To get the results corresponding to Normal SVM use ``USE_directMatVec``.
#. Further run the file ``FSVM/.cpp`` using the ``FSVM/Makefile3D.mk``. After successful compilation, an executable by name ``FSVM/testFSVM`` gets generated. Then input the following at run-time

  ``numPoints`` : 8 through 11

   For example, run the following command::

      FSVM/testFSVM 8
