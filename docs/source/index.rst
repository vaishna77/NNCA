.. role:: underline
    :class: underline

Welcome to NNCAlib's Documentation!
**************************************

About :math:`\texttt{NNCAlib}`:
==================================

NNCA (New Nested Cross Approximation) is an :math:`\mathcal{O}(N)` complexity algorithm to construct :math:`\mathcal{H}^{2}` matrices. :math:`\texttt{NNCAlib}` is a library to compute fast matrix-vector products using NNCA.

NCA is a technique developed to construct the :math:`\mathcal{H}^{2}` matrix approximation of a given matrix arising out of :math:`N` body problems. NNCA is a new variant of the Nested Cross Approximation (NCA). NNCA differs from the existing NCAs in the following ways. The existing NCA by `Bebendorf et al. <http://webdoc.sub.gwdg.de/ebook/serien/e/sfb611/491.pdf>`_ is not an analytic method, whereas NNCA is an algebraic method. Another existing NNCA by `Zhao et al. <https://engineering.purdue.edu/~djiao/publications/Zhao_NCA.pdf>`_ is an algebraic method and uses two tree traversals in its algorithm. Whereas NNCA uses a single tree traversal in its algorithm. It has been numerically shown in the NNCA article that NNCA performs computationally better than the existing NCAs.

The code is written in C++ and features an easy-to-use interface.

More about NNCA can be found in the `article <https://arxiv.org/pdf/2203.14832.pdf>`_.

Doc Contents
============
.. toctree::
   :maxdepth: 2
   :caption: Contents:

   installation
   tutorial
   reproducibility

Other Links
===========

Learn more about :math:`\texttt{NNCAlib}` by visiting the

* Code Repository: https://github.com/SAFRAN-LAB/NNCA
* Documentation: https://nnca.readthedocs.io/en/latest/
* Article: https://arxiv.org/pdf/2203.14832.pdf
