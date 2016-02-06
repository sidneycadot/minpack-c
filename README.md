MINPACK in C
============

This repository contains a conversion of MINPACK, the venerable Fortran package that provides
functionality for function minimization, to C.

We currently focus on getting the lmdif() function to work, which executes the Levenberg-Marquard nonlinear least-squares
algorithm without the need for an analytical function to evaluate the Jacobian (it uses a finite-difference approximation instead).

In a following iteration we will add support for the lmder() function, which supports an analytical Jacobian.
This should be more robust and possibly faster.

Directories
-----------

* docs : Contains original documentation of the MINPACK package (see External Links).
* v0_original_fortran : Contains the original Fortran-66 code as downloaded from the Netlib website.
* v1_reference_f2c : An automated translation of the Fortran code to C, (using "f2c"). Includes minimal changes to get it to run independently of libf2c.
* v2_cleanup_from_reference_f2c : A cleaned-up, but functionally identical, version of the reference_f2c version.
* v3_simplified_api : A version of the C code that omits obscure functionality (such as printing during evaluation) and adds API improvements, including enum types and a user function context parameter.
* v4_python : A version of the solver in Python, aimed at understanding what is going on mathematically.
* testdata : The codes from v1 onward come with a testcase that executes a fit. This directory contains the corresponding test data.
* misc : Miscellaneous items.

External links
--------------

* https://en.wikipedia.org/wiki/MINPACK
* http://www.netlib.org/minpack/
* http://www.mcs.anl.gov/~more/ANL8074a.pdf
* http://www.mcs.anl.gov/~more/ANL8074b.pdf
