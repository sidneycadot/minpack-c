/* dpmpar.f -- translated by f2c (version 20100827).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/

#ifndef dpmpar_h
#define dpmpar_h

#include <assert.h>
#include "minpack_c.h"

//     Provides double precision machine parameters
//     when the appropriate set of data statements is activated (by
//     removing the c from column 1) and all other data statements are
//     rendered inactive. Most of the parameter values were obtained
//     from the corresponding Bell Laboratories Port Library function.

//     The value is

//       dpmpar_<i>

//     where

//       i is an integer input variable set to 1, 2, or 3 which
//         selects the desired machine parameter. If the machine has
//         t base b digits and its smallest and largest exponents are
//         emin and emax, respectively, then these parameters are

//         dpmpar_1 = b**(1 - t), the machine precision,
//         dpmpar_2 = b**(emin - 1), the smallest magnitude, 
//         dpmpar_3 = b**emax*(1 - b**(-t)), the largest magnitude.

//     Argonne National Laboratory. MINPACK Project. November 1996.
//     Burton S. Garbow, Kenneth E. Hillstrom, Jorge J. More'

static const doublereal dpmpar_1 = 2.22044604926e-16;
static const doublereal dpmpar_2 = 2.22507385852e-308;
static const doublereal dpmpar_3 = 1.79769313485e308;

#endif // dpmpar_h
