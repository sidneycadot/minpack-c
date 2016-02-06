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

#include "minpack_c.h"

doublereal dpmpar_(integer *i)
{
    /* Initialized data */

  const doublereal constants[3] = {
    2.22044604926e-16,
    2.22507385852e-308, 
    1.79769313485e308
  };

  return constants[*i - 1];

} /* dpmpar_ */

