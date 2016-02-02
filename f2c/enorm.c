/* enorm.f -- translated by f2c (version 20100827).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/

#include <math.h>
#include "minpack_c.h"

inline doublereal square(const doublereal x)
{
    return x * x;
}

doublereal enorm(const integer n, const doublereal *x)
{
    //     function enorm */

    //     given an n-vector x, this function calculates the
    //     euclidean norm of x.

    //     the euclidean norm is computed by accumulating the sum of
    //     squares in three different sums. the sums of squares for the
    //     small and large components are scaled so that no overflows
    //     occur. non-destructive underflows are permitted. underflows
    //     and overflows do not occur in the computation of the unscaled
    //     sum of squares for the intermediate components.
    //     the definitions of small, intermediate and large components
    //     depend on two constants, rdwarf and rgiant. the main
    //     restrictions on these constants are that rdwarf**2 not
    //     underflow and rgiant**2 not overflow. the constants
    //     given here are suitable for every known computer.

    //     the function statement is

    //       double precision function enorm(n,x)

    //     where

    //       n is a positive integer input variable.

    //       x is an input array of length n.

    //     subprograms called

    //       fortran-supplied ... dabs,dsqrt

    //     argonne national laboratory. minpack project. march 1980.
    //     burton s. garbow, kenneth e. hillstrom, jorge j. more

    // Initialized data

    const doublereal rdwarf = 3.834e-20;
    const doublereal rgiant = 1.304e19;

    // Function body

    doublereal s1    = 0.0;
    doublereal s2    = 0.0;
    doublereal s3    = 0.0;
    doublereal x1max = 0.0;
    doublereal x3max = 0.0;

    const doublereal agiant = rgiant / n;

    for (int i = 0; i < n; ++i)
    {
        doublereal xabs = fabs(x[i]);

        if (xabs <= rdwarf)
        {
            // Sum for small components.

            if (xabs <= x3max)
            {
                if (xabs != 0.0)
                {
                    s3 += square(xabs / x3max);
                }
            }
            else
            {
                s3 = 1.0 + s3 * square(x3max / xabs);
                x3max = xabs;
            }
        }
        else if (xabs < agiant)
        {
            // Sum for intermediate components.
            s2 += square(xabs);
        }
        else
        {
            // Sum for large components.

            if (xabs <= x1max)
            {
                s1 += square(xabs / x1max);
            }
            else
            {
                s1 = 1.0 + s1 * square(x1max / xabs);
                x1max = xabs;
            }
        }
    }

    // Calculation of norm.

    doublereal norm;

    if (s1 != 0.0)
    {
        norm = x1max * sqrt(s1 + s2 / x1max / x1max);
    }
    else if (s2 != 0.0)
    {
        if (s2 >= x3max)
        {
            norm = sqrt(s2 * (1.0 + x3max / s2 * (x3max * s3)));
        }
        else
        {
            norm = sqrt(x3max * (s2 / x3max + x3max * s3));
        }
    }
    else
    {
        norm = x3max * sqrt(s3);
    }

    return norm;
}
