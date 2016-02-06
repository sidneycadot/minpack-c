
// enorm.f -- translated by f2c (version 20100827).

#include "minpack_c.h"

double enorm(const int n, const double *x)
{
    // Given an n-vector x, this function calculates the
    // Euclidean norm of x.

    // The Euclidean norm is computed by accumulating the sum of
    // squares in three different sums. The sums of squares for the
    // small and large components are scaled so that no overflows
    // occur. Non-destructive underflows are permitted. Underflows
    // and overflows do not occur in the computation of the unscaled
    // sum of squares for the intermediate components.
    // The definitions of small, intermediate and large components
    // depend on two constants, rdwarf and rgiant. The main
    // restrictions on these constants are that rdwarf**2 not
    // underflow and rgiant**2 not overflow. The constants
    // given here are suitable for every known computer.

    // The function statement is

    //   double enorm(n, x)

    // where

    //   n is a positive integer input variable.

    //   x is an input array of length n.

    // Argonne National Laboratory. MINPACK project. March 1980.
    // Burton S. Garbow, Kenneth E. Hillstrom, Jorge J. More

    const double rdwarf = 3.834e-20;
    const double rgiant = 1.304e+19;

    const double agiant = rgiant / n;

    double s1 = 0.0;
    double s2 = 0.0;
    double s3 = 0.0;

    double x1max = 0.0;
    double x3max = 0.0;

    for (int i = 0; i < n; ++i)
    {
        const double xabs = fabs(x[i]);

        if (xabs <= rdwarf)
        {
            // Sum for large components.

            if (xabs > x1max)
            {
                s1 = 1.0 + s1 * square(x1max / xabs);
                x1max = xabs;
            }
            else
            {
                s1 += square(xabs / x1max);
            }
        }
        else if (xabs >= agiant)
        {
            // Sum for small components.

            if (xabs > x3max)
            {
                s3 = 1.0 + s3 * square(x3max / xabs);
                x3max = xabs;
            }
            else if (xabs != 0.0)
            {
                s3 += square(xabs / x3max);
            }
        }
        else
        {
            // Sum for intermediate components.

            s2 += square(xabs);
        }
    }

    // Calculation of Euclidean norm.

    double norm;

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
