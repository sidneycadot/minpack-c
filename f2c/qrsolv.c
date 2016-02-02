
// qrsolv.f -- translated by f2c (version 20100827).

#include <math.h>
#include "minpack_c.h"

void qrsolv(
        const int n,
        double *r,
        const int ldr,
        const int *ipvt,
        double *diag,
        double *qtb,
        double *x,
        double *sdiag,
        double *wa
    )
{
    // Given an m by n matrix a, an n by n diagonal matrix d,
    // and an m-vector b, the problem is to determine an x which
    // solves the system

    //       a*x = b ,     d*x = 0 ,

    // in the least squares sense.

    // This subroutine completes the solution of the problem
    // if it is provided with the necessary information from the
    // qr factorization, with column pivoting, of a. That is, if
    // a*p = q*r, where p is a permutation matrix, q has orthogonal
    // columns, and r is an upper triangular matrix with diagonal
    // elements of nonincreasing magnitude, then qrsolv expects
    // the full upper triangle of r, the permutation matrix p,
    // and the first n components of (q transpose)*b. The system
    // a*x = b, d*x = 0, is then equivalent to

    //              t       t
    //       r*z = q *b ,  p *d*p*z = 0 ,

    // where x = p*z. If this system does not have full rank,
    // then a least squares solution is obtained. On output qrsolv
    // also provides an upper triangular matrix s such that

    //        t   t               t
    //       p *(a *a + d*d)*p = s *s .

    // s is computed within qrsolv and may be of separate interest.

    // The subroutine statement is

    //   subroutine qrsolv(n,r,ldr,ipvt,diag,qtb,x,sdiag,wa)

    // where

    //   n is a positive integer input variable set to the order of r.

    //   r is an n by n array. on input the full upper triangle
    //     must contain the full upper triangle of the matrix r.
    //     On output the full upper triangle is unaltered, and the
    //     strict lower triangle contains the strict upper triangle
    //     (transposed) of the upper triangular matrix s.

    //   ldr is a positive integer input variable not less than n
    //     which specifies the leading dimension of the array r.

    //   ipvt is an integer input array of length n which defines the
    //     permutation matrix p such that a*p = q*r. column j of p
    //     is column ipvt(j) of the identity matrix.

    //   diag is an input array of length n which must contain the
    //     diagonal elements of the matrix d.

    //   qtb is an input array of length n which must contain the first
    //     n elements of the vector (q transpose)*b.

    //   x is an output array of length n which contains the least
    //     squares solution of the system a*x = b, d*x = 0.

    //   sdiag is an output array of length n which contains the
    //     diagonal elements of the upper triangular matrix s.

    //   wa is a work array of length n.

    // Argonne National Laboratory. MINPACK project. March 1980.
    // Burton S. Garbow, Kenneth E. Hillstrom, Jorge J. More

    // Copy r and (q transpose) * b to preserve input and initialize s.
    // In particular, save the diagonal elements of r in x.

    for (int j = 0; j < n; ++j)
    {
        for (int i = j; i < n; ++i)
        {
            r[i + j * ldr] = r[j + i * ldr];
        }
        x[j] = r[j + j * ldr];
        wa[j] = qtb[j];
    }

    // Eliminate the diagonal matrix d using a givens rotation.

    for (int j = 0; j < n; ++j)
    {
        // Prepare the row of d to be eliminated, locating the
        // diagonal element using p from the qr factorization.

        const int l = ipvt[j] - 1;

        if (diag[l] != 0.0)
        {
            for (int k = j; k < n; ++k)
            {
                sdiag[k] = 0.0;
            }

            sdiag[j] = diag[l];

            // The transformations to eliminate the row of d
            // modify only a single element of (q transpose) * b
            // beyond the first n, which is initially zero.

            double qtbpj = 0.0;

            for (int k = j; k < n; ++k)
            {
                // Determine a givens rotation which eliminates the
                // appropriate element in the current row of d.

                if (sdiag[k] != 0.0)
                {
                    double sin__, cos__;

                    if (fabs(r[k + k * ldr]) < fabs(sdiag[k]))
                    {

                        const double cotan = r[k + k * ldr] / sdiag[k];

                        sin__ = 1.0 / sqrt(1.0 + cotan * cotan);
                        cos__ = sin__ * cotan;
                    }
                    else
                    {
                        const double tan__ = sdiag[k] / r[k + k * ldr];
                        cos__ = 1.0 / sqrt(1.0 + tan__ * tan__);
                        sin__ = cos__ * tan__;
                    }

                    // Compute the modified diagonal element of r and
                    // the modified element of ((q transpose)*b,0).

                    r[k + k * ldr] = cos__ * r[k + k * ldr] + sin__ * sdiag[k];

                    const double ww = wa[k];

                    wa[k] =  cos__ * ww + sin__ * qtbpj;
                    qtbpj = -sin__ * ww + cos__ * qtbpj;

                    // Accumulate the tranformation in the row of s.

                    for (int i = k+1; i < n; ++i)
                    {
                        const double rr = r[i + k * ldr];

                        r[i + k * ldr] =  cos__ * rr + sin__ * sdiag[i];
                        sdiag[i]       = -sin__ * rr + cos__ * sdiag[i];
                    }
                }
            }
        }

        // Store the diagonal element of s and restore
        // the corresponding diagonal element of r.

        sdiag[j] = r[j + j * ldr];
        r[j + j * ldr] = x[j];
    }

    // Solve the triangular system for z. if the system is
    // singular, then obtain a least squares solution.

    int nsing = n;

    for (int j = 0; j < n; ++j)
    {
        if (sdiag[j] == 0.0 && nsing == n)
        {
            nsing = j;
        }

        if (nsing < n)
        {
            wa[j] = 0.0;
        }
    }

    for (int k = 0; k < nsing; ++k)
    {
        int j = nsing - k;

        double sum = 0.0;

        for (int i = j; i < nsing; ++i)
        {
            sum += r[i + j * ldr] * wa[i];
        }

        wa[j] = (wa[j] - sum) / sdiag[j];
    }

    // Permute the components of z back to components of x.

    for (int j = 0; j < n; ++j)
    {
        int l = ipvt[j] - 1;
        x[l] = wa[j];
    }
}
