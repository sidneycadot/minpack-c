
// qrfac.f -- translated by f2c (version 20100827).

#include <stdbool.h>
#include <math.h>
#include "minpack_c.h"

inline doublereal square(const doublereal x)
{
    return x * x;
}

void qrfac(
        const integer  m,
        const integer  n,
        doublereal    *a,
        const integer  lda,
        const bool     pivot,
        integer       *ipvt,
        const integer  lipvt,
        doublereal    *rdiag,
        doublereal    *acnorm,
        doublereal    *wa
    )
{
    //     this subroutine uses householder transformations with column */
    //     pivoting (optional) to compute a qr factorization of the */
    //     m by n matrix a. that is, qrfac determines an orthogonal */
    //     matrix q, a permutation matrix p, and an upper trapezoidal */
    //     matrix r with diagonal elements of nonincreasing magnitude, */
    //     such that a*p = q*r. the householder transformation for */
    //     column k, k = 1,2,...,min(m,n), is of the form */

    //                           t */
    //           i - (1/u(k))*u*u */

    //     where u has zeros in the first k-1 positions. the form of */
    //     this transformation and the method of pivoting first */
    //     appeared in the corresponding linpack subroutine. */

    //     the subroutine statement is */

    //       subroutine qrfac(m,n,a,lda,pivot,ipvt,lipvt,rdiag,acnorm,wa) */

    //     where */

    //       m is a positive integer input variable set to the number */
    //         of rows of a. */

    //       n is a positive integer input variable set to the number */
    //         of columns of a. */

    //       a is an m by n array. on input a contains the matrix for */
    //         which the qr factorization is to be computed. on output */
    //         the strict upper trapezoidal part of a contains the strict */
    //         upper trapezoidal part of r, and the lower trapezoidal */
    //         part of a contains a factored form of q (the non-trivial */
    //         elements of the u vectors described above). */

    //       lda is a positive integer input variable not less than m */
    //         which specifies the leading dimension of the array a. */

    //       pivot is a logical input variable. if pivot is set true, */
    //         then column pivoting is enforced. if pivot is set false, */
    //         then no column pivoting is done. */

    //       ipvt is an integer output array of length lipvt. ipvt */
    //         defines the permutation matrix p such that a*p = q*r. */
    //         column j of p is column ipvt(j) of the identity matrix. */
    //         if pivot is false, ipvt is not referenced. */

    //       lipvt is a positive integer input variable. if pivot is false, */
    //         then lipvt may be as small as 1. if pivot is true, then */
    //         lipvt must be at least n. */

    //       rdiag is an output array of length n which contains the */
    //         diagonal elements of r. */

    //       acnorm is an output array of length n which contains the */
    //         norms of the corresponding columns of the input matrix a. */
    //         if this information is not needed, then acnorm can coincide */
    //         with rdiag. */

    //       wa is a work array of length n. if pivot is false, then wa */
    //         can coincide with rdiag. */

    //     subprograms called */

    //       minpack-supplied ... dpmpar,enorm */

    //       fortran-supplied ... dmax1,dsqrt,min0 */

    //     argonne national laboratory. minpack project. march 1980. */
    //     burton s. garbow, kenneth e. hillstrom, jorge j. more */

    (void)lipvt; // unused

    // Initialized data.

    const doublereal p05  = 0.05;

    // epsmch is the machine precision.

    const doublereal epsmch = dpmpar_1;

    // Compute the initial column norms and initialize several arrays.

    for (int j = 0; j < n; ++j)
    {
        wa[j] = rdiag[j] = acnorm[j] = enorm(m, &a[j * lda]);
        if (pivot)
        {
            ipvt[j] = j + 1;
        }
    }

    // Reduce a to r with householder transformations.

    const integer minmn = (m < n) ? m : n; // minimum of m and n.

    for (int j = 0; j < minmn; ++j)
    {
        if (pivot)
        {
            // Bring the column of largest norm into the pivot position.

            integer kmax = j;

            for (int k = j; k < n; ++k)
            {
                if (rdiag[k] > rdiag[kmax])
                {
                    kmax = k;
                }
            }

            if (kmax != j) // Exchange columns (kmax) and (j).
            {
                for (int i = 0; i < m; ++i)
                {
                    doublereal swap = a[i + j * lda];
                    a[i + j * lda] = a[i + kmax * lda];
                    a[i + kmax * lda] = swap;
                }

                rdiag[kmax] = rdiag[j];
                wa[kmax] = wa[j];

                int swap = ipvt[j];
                ipvt[j] = ipvt[kmax];
                ipvt[kmax] = swap;
            }
        }

        // Compute the householder transformation to reduce the
        // j-th column of a to a multiple of the j-th unit vector.

        doublereal ajnorm = enorm(m - j, &a[j + j * lda]);

        if (ajnorm != 0.0)
        {
            if (a[j + j * lda] < 0.0)
            {
                ajnorm = -ajnorm;
            }

            for (int i = j; i < m; ++i)
            {
                a[i + j * lda] /= ajnorm;
            }

            a[j + j * lda] += 1.0;

            // Apply the transformation to the remaining columns and update the norms.

            for (int k = (j+1); k < n; ++k)
            {
                doublereal sum = 0.0;
                for (int i = j; i < m; ++i)
                {
                    sum += a[i + j * lda] * a[i + k * lda];
                }
                doublereal temp = sum / a[j + j * lda];

                for (int i = j; i < m; ++i)
                {
                    a[i + k * lda] -= temp * a[i + j * lda];
                }

                if (pivot && rdiag[k] != 0.0)
                {
                    doublereal temp = a[j + k * lda] / rdiag[k];
                    rdiag[k] *= sqrt((fmax(0.0, 1.0 - square(temp))));

                    if (p05 * square(rdiag[k] / wa[k]) <= epsmch)
                    {
                        rdiag[k] = enorm(m - j + 1, &a[j + 1 + k * lda]);
                        wa[k] = rdiag[k];
                    }
                }
            }
        }
        rdiag[j] = -ajnorm;
    }
}
