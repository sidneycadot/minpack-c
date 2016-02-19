
// lmpar.f -- translated by f2c (version 20100827).

#include "minpack_c.h"

void lmpar(const int n, double *r, const int ldr, const int *ipvt, const double *diag, const double *qtb, const double delta, double *par, double *x, double *sdiag, double *wa1, double *wa2)
{
    // Given an m by n matrix a, an n by n nonsingular diagonal
    // matrix d, an m-vector b, and a positive number delta,
    // the problem is to determine a value for the parameter
    // par such that if x solves the system

    //       a*x = b ,     sqrt(par)*d*x = 0 ,

    // in the least squares sense, and dxnorm is the Euclidean
    // norm of d*x, then either par is zero and

    //       (dxnorm-delta) .le. 0.1*delta ,

    // or par is positive and

    //       abs(dxnorm-delta) .le. 0.1*delta .

    // This subroutine completes the solution of the problem
    // if it is provided with the necessary information from the
    // QR factorization, with column pivoting, of a. That is, if
    // a*p = q*r, where p is a permutation matrix, q has orthogonal
    // columns, and r is an upper triangular matrix with diagonal
    // elements of nonincreasing magnitude, then lmpar expects
    // the full upper triangle of r, the permutation matrix p,
    // and the first n components of (q transpose)*b. On output
    // lmpar also provides an upper triangular matrix s such that

    //        t   t                   t
    //       p *(a *a + par*d*d)*p = s *s .

    // s is employed within lmpar and may be of separate interest.

    // Only a few iterations are generally needed for convergence
    // of the algorithm. If, however, the limit of 10 iterations
    // is reached, then the output par will contain the best
    // value obtained so far.

    // The subroutine statement is

    //   subroutine lmpar(n,r,ldr,ipvt,diag,qtb,delta,par,x,sdiag,
    //                    wa1,wa2)

    // where

    //   n is a positive integer input variable set to the order of r.

    //   r is an n by n array. On input the full upper triangle
    //     must contain the full upper triangle of the matrix r.
    //     On output the full upper triangle is unaltered, and the
    //     strict lower triangle contains the strict upper triangle
    //     (transposed) of the upper triangular matrix s.

    //   ldr is a positive integer input variable not less than n
    //     which specifies the leading dimension of the array r.

    //   ipvt is an integer input array of length n which defines the
    //     permutation matrix p such that a*p = q*r. Column j of p
    //     is column ipvt(j) of the identity matrix.

    //   diag is an input array of length n which must contain the
    //     diagonal elements of the matrix d.

    //   qtb is an input array of length n which must contain the first
    //     n elements of the vector (q transpose)*b.

    //   delta is a positive input variable which specifies an upper
    //     bound on the Euclidean norm of d*x.

    //   par is a nonnegative variable. On input par contains an
    //     initial estimate of the Levenberg-Marquardt parameter.
    //     on output par contains the final estimate.

    //   x is an output array of length n which contains the least
    //     squares solution of the system a*x = b, sqrt(par)*d*x = 0,
    //     for the output par.

    //   sdiag is an output array of length n which contains the
    //     diagonal elements of the upper triangular matrix s.

    //   wa1 and wa2 are work arrays of length n.

    // Subprograms called:

    //   MINPACK-supplied ... enorm, qrsolv

    // Argonne National Laboratory. MINPACK project. March 1980.
    // Burton S. Garbow, Kenneth E. Hillstrom, Jorge J. More

    const double p1   = 0.1;
    const double p001 = 0.001;

    // Compute and store in x the Gauss-Newton direction. If the
    // Jacobian is rank-deficient, obtain a least squares solution.

    int nsing = n;

    for (int j = 0; j < n; ++j)
    {
        wa1[j] = qtb[j];

        if (r[j + j * ldr] == 0.0 && nsing == n)
        {
            nsing = j;
        }

        if (nsing < n)
        {
            wa1[j] = 0.0;
        }
    }

    for (int k = 1; k <= nsing; ++k)
    {
        const int j = nsing - k;

        wa1[j] /= r[j + j * ldr];

        const double temp = wa1[j];

        for (int i = 0; i < j; ++i)
        {
            wa1[i] -= r[i + j * ldr] * temp;
        }
    }

    for (int j = 0; j < n; ++j)
    {
        const int l = ipvt[j];
        x[l] = wa1[j];
    }

    // Initialize the iteration counter.
    // Evaluate the function at the origin, and test
    // for acceptance of the Gauss-Newton direction.

    int iter = 0;

    for (int j = 0; j < n; ++j)
    {
        wa2[j] = diag[j] * x[j];
    }

    double dxnorm = enorm(n, wa2);
    double fp = dxnorm - delta;

    if (fp > p1 * delta)
    {
        // If the Jacobian is not rank deficient, the Newton
        // step provides a lower bound, parl, for the zero of
        // the function. Otherwise set this bound to zero.

        double parl = 0.0;

        if (nsing >= n)
        {
            for (int j = 0; j < n; ++j)
            {
                const int l = ipvt[j];

                wa1[j] = diag[l] * (wa2[l] / dxnorm);
            }

            for (int j = 0; j < n; ++j)
            {
                double sum = 0.0;

                for (int i = 0; i < j; ++i)
                {
                    sum += r[i + j * ldr] * wa1[i];
                }

                wa1[j] = (wa1[j] - sum) / r[j + j * ldr];
            }

            const double temp = enorm(n, wa1);

            parl = fp / delta / temp / temp;
        }

        // Calculate an upper bound, paru, for the zero of the function.

        for (int j = 0; j < n; ++j)
        {
            double sum = 0.0;

            for (int i = 0; i < j + 1; ++i)
            {
                sum += r[i + j * ldr] * qtb[i];
            }

            const int l = ipvt[j];

            wa1[j] = sum / diag[l];
        }

        const double gnorm = enorm(n, wa1);

        double paru = gnorm / delta;

        if (paru == 0.0)
        {
            paru = MACHINE_MINPOS / fmin(delta, p1);
        }

        // If the input par lies outside of the interval (parl, paru),
        // set par to the closer endpoint.

        *par = fmax(*par, parl);
        *par = fmin(*par, paru);

        if (*par == 0.0)
        {
            *par = gnorm / dxnorm;
        }

        for (;;)
        {
            // Beginning of an iteration.

            ++iter;

            // Evaluate the function at the current value of par.

            if (*par == 0.0)
            {
                *par = fmax(MACHINE_MINPOS, p001 * paru);
            }

            double temp = sqrt(*par);

            for (int j = 0; j < n; ++j)
            {
                wa1[j] = temp * diag[j];
            }

            qrsolv(n, r, ldr, ipvt, wa1, qtb, x, sdiag, wa2);

            for (int j = 0; j < n; ++j)
            {
                wa2[j] = diag[j] * x[j];
            }

            dxnorm = enorm(n, wa2);
            temp = fp;
            fp = dxnorm - delta;

            // If the function is small enough, accept the current value
            // of par. Also test for the exceptional cases where parl
            // is zero or the number of iterations has reached 10.

            if (fabs(fp) <= p1 * delta || (parl == 0.0 && fp <= temp && temp < 0.0) || iter == 10)
            {
                break;
            }

            // Compute the Newton correction.

            for (int j = 0; j < n; ++j)
            {
                const int l = ipvt[j];

                wa1[j] = diag[l] * (wa2[l] / dxnorm);
            }

            for (int j = 0; j < n; ++j)
            {
                const double temp = wa1[j] /= sdiag[j];

                for (int i = j + 1; i < n; ++i)
                {
                    wa1[i] -= r[i + j * ldr] * temp;
                }
            }

            temp = enorm(n, wa1);

            const double parc = fp / delta / temp / temp;

            // Depending on the sign of the function, update parl or paru.

            if (fp > 0.0)
            {
                parl = fmax(parl, *par);
            }
            else
            {
                paru = fmin(paru, *par);
            }

            // Compute an improved estimate for par.

            *par = fmax(parl, *par + parc);

            // End of an iteration.
        }
    }

    // Termination.

    if (iter == 0)
    {
        *par = 0.0;
    }
}
