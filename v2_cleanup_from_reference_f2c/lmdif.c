
// lmdif.f -- translated by f2c (version 20100827).

#include "minpack_c.h"

void lmdif(lmdif_fcn_ptr fcn, const int m, const int n, double *x, double *fvec, const double ftol, const double xtol, const double gtol, const int maxfev,
           const double epsfcn, double *diag, const int mode, const double factor, const int nprint, int *info, int *nfev,
           double *fjac, const int ldfjac, int *ipvt, double *qtf, double *wa1, double *wa2, double *wa3, double *wa4)
{
    // The purpose of lmdif is to minimize the sum of the squares of
    // m nonlinear functions in n variables by a modification of
    // the Levenberg-Marquardt algorithm. The user must provide a
    // subroutine which calculates the functions. The Jacobian is
    // then calculated by a forward-difference approximation.

    // The subroutine statement is

    //   subroutine lmdif(fcn,m,n,x,fvec,ftol,xtol,gtol,maxfev,epsfcn,
    //                    diag,mode,factor,nprint,info,nfev,fjac,
    //                    ldfjac,ipvt,qtf,wa1,wa2,wa3,wa4)

    // where

    //   fcn is the name of the user-supplied subroutine which
    //     calculates the functions. fcn must be declared
    //     in an external statement in the user calling
    //     program, and should be written as follows.

    //     subroutine fcn(m,n,x,fvec,iflag)
    //     integer m,n,iflag
    //     double precision x(n),fvec(m)
    //     ----------
    //     calculate the functions at x and
    //     return this vector in fvec.
    //     ----------
    //     return
    //     end

    //     The value of iflag should not be changed by fcn unless
    //     the user wants to terminate execution of lmdif.
    //     In this case set iflag to a negative integer.

    //   m is a positive integer input variable set to the number
    //     of functions.

    //   n is a positive integer input variable set to the number
    //     of variables. n must not exceed m.

    //   x is an array of length n. On input x must contain
    //     an initial estimate of the solution vector. On output x
    //     contains the final estimate of the solution vector.

    //   fvec is an output array of length m which contains
    //     the functions evaluated at the output x.

    //   ftol is a nonnegative input variable. Termination
    //     occurs when both the actual and predicted relative
    //     reductions in the sum of squares are at most ftol.
    //     Therefore, ftol measures the relative error desired
    //     in the sum of squares.

    //   xtol is a nonnegative input variable. Termination
    //     occurs when the relative error between two consecutive
    //     iterates is at most xtol. Therefore, xtol measures the
    //     relative error desired in the approximate solution.

    //   gtol is a nonnegative input variable. Termination
    //     occurs when the cosine of the angle between fvec and
    //     any column of the Jacobian is at most gtol in absolute
    //     value. Therefore, gtol measures the orthogonality
    //     desired between the function vector and the columns
    //     of the Jacobian.

    //   maxfev is a positive integer input variable. Termination
    //     occurs when the number of calls to fcn is at least
    //     maxfev by the end of an iteration.

    //   epsfcn is an input variable used in determining a suitable
    //     step length for the forward-difference approximation. This
    //     approximation assumes that the relative errors in the
    //     functions are of the order of epsfcn. If epsfcn is less
    //     than the machine precision, it is assumed that the relative
    //     errors in the functions are of the order of the machine
    //     precision.

    //   diag is an array of length n. If mode = 1 (see
    //     below), diag is internally set. If mode = 2, diag
    //     must contain positive entries that serve as
    //     multiplicative scale factors for the variables.

    //   mode is an integer input variable. If mode = 1, the
    //     variables will be scaled internally. If mode = 2,
    //     the scaling is specified by the input diag. Other
    //     values of mode are equivalent to mode = 1.

    //   factor is a positive input variable used in determining the
    //     initial step bound. This bound is set to the product of
    //     factor and the Euclidean norm of diag*x if nonzero, or else
    //     to factor itself. In most cases factor should lie in the
    //     interval (.1, 100.). 100. is a generally recommended value.

    //   nprint is an integer input variable that enables controlled
    //     printing of iterates if it is positive. In this case,
    //     fcn is called with iflag = 0 at the beginning of the first
    //     iteration and every nprint iterations thereafter and
    //     immediately prior to return, with x and fvec available
    //     for printing. If nprint is not positive, no special calls
    //     of fcn with iflag = 0 are made.

    //   info is an integer output variable. If the user has
    //     terminated execution, info is set to the (negative)
    //     value of iflag. See description of fcn. Otherwise,
    //     info is set as follows.

    //     info = 0  improper input parameters.

    //     info = 1  both actual and predicted relative reductions
    //               in the sum of squares are at most ftol.

    //     info = 2  relative error between two consecutive iterates
    //               is at most xtol.

    //     info = 3  conditions for info = 1 and info = 2 both hold.

    //     info = 4  the cosine of the angle between fvec and any
    //               column of the Jacobian is at most gtol in
    //               absolute value.

    //     info = 5  number of calls to fcn has reached or
    //               exceeded maxfev.

    //     info = 6  ftol is too small. No further reduction in
    //               the sum of squares is possible.

    //     info = 7  xtol is too small. No further improvement in
    //               the approximate solution x is possible.

    //     info = 8  gtol is too small. fvec is orthogonal to the
    //               columns of the Jacobian to machine precision.

    //   nfev is an integer output variable set to the number of
    //     calls to fcn.

    //   fjac is an output m by n array. The upper n by n submatrix
    //     of fjac contains an upper triangular matrix r with
    //     diagonal elements of nonincreasing magnitude such that

    //            t     t           t
    //           p *(jac *jac)*p = r *r,

    //     where p is a permutation matrix and jac is the final
    //     calculated Jacobian. Column j of p is column ipvt(j)
    //     (see below) of the identity matrix. The lower trapezoidal
    //     part of fjac contains information generated during
    //     the computation of r.

    //   ldfjac is a positive integer input variable not less than m
    //     which specifies the leading dimension of the array fjac.

    //   ipvt is an integer output array of length n. ipvt
    //     defines a permutation matrix p such that jac*p = q*r,
    //     where jac is the final calculated Jacobian, q is
    //     orthogonal (not stored), and r is upper triangular
    //     with diagonal elements of nonincreasing magnitude.
    //     column j of p is column ipvt(j) of the identity matrix.

    //   qtf is an output array of length n which contains
    //     the first n elements of the vector (q transpose)*fvec.

    //   wa1, wa2, and wa3 are work arrays of length n.

    //   wa4 is a work array of length m.

    // Subprograms called:

    //   user-supplied ...... fcn
    //   MINPACK-supplied ... enorm, fdjac2, lmpar, qrfac

    // Argonne National Laboratory. MINPACK project. March 1980.
    // Burton S. Garbow, Kenneth E. Hillstrom, Jorge J. More

    const double p1    = 0.1;
    const double p5    = 0.5;
    const double p25   = 0.25;
    const double p75   = 0.75;
    const double p0001 = 1e-4;

    // Local variables.

    double delta = 0.0; // Initialize to suppress warning.
    double xnorm = 0.0; // Initialize to suppress warning.

    // Parameter adjustments.

    *info = 0;
    *nfev = 0;

    int iflag = 0;

    // Check the input parameters for errors.

    if (n <= 0 || m < n || ldfjac < m || ftol < 0.0 || xtol < 0.0 || gtol < 0.0 || maxfev <= 0 || factor <= 0.0)
    {
        goto L_TERMINATE;
    }

    if (mode == 2)
    {
        for (int j = 0; j < n; ++j)
        {
            if (diag[j] <= 0.0)
            {
                goto L_TERMINATE;
            }
        }
    }

    // Evaluate the function at the starting point
    // and calculate its norm.

    iflag = 1;
    (*fcn)(m, n, x, fvec, &iflag);
    *nfev = 1;

    if (iflag < 0)
    {
        goto L_TERMINATE;
    }

    double fnorm = enorm(m, fvec);

    // Initialize Levenberg-Marquardt parameter and iteration counter.

    double par = 0.0;
    int iter = 1;

    for (;;) // Outer loop.
    {
        // Beginning of the outer loop.

        // Calculate the Jacobian matrix.

        iflag = 2;
        fdjac2(fcn, m, n, x, fvec, fjac, ldfjac, &iflag, epsfcn, wa4);
        *nfev += n;

        if (iflag < 0)
        {
            goto L_TERMINATE;
        }

        // If requested, call fcn to enable printing of iterates.

        if (nprint > 0)
        {
            iflag = 0;
            if ((iter - 1) % nprint == 0)
            {
                (*fcn)(m, n, x, fvec, &iflag);
            }

            if (iflag < 0)
            {
                goto L_TERMINATE;
            }
        }

        // Compute the QR factorization of the Jacobian.

        qrfac(m, n, fjac, ldfjac, true, ipvt, n, wa1, wa2, wa3);

        // On the first iteration and if mode is 1, scale according
        // to the norms of the columns of the initial Jacobian.

        if (iter == 1)
        {
            if (mode != 2)
            {
                for (int j = 0; j < n; ++j)
                {
                    diag[j] = wa2[j];

                    if (wa2[j] == 0.0)
                    {
                        diag[j] = 1.0;
                    }
                }
            }

            // On the first iteration, calculate the norm of the scaled x
            // and initialize the step bound delta.

            for (int j = 0; j < n; ++j)
            {
                wa3[j] = diag[j] * x[j];
            }

            xnorm = enorm(n, wa3);
            delta = factor * xnorm;

            if (delta == 0.0)
            {
                delta = factor;
            }
        }

        // Form (q transpose)*fvec and store the first n components in
        // qtf.

        for (int i = 0; i < m; ++i)
        {
            wa4[i] = fvec[i];
        }

        for (int j = 0; j < n; ++j)
        {
            if (fjac[j + j * ldfjac] != 0.0)
            {
                double sum = 0.0;

                for (int i = j; i < m; ++i)
                {
                    sum += fjac[i + j * ldfjac] * wa4[i];
                }

                const double temp = -sum / fjac[j + j * ldfjac];

                for (int i = j; i < m; ++i)
                {
                    wa4[i] += fjac[i + j * ldfjac] * temp;
                }
            }

            fjac[j + j * ldfjac] = wa1[j];
            qtf[j] = wa4[j];
        }

        // Compute the norm of the scaled gradient.

        double gnorm = 0.0;

        if (fnorm != 0.0)
        {
            for (int j = 0; j < n; ++j)
            {
                const int l = ipvt[j] - PIVOT_OFFSET;

                if (wa2[l] != 0.0)
                {
                    double sum = 0.0;

                    for (int i = 0; i < j + 1; ++i)
                    {
                        sum += fjac[i + j * ldfjac] * (qtf[i] / fnorm);
                    }
                    gnorm = fmax(gnorm, fabs(sum / wa2[l]));
                }
            }
        }

        // Test for convergence of the gradient norm.

        if (gnorm <= gtol)
        {
            *info = 4;
        }

        if (*info != 0)
        {
            goto L_TERMINATE;
        }

        // Rescale if necessary.

        if (mode != 2)
        {
            for (int j = 0; j < n; ++j)
            {
                diag[j] = fmax(diag[j], wa2[j]);
            }
        }

        for (;;) // Inner loop.
        {
            // Beginning of the inner loop.

            // Determine the Levenberg-Marquardt parameter.

            lmpar(n, fjac, ldfjac, ipvt, diag, qtf, delta, &par, wa1, wa2, wa3, wa4);

            // Store the direction p and x + p. calculate the norm of p.

            for (int j = 0; j < n; ++j)
            {
                wa1[j] = -wa1[j];
                wa2[j] = x[j] + wa1[j];
                wa3[j] = diag[j] * wa1[j];
            }

            const double pnorm = enorm(n, wa3);

            // On the first iteration, adjust the initial step bound.

            if (iter == 1)
            {
                delta = fmin(delta, pnorm);
            }

            // Evaluate the function at x + p and calculate its norm.

            iflag = 1;
            (*fcn)(m, n, wa2, wa4, &iflag);
            ++(*nfev);

            if (iflag < 0)
            {
                goto L_TERMINATE;
            }

            const double fnorm1 = enorm(m, wa4);

            // Compute the scaled actual reduction.

            double actred = -1.0;

            if (p1 * fnorm1 < fnorm)
            {
                actred = 1.0 - square(fnorm1 / fnorm);
            }

            // Compute the scaled predicted reduction and
            // the scaled directional derivative.

            for (int j = 0; j < n; ++j)
            {
                wa3[j] = 0.0;

                const int l = ipvt[j] - PIVOT_OFFSET;

                const double temp = wa1[l];

                for (int i = 0; i < j + 1; ++i)
                {
                    wa3[i] += fjac[i + j * ldfjac] * temp;
                }
            }

            const double temp1 = enorm(n, wa3) / fnorm;
            const double temp2 = sqrt(par) * pnorm / fnorm;

            const double prered = square(temp1) + 2 * square(temp2);
            const double dirder = -(square(temp1) + square(temp2));

            // Compute the ratio of the actual to the predicted
            // reduction.

            double ratio = 0.0;

            if (prered != 0.0)
            {
                ratio = actred / prered;
            }

            // Update the step bound.

            if (ratio <= p25)
            {
                double temp;

                if (actred >= 0.0)
                {
                    temp = p5;
                }
                else
                {
                    temp = p5 * dirder / (dirder + p5 * actred);
                }

                if (p1 * fnorm1 >= fnorm || temp < p1)
                {
                    temp = p1;
                }

                delta = temp * fmin(delta, pnorm / p1);
                par /= temp;
            }
            else if (par == 0.0 || ratio >= p75)
            {
                delta = pnorm / p5;
                par = p5 * par;
            }

            // Test for successful iteration.

            if (ratio >= p0001)
            {
                // Successful iteration. Update x, fvec, and their norms.

                for (int j = 0; j < n; ++j)
                {
                    x[j] = wa2[j];
                    wa2[j] = diag[j] * x[j];
                }

                for (int i = 0; i < m; ++i)
                {
                    fvec[i] = wa4[i];
                }

                xnorm = enorm(n, wa2);
                fnorm = fnorm1;

                ++iter;
            }

            // Tests for convergence.

            if (fabs(actred) <= ftol && prered <= ftol && p5 * ratio <= 1.0)
            {
                *info = 1;
            }

            if (delta <= xtol * xnorm)
            {
                *info = 2;
            }

            if (fabs(actred) <= ftol && prered <= ftol && p5 * ratio <= 1.0 && *info == 2)
            {
                *info = 3;
            }

            if (*info != 0)
            {
                goto L_TERMINATE;
            }

            // Tests for termination and stringent tolerances.

            if (*nfev >= maxfev)
            {
                *info = 5;
            }

            if (fabs(actred) <= MACHINE_EPSILON && prered <= MACHINE_EPSILON && p5 * ratio <= 1.0)
            {
                *info = 6;
            }

            if (delta <= MACHINE_EPSILON * xnorm)
            {
                *info = 7;
            }

            if (gnorm <= MACHINE_EPSILON)
            {
                *info = 8;
            }

            if (*info != 0)
            {
                goto L_TERMINATE;
            }

            // End of the inner loop. Repeat if iteration unsuccessful.

            if (ratio >= p0001)
            {
                break;
            }
        }
        // End of the outer loop.
    }

L_TERMINATE:

    // Termination, either normal or user imposed.

    if (iflag < 0)
    {
        *info = iflag;
    }

    iflag = 0;
    if (nprint > 0)
    {
        (*fcn)(m, n, x, fvec, &iflag);
    }
}
