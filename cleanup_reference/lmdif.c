// lmdif.f -- translated by f2c (version 20100827).

#include "my_include.h"

void lmdif(U_fp fcn, const int m, const int n, double *x, double *fvec, double *ftol, double *xtol, double *gtol, int *maxfev,
           double *epsfcn, double *diag, int *mode, double *factor, int *nprint, int *info, int *nfev,
           double *fjac, const int ldfjac, int *ipvt, double *qtf, double *wa1, double *wa2, double *wa3, double *wa4)
{
    // the purpose of lmdif is to minimize the sum of the squares of
    // m nonlinear functions in n variables by a modification of
    // the levenberg-marquardt algorithm. the user must provide a
    // subroutine which calculates the functions. the jacobian is
    // then calculated by a forward-difference approximation.

    // the subroutine statement is

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

    //     the value of iflag should not be changed by fcn unless
    //     the user wants to terminate execution of lmdif.
    //     in this case set iflag to a negative integer.

    //   m is a positive integer input variable set to the number
    //     of functions.

    //   n is a positive integer input variable set to the number
    //     of variables. n must not exceed m.

    //   x is an array of length n. on input x must contain
    //     an initial estimate of the solution vector. on output x
    //     contains the final estimate of the solution vector.

    //   fvec is an output array of length m which contains
    //     the functions evaluated at the output x.

    //   ftol is a nonnegative input variable. termination
    //     occurs when both the actual and predicted relative
    //     reductions in the sum of squares are at most ftol.
    //     therefore, ftol measures the relative error desired
    //     in the sum of squares.

    //   xtol is a nonnegative input variable. termination
    //     occurs when the relative error between two consecutive
    //     iterates is at most xtol. therefore, xtol measures the
    //     relative error desired in the approximate solution.

    //   gtol is a nonnegative input variable. termination
    //     occurs when the cosine of the angle between fvec and
    //     any column of the jacobian is at most gtol in absolute
    //     value. therefore, gtol measures the orthogonality
    //     desired between the function vector and the columns
    //     of the jacobian.

    //   maxfev is a positive integer input variable. termination
    //     occurs when the number of calls to fcn is at least
    //     maxfev by the end of an iteration.

    //   epsfcn is an input variable used in determining a suitable
    //     step length for the forward-difference approximation. this
    //     approximation assumes that the relative errors in the
    //     functions are of the order of epsfcn. if epsfcn is less
    //     than the machine precision, it is assumed that the relative
    //     errors in the functions are of the order of the machine
    //     precision.

    //   diag is an array of length n. if mode = 1 (see
    //     below), diag is internally set. if mode = 2, diag
    //     must contain positive entries that serve as
    //     multiplicative scale factors for the variables.

    //   mode is an integer input variable. if mode = 1, the
    //     variables will be scaled internally. if mode = 2,
    //     the scaling is specified by the input diag. other
    //     values of mode are equivalent to mode = 1.

    //   factor is a positive input variable used in determining the
    //     initial step bound. this bound is set to the product of
    //     factor and the euclidean norm of diag*x if nonzero, or else
    //     to factor itself. in most cases factor should lie in the
    //     interval (.1,100.). 100. is a generally recommended value.

    //   nprint is an integer input variable that enables controlled
    //     printing of iterates if it is positive. in this case,
    //     fcn is called with iflag = 0 at the beginning of the first
    //     iteration and every nprint iterations thereafter and
    //     immediately prior to return, with x and fvec available
    //     for printing. if nprint is not positive, no special calls
    //     of fcn with iflag = 0 are made.

    //   info is an integer output variable. if the user has
    //     terminated execution, info is set to the (negative)
    //     value of iflag. see description of fcn. otherwise,
    //     info is set as follows.

    //     info = 0  improper input parameters.

    //     info = 1  both actual and predicted relative reductions
    //               in the sum of squares are at most ftol.

    //     info = 2  relative error between two consecutive iterates
    //               is at most xtol.

    //     info = 3  conditions for info = 1 and info = 2 both hold.

    //     info = 4  the cosine of the angle between fvec and any
    //               column of the jacobian is at most gtol in
    //               absolute value.

    //     info = 5  number of calls to fcn has reached or
    //               exceeded maxfev.

    //     info = 6  ftol is too small. no further reduction in
    //               the sum of squares is possible.

    //     info = 7  xtol is too small. no further improvement in
    //               the approximate solution x is possible.

    //     info = 8  gtol is too small. fvec is orthogonal to the
    //               columns of the jacobian to machine precision.

    //   nfev is an integer output variable set to the number of
    //     calls to fcn.

    //   fjac is an output m by n array. the upper n by n submatrix
    //     of fjac contains an upper triangular matrix r with
    //     diagonal elements of nonincreasing magnitude such that

    //            t     t           t
    //           p *(jac *jac)*p = r *r,

    //     where p is a permutation matrix and jac is the final
    //     calculated jacobian. column j of p is column ipvt(j)
    //     (see below) of the identity matrix. the lower trapezoidal
    //     part of fjac contains information generated during
    //     the computation of r.

    //   ldfjac is a positive integer input variable not less than m
    //     which specifies the leading dimension of the array fjac.

    //   ipvt is an integer output array of length n. ipvt
    //     defines a permutation matrix p such that jac*p = q*r,
    //     where jac is the final calculated jacobian, q is
    //     orthogonal (not stored), and r is upper triangular
    //     with diagonal elements of nonincreasing magnitude.
    //     column j of p is column ipvt(j) of the identity matrix.

    //   qtf is an output array of length n which contains
    //     the first n elements of the vector (q transpose)*fvec.

    //   wa1, wa2, and wa3 are work arrays of length n.

    //   wa4 is a work array of length m.

    // subprograms called

    //   user-supplied ...... fcn

    //   minpack-supplied ... dpmpar,enorm,fdjac2,lmpar,qrfac

    //   fortran-supplied ... dabs,dmax1,dmin1,dsqrt,mod

    // argonne national laboratory. minpack project. march 1980.
    // burton s. garbow, kenneth e. hillstrom, jorge j. more

    // Initialized data

    const double one   = 1.0;
    const double p1    = 0.1;
    const double p5    = 0.5;
    const double p25   = 0.25;
    const double p75   = 0.75;
    const double p0001 = 1e-4;
    const double zero  = 0.0;

    // System generated locals

    // Local variables
    double par, sum;
    int iter;
    double temp, temp1, temp2;
    int iflag;
    double delta;
    double ratio;
    double fnorm, gnorm;
    double pnorm, xnorm, fnorm1, actred, dirder, prered;

    // Parameter adjustments

    --wa4;
    --fvec;
    --wa3;
    --wa2;
    --wa1;
    --qtf;
    --ipvt;
    --diag;
    --x;

    const int fjac_dim1 = ldfjac;
    const int fjac_offset = 1 + fjac_dim1;

    fjac -= fjac_offset;

    *info = 0;
    iflag = 0;
    *nfev = 0;

    // Check the input parameters for errors.

    if (n <= 0 || m < n || ldfjac < m || *ftol < zero || *xtol < zero || *gtol < zero || *maxfev <= 0 || *factor <= zero)
    {
        goto L300;
    }

    if (*mode != 2)
    {
        goto L20;
    }

    for (int j = 1; j <= n; ++j)
    {
        if (diag[j] <= zero)
        {
            goto L300;
        }
    }

L20:

    // Evaluate the function at the starting point
    // and calculate its norm.

    iflag = 1;
    (*fcn)(m, n, &x[1], &fvec[1], &iflag);
    *nfev = 1;

    if (iflag < 0)
    {
        goto L300;
    }

    fnorm = enorm(m, &fvec[1]);

    // Initialize levenberg-marquardt parameter and iteration counter.

    par = zero;
    iter = 1;

    // Beginning of the outer loop.

L30:

    // Calculate the jacobian matrix.

    iflag = 2;
    fdjac2(fcn, m, n, &x[1], &fvec[1], &fjac[fjac_offset], ldfjac, &iflag, *epsfcn, &wa4[1]);
    *nfev += n;

    if (iflag < 0)
    {
        goto L300;
    }

    // If requested, call fcn to enable printing of iterates.

    if (*nprint <= 0)
    {
        goto L40;
    }

    iflag = 0;
    if ((iter - 1) % *nprint == 0)
    {
        (*fcn)(m, n, &x[1], &fvec[1], &iflag);
    }

    if (iflag < 0)
    {
        goto L300;
    }

L40:

    // Compute the qr factorization of the Jacobian.

    qrfac(m, n, &fjac[fjac_offset], ldfjac, true, &ipvt[1], n, &wa1[1], &wa2[1], &wa3[1]);

    // On the first iteration and if mode is 1, scale according
    // to the norms of the columns of the initial Jacobian.

    if (iter != 1)
    {
        goto L80;
    }

    if (*mode == 2)
    {
        goto L60;
    }

    for (int j = 1; j <= n; ++j)
    {
        diag[j] = wa2[j];
        if (wa2[j] == zero)
        {
            diag[j] = one;
        }
    }

L60:

    // On the first iteration, calculate the norm of the scaled x
    // and initialize the step bound delta.

    for (int j = 1; j <= n; ++j)
    {
        wa3[j] = diag[j] * x[j];
    }

    xnorm = enorm(n, &wa3[1]);
    delta = *factor * xnorm;

    if (delta == zero)
    {
        delta = *factor;
    }

L80:

    // form (q transpose)*fvec and store the first n components in
    // qtf.

    for (int i = 1; i <= m; ++i)
    {
        wa4[i] = fvec[i];
    }

    for (int j = 1; j <= n; ++j)
    {
        if (fjac[j + j * fjac_dim1] == zero)
        {
            goto L120;
        }

        sum = zero;

        for (int i = j; i <= m; ++i)
        {
            sum += fjac[i + j * fjac_dim1] * wa4[i];
        }

        temp = -sum / fjac[j + j * fjac_dim1];

        for (int i = j; i <= m; ++i)
        {
            wa4[i] += fjac[i + j * fjac_dim1] * temp;
        }
L120:
        fjac[j + j * fjac_dim1] = wa1[j];
        qtf[j] = wa4[j];
    }

    // Compute the norm of the scaled gradient.

    gnorm = zero;
    if (fnorm == zero)
    {
        goto L170;
    }

    for (int j = 1; j <= n; ++j)
    {
        const int l = ipvt[j];
        if (wa2[l] == zero)
        {
            goto L150;
        }
        sum = zero;

        for (int i = 1; i <= j; ++i)
        {
            sum += fjac[i + j * fjac_dim1] * (qtf[i] / fnorm);
        }
        gnorm = fmax(gnorm, fabs(sum / wa2[l]));
L150:
        ;
    }

L170:

    // Test for convergence of the gradient norm.

    if (gnorm <= *gtol)
    {
        *info = 4;
    }

    if (*info != 0)
    {
        goto L300;
    }

    // Rescale if necessary.

    if (*mode == 2)
    {
        goto L190;
    }

    for (int j = 1; j <= n; ++j)
    {
        diag[j] = fmax(diag[j], wa2[j]);
    }

L190:

    // Beginning of the inner loop.

L200:

    // Determine the Levenberg-Marquardt parameter.

    lmpar(n, &fjac[fjac_offset], ldfjac, &ipvt[1], &diag[1], &qtf[1], delta, &par, &wa1[1], &wa2[1], &wa3[1], &wa4[1]);

    // Store the direction p and x + p. calculate the norm of p.

    for (int j = 1; j <= n; ++j)
    {
        wa1[j] = -wa1[j];
        wa2[j] = x[j] + wa1[j];
        wa3[j] = diag[j] * wa1[j];
    }

    pnorm = enorm(n, &wa3[1]);

    // On the first iteration, adjust the initial step bound.

    if (iter == 1)
    {
        delta = fmin(delta, pnorm);
    }

    // Evaluate the function at x + p and calculate its norm.

    iflag = 1;
    (*fcn)(m, n, &wa2[1], &wa4[1], &iflag);
    ++(*nfev);

    if (iflag < 0)
    {
        goto L300;
    }

    fnorm1 = enorm(m, &wa4[1]);

    // Compute the scaled actual reduction.

    actred = -one;
    if (p1 * fnorm1 < fnorm)
    {
        actred = one - square(fnorm1 / fnorm);
    }

    // Compute the scaled predicted reduction and
    // The scaled directional derivative.

    for (int j = 1; j <= n; ++j)
    {
        wa3[j] = zero;
        const int l = ipvt[j];
        temp = wa1[l];

        for (int i = 1; i <= j; ++i)
        {
            wa3[i] += fjac[i + j * fjac_dim1] * temp;
        }
    }

    temp1 = enorm(n, &wa3[1]) / fnorm;
    temp2 = sqrt(par) * pnorm / fnorm;
    prered = square(temp1) + square(temp2) / p5;
    dirder = -(square(temp1) + square(temp2));

    // Compute the ratio of the actual to the predicted
    // Reduction.

    ratio = zero;
    if (prered != zero)
    {
        ratio = actred / prered;
    }

    // Update the step bound.

    if (ratio > p25)
    {
        goto L240;
    }

    if (actred >= zero)
    {
        temp = p5;
    }

    if (actred < zero)
    {
        temp = p5 * dirder / (dirder + p5 * actred);
    }

    if (p1 * fnorm1 >= fnorm || temp < p1)
    {
        temp = p1;
    }

    delta = temp * fmin(delta, pnorm / p1);
    par /= temp;
    goto L260;

L240:

    if (par != zero && ratio < p75)
    {
        goto L250;
    }

    delta = pnorm / p5;
    par = p5 * par;

L250:
L260:

    // Test for successful iteration.

    if (ratio < p0001)
    {
        goto L290;
    }

    // Successful iteration. update x, fvec, and their norms.

    for (int j = 1; j <= n; ++j)
    {
        x[j] = wa2[j];
        wa2[j] = diag[j] * x[j];
    }

    for (int i = 1; i <= m; ++i)
    {
        fvec[i] = wa4[i];
    }

    xnorm = enorm(n, &wa2[1]);
    fnorm = fnorm1;
    ++iter;

L290:

    // Tests for convergence.

    if (fabs(actred) <= *ftol && prered <= *ftol && p5 * ratio <= one)
    {
        *info = 1;
    }

    if (delta <= *xtol * xnorm)
    {
        *info = 2;
    }

    if (fabs(actred) <= *ftol && prered <= *ftol && p5 * ratio <= one && *info == 2)
    {
        *info = 3;
    }

    if (*info != 0)
    {
        goto L300;
    }

    // Tests for termination and stringent tolerances.

    if (*nfev >= *maxfev)
    {
        *info = 5;
    }

    if (fabs(actred) <= MACHINE_EPSILON && prered <= MACHINE_EPSILON && p5 * ratio <= one)
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
        goto L300;
    }

    // End of the inner loop. Repeat if iteration unsuccessful.

    if (ratio < p0001)
    {
        goto L200;
    }

    // End of the outer loop.

    goto L30;

L300:

    // Termination, either normal or user imposed.

    if (iflag < 0)
    {
        *info = iflag;
    }

    iflag = 0;
    if (*nprint > 0)
    {
        (*fcn)(m, n, &x[1], &fvec[1], &iflag);
    }
}
