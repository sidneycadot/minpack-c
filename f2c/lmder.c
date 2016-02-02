// lmder.f -- translated by f2c (version 20100827).

#include <math.h>
#include <float.h>
#include <stdbool.h>
#include "minpack_c.h"

void lmder(
        lmder_fcn     fcn,
        const int     m,
        const int     n,
        double       *x,
        double       *fvec,
        double       *fjac,
        const int     ldfjac,
        const double  ftol,
        const double  xtol,
        const double  gtol,
        const int     maxfev,
        double       *diag,
        const int     mode,
        const double  factor,
        const int     nprint,
        int          *info,
        int          *nfev,
        int          *njev,
        int          *ipvt,
        double       *qtf,
        double       *wa1,
        double       *wa2,
        double       *wa3,
        double       *wa4
    )
{
    // The purpose of lmder is to minimize the sum of the squares of
    // m nonlinear functions in n variables by a modification of
    // the Levenberg-Marquardt algorithm. The user must provide a
    // subroutine which calculates the functions and the Jacobian.

    // The subroutine statement is

    //   subroutine lmder(fcn,m,n,x,fvec,fjac,ldfjac,ftol,xtol,gtol,
    //                    maxfev,diag,mode,factor,nprint,info,nfev,
    //                    njev,ipvt,qtf,wa1,wa2,wa3,wa4)

    // where

    //   fcn is the name of the user-supplied subroutine which
    //     calculates the functions and the Jacobian. fcn must
    //     be declared in an external statement in the user
    //     calling program, and should be written as follows.

    //     subroutine fcn(m,n,x,fvec,fjac,ldfjac,iflag)
    //     integer m,n,ldfjac,iflag
    //     double precision x(n),fvec(m),fjac(ldfjac,n)
    //     ----------
    //     if iflag = 1 calculate the functions at x and
    //     return this vector in fvec. do not alter fjac.
    //     if iflag = 2 calculate the Jacobian at x and
    //     return this matrix in fjac. do not alter fvec.
    //     ----------
    //     return
    //     end

    //     The value of iflag should not be changed by fcn unless
    //     the user wants to terminate execution of lmder.
    //     In this case set iflag to a negative integer.

    //   m is a positive integer input variable set to the number
    //     of functions.

    //   n is a positive integer input variable set to the number
    //     of variables. n must not exceed m.

    //   x is an array of length n. on input x must contain
    //     an initial estimate of the solution vector. On output x
    //     contains the final estimate of the solution vector.

    //   fvec is an output array of length m which contains
    //     the functions evaluated at the output x.

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

    //   maxfev is a positive integer input variable. termination
    //     occurs when the number of calls to fcn with iflag = 1
    //     has reached maxfev.

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
    //     interval (.1,100.).100. is a generally recommended value.

    //   nprint is an integer input variable that enables controlled
    //     printing of iterates if it is positive. in this case,
    //     fcn is called with iflag = 0 at the beginning of the first
    //     iteration and every nprint iterations thereafter and
    //     immediately prior to return, with x, fvec, and fjac
    //     available for printing. fvec and fjac should not be
    //     altered. if nprint is not positive, no special calls
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
    //               column of the Jacobian is at most gtol in
    //               absolute value.

    //     info = 5  number of calls to fcn with iflag = 1 has
    //               reached maxfev.

    //     info = 6  ftol is too small. no further reduction in
    //               the sum of squares is possible.

    //     info = 7  xtol is too small. no further improvement in
    //               the approximate solution x is possible.

    //     info = 8  gtol is too small. fvec is orthogonal to the
    //               columns of the Jacobian to machine precision.

    //   nfev is an integer output variable set to the number of
    //     calls to fcn with iflag = 1.

    //   njev is an integer output variable set to the number of
    //     calls to fcn with iflag = 2.

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

    // subprograms called

    //   user-supplied ...... fcn

    //   minpack-supplied ... dpmpar,enorm,lmpar,qrfac

    //   fortran-supplied ... dabs,dmax1,dmin1,dsqrt,mod

    // Argonne National Laboratory. MINPACK project. March 1980.
    // Burton S. Garbow, Kenneth E. Hillstrom, Jorge J. More

    // Initialized data

    const double one = 1.;
    const double p1 = .1;
    const double p5 = .5;
    const double p25 = .25;
    const double p75 = .75;
    const double p0001 = 1e-4;
    const double zero = 0.;

    // System generated locals
    int fjac_dim1, fjac_offset, i__1, i__2;
    double d__1, d__2, d__3;

    // Local variables

    int i__, j, l;
    double par, sum;
    int iter;
    double temp, temp1, temp2;
    int iflag;
    double delta;
    double ratio;
    double fnorm, gnorm, pnorm, xnorm = 0.0, fnorm1, actred, dirder, prered;

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
    fjac_dim1 = ldfjac;
    fjac_offset = 1 + fjac_dim1;
    fjac -= fjac_offset;

    // Function Body

    // epsmch is the machine precision.

    const double epsmch = DBL_EPSILON;

    *info = 0;
    iflag = 0;
    *nfev = 0;
    *njev = 0;

    // check the input parameters for errors.

    if (n <= 0 || m < n || ldfjac < m || ftol < zero || xtol < zero || gtol < zero || maxfev <= 0 || factor <= zero)
    {
        goto L300;
    }

    if (mode == 2)
    {

        i__1 = n;

        for (j = 1; j <= i__1; ++j)
        {
            if (diag[j] <= zero)
            {
                goto L300;
            }
        }
    }

    // evaluate the function at the starting point 
    // and calculate its norm.

    iflag = 1;
    (*fcn)(m, n, &x[1], &fvec[1], &fjac[fjac_offset], ldfjac, &iflag);
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

    // calculate the jacobian matrix.

    iflag = 2;
    (*fcn)(m, n, &x[1], &fvec[1], &fjac[fjac_offset], ldfjac, &iflag);
    ++(*njev);

    if (iflag < 0)
    {
        goto L300;
    }

    // if requested, call fcn to enable printing of iterates.

    if (nprint <= 0)
    {
        goto L40;
    }

    iflag = 0;
    if ((iter - 1) % nprint == 0)
    {
        (*fcn)(m, n, &x[1], &fvec[1], &fjac[fjac_offset], ldfjac, &iflag);
    }

    if (iflag < 0)
    {
        goto L300;
    }

L40:

    // Compute the qr factorization of the jacobian.

    qrfac(m, n, &fjac[fjac_offset], ldfjac, true, &ipvt[1], n, &wa1[1], &wa2[1], &wa3[1]);

    // on the first iteration and if mode is 1, scale according
    // to the norms of the columns of the initial jacobian.

    if (iter != 1)
    {
        goto L80;
    }

    if (mode == 2)
    {
        goto L60;
    }

    i__1 = n;
    for (j = 1; j <= i__1; ++j)
    {
        diag[j] = wa2[j];
        if (wa2[j] == zero)
        {
            diag[j] = one;
        }
    }

L60:

    // on the first iteration, calculate the norm of the scaled x
    // and initialize the step bound delta.

    i__1 = n;
    for (j = 1; j <= i__1; ++j)
    {
        wa3[j] = diag[j] * x[j];
    }
    xnorm = enorm(n, &wa3[1]);
    delta = factor * xnorm;
    if (delta == zero)
    {
        delta = factor;
    }

L80:

    // form (q transpose)*fvec and store the first n components in
    // qtf.

    i__1 = m;
    for (i__ = 1; i__ <= i__1; ++i__)
    {
        wa4[i__] = fvec[i__];
    }

    i__1 = n;
    for (j = 1; j <= i__1; ++j)
    {
        if (fjac[j + j * fjac_dim1] == zero)
        {
            goto L120;
        }
        sum = zero;
        i__2 = m;
        for (i__ = j; i__ <= i__2; ++i__)
        {
            sum += fjac[i__ + j * fjac_dim1] * wa4[i__];
        }
        temp = -sum / fjac[j + j * fjac_dim1];
        i__2 = m;
        for (i__ = j; i__ <= i__2; ++i__)
        {
            wa4[i__] += fjac[i__ + j * fjac_dim1] * temp;
        }

L120:
        fjac[j + j * fjac_dim1] = wa1[j];
        qtf[j] = wa4[j];
    }

    // compute the norm of the scaled gradient.

    gnorm = zero;

    if (fnorm == zero)
    {
        goto L170;
    }

    i__1 = n;
    for (j = 1; j <= i__1; ++j)
    {
        l = ipvt[j];
        if (wa2[l] == zero)
        {
            goto L150;
        }

        sum = zero;
        i__2 = j;
        for (i__ = 1; i__ <= i__2; ++i__)
        {
            sum += fjac[i__ + j * fjac_dim1] * (qtf[i__] / fnorm);
        }
        // Computing MAX
        d__2 = gnorm, d__3 = (d__1 = sum / wa2[l], fabs(d__1));
        gnorm = fmax(d__2, d__3);

L150: ;
    }

L170:

    // test for convergence of the gradient norm.

    if (gnorm <= gtol)
    {
        *info = 4;
    }

    if (*info != 0)
    {
        goto L300;
    }

    // rescale if necessary.

    if (mode == 2)
    {
        goto L190;
    }

    i__1 = n;
    for (j = 1; j <= i__1; ++j)
    {
        // Computing MAX
        d__1 = diag[j], d__2 = wa2[j];
        diag[j] = fmax(d__1,d__2);
    }

L190:

    // beginning of the inner loop.

L200:

    // determine the levenberg-marquardt parameter.

    lmpar(n, &fjac[fjac_offset], ldfjac, &ipvt[1], &diag[1], &qtf[1], &delta, &par, &wa1[1], &wa2[1], &wa3[1], &wa4[1]);

    // store the direction p and x + p. calculate the norm of p.

    i__1 = n;
    for (j = 1; j <= i__1; ++j)
    {
        wa1[j] = -wa1[j];
        wa2[j] = x[j] + wa1[j];
        wa3[j] = diag[j] * wa1[j];
    }
    pnorm = enorm(n, &wa3[1]);

    // on the first iteration, adjust the initial step bound.

    if (iter == 1)
    {
        delta = fmin(delta, pnorm);
    }

    // evaluate the function at x + p and calculate its norm.

    iflag = 1;
    (*fcn)(m, n, &wa2[1], &wa4[1], &fjac[fjac_offset], ldfjac, &iflag);
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
        // Computing 2nd power
        d__1 = fnorm1 / fnorm;
        actred = one - d__1 * d__1;
    }

    // Compute the scaled predicted reduction and
    // the scaled directional derivative.

    i__1 = n;
    for (j = 1; j <= i__1; ++j)
    {
        wa3[j] = zero;
        l = ipvt[j];
        temp = wa1[l];
        i__2 = j;
        for (i__ = 1; i__ <= i__2; ++i__)
        {
            wa3[i__] += fjac[i__ + j * fjac_dim1] * temp;
	}
    }

    temp1 = enorm(n, &wa3[1]) / fnorm;
    temp2 = sqrt(par) * pnorm / fnorm;

    // Computing 2nd power
    d__1 = temp1;

    // Computing 2nd power
    d__2 = temp2;

    prered = d__1 * d__1 + d__2 * d__2 / p5;

    // Computing 2nd power
    d__1 = temp1;

    // Computing 2nd power
    d__2 = temp2;

    dirder = -(d__1 * d__1 + d__2 * d__2);

    // Compute the ratio of the actual to the predicted
    // reduction.

    ratio = zero;
    if (prered != zero)
    {
        ratio = actred / prered;
    }

    // update the step bound.

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

    // Computing MIN

    d__1 = delta, d__2 = pnorm / p1;
    delta = temp * fmin(d__1, d__2);
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

    i__1 = n;
    for (j = 1; j <= i__1; ++j)
    {
        x[j] = wa2[j];
        wa2[j] = diag[j] * x[j];
    }

    i__1 = m;
    for (i__ = 1; i__ <= i__1; ++i__)
    {
        fvec[i__] = wa4[i__];
    }

    xnorm = enorm(n, &wa2[1]);
    fnorm = fnorm1;
    ++iter;

L290:

    // Tests for convergence.

    if (fabs(actred) <= ftol && prered <= ftol && p5 * ratio <= one)
    {
        *info = 1;
    }

    if (delta <= xtol * xnorm)
    {
        *info = 2;
    }

    if (fabs(actred) <= ftol && prered <= ftol && p5 * ratio <= one && *info == 2)
    {
        *info = 3;
    }

    if (*info != 0)
    {
        goto L300;
    }

    // Tests for termination and stringent tolerances.

    if (*nfev >= maxfev)
    {
        *info = 5;
    }

    if (fabs(actred) <= epsmch && prered <= epsmch && p5 * ratio <= one)
    {
        *info = 6;
    }

    if (delta <= epsmch * xnorm)
    {
        *info = 7;
    }

    if (gnorm <= epsmch)
    {
        *info = 8;
    }

    if (*info != 0)
    {
        goto L300;
    }

    // End of the inner loop. repeat if iteration unsuccessful.

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
    if (nprint > 0)
    {
        (*fcn)(m, n, &x[1], &fvec[1], &fjac[fjac_offset], ldfjac, &iflag);
    }
}
