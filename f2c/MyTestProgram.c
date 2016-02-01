
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "minpack_c.h"

double goodfunc(double x)
{
    x = (x - 11.234);
    x /= 2.345;
    return exp(-x * x);
}

double fitfunc(double x, double offset, double scale)
{
    x = (x - offset);
    x /= scale;
    return exp(-x * x);
}

double PX[20];
double PY[20];

/*       fcn is the name of the user-supplied subroutine which */
/*         calculates the functions and the jacobian. fcn must */
/*         be declared in an external statement in the user */
/*         calling program, and should be written as follows. */

/*         subroutine fcn(m,n,x,fvec,fjac,ldfjac,iflag) */
/*         integer m,n,ldfjac,iflag */
/*         double precision x(n),fvec(m),fjac(ldfjac,n) */
/*         ---------- */
/*         if iflag = 1 calculate the functions at x and */
/*         return this vector in fvec. do not alter fjac. */
/*         if iflag = 2 calculate the jacobian at x and */
/*         return this matrix in fjac. do not alter fvec. */
/*         ---------- */
/*         return */
/*         end */

/*         the value of iflag should not be changed by fcn unless */
/*         the user wants to terminate execution of lmder1. */
/*         in this case set iflag to a negative integer. */

void fcn(integer *m, integer *n, doublereal * x, doublereal *fvec, doublereal *fjac, integer *ldfjac, integer *iflag)
{
    printf("invocation fcn(m = %d, n = %d, x[], fvec[], fjac[], ldfac = %d, iflag = %d)\n", *m, *n, (*ldfjac), (*iflag));
    for (int i = 0; i <= (*n); ++i)
    {
        printf("x[%d] = %f\n", i, x[i]);
    }

    if (iflag == 1)
    {
        for (int i = 0; i < (*m); ++i)
        {
            fvec[i] = fitfunc(PX[i], x[0], x[1]) - PY[i];
        }
    }

    exit(1);
}

/*       m is a positive integer input variable set to the number */
/*         of functions. */

/*       n is a positive integer input variable set to the number */
/*         of variables. n must not exceed m. */

/*       x is an array of length n. on input x must contain */
/*         an initial estimate of the solution vector. on output x */
/*         contains the final estimate of the solution vector. */

/*       fvec is an output array of length m which contains */
/*         the functions evaluated at the output x. */

/*       fjac is an output m by n array. the upper n by n submatrix */
/*         of fjac contains an upper triangular matrix r with */
/*         diagonal elements of nonincreasing magnitude such that */

/*                t     t           t */
/*               p *(jac *jac)*p = r *r, */

/*         where p is a permutation matrix and jac is the final */
/*         calculated jacobian. column j of p is column ipvt(j) */
/*         (see below) of the identity matrix. the lower trapezoidal */
/*         part of fjac contains information generated during */
/*         the computation of r. */

/*       ldfjac is a positive integer input variable not less than m */
/*         which specifies the leading dimension of the array fjac. */

/*       tol is a nonnegative input variable. termination occurs */
/*         when the algorithm estimates either that the relative */
/*         error in the sum of squares is at most tol or that */
/*         the relative error between x and the solution is at */
/*         most tol. */

/*       info is an integer output variable. if the user has */
/*         terminated execution, info is set to the (negative) */
/*         value of iflag. see description of fcn. otherwise, */
/*         info is set as follows. */

/*         info = 0  improper input parameters. */

/*         info = 1  algorithm estimates that the relative error */
/*                   in the sum of squares is at most tol. */

/*         info = 2  algorithm estimates that the relative error */
/*                   between x and the solution is at most tol. */

/*         info = 3  conditions for info = 1 and info = 2 both hold. */

/*         info = 4  fvec is orthogonal to the columns of the */
/*                   jacobian to machine precision. */

/*         info = 5  number of calls to fcn with iflag = 1 has */
/*                   reached 100*(n+1). */

/*         info = 6  tol is too small. no further reduction in */
/*                   the sum of squares is possible. */

/*         info = 7  tol is too small. no further improvement in */
/*                   the approximate solution x is possible. */

/*       ipvt is an integer output array of length n. ipvt */
/*         defines a permutation matrix p such that jac*p = q*r, */
/*         where jac is the final calculated jacobian, q is */
/*         orthogonal (not stored), and r is upper triangular */
/*         with diagonal elements of nonincreasing magnitude. */
/*         column j of p is column ipvt(j) of the identity matrix. */

/*       wa is a work array of length lwa. */

/*       lwa is a positive integer input variable not less than 5*n+m. */

int main(void)
{
    printf("hello!\n");

    integer m = 20; // number of points?
    integer n = 2; // number of parameters to solve for (2).

    doublereal x[3];           // initial parameter guess and solution
    doublereal fvec[1001];
    doublereal fjac[1001];
    integer    ldfjac = m;
    doublereal tol = 1e-6;
    integer    info;
    integer    ipvt[1001];
    doublereal wa[10001];
    integer    lwa = 5000; // length of wa buffer

    x[0] = 101;
    x[1] = 102;
    x[2] = 103;

    for (int i = 1; i <= 20; ++i)
    {
        PX[i] = 0.3 + 1.1 * i;
        PY[i] = goodfunc(PX[i]);
    }

    int rc = lmder1(
        fcn,
        &m,
        &n,
        x,
        fvec,
        fjac,
        &ldfjac,
        &tol,
        &info,
        ipvt,
        wa,
        &lwa
    );

    printf("rc: %d\n", rc);
    printf("info: %d\n", info);

    return 0;
}
