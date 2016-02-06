
// fdjac2.f

#include "minpack_c.h"

void fdjac2(lmdif_fcn_ptr fcn, void * fcn_context,
            const int m, const int n,
            double *x, const double *fvec,
            double *fjac, const int ldfjac,
            const double epsfcn, double *wa)
{
    // This subroutine computes a forward-difference approximation
    // to the m by n Jacobian matrix associated with a specified
    // problem of m functions in n variables.

    // The subroutine statement is

    //   fdjac2(fcn, m, n, x, fvec, fjac, ldfjac, iflag, epsfcn, wa)

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
    //     the user wants to terminate execution of fdjac2.
    //     In this case set iflag to a negative integer.

    //   m is a positive integer input variable set to the number
    //     of functions.

    //   n is a positive integer input variable set to the number
    //     of variables. n must not exceed m.

    //   x is an input array of length n.

    //   fvec is an input array of length m which must contain the
    //     functions evaluated at x.

    //   fjac is an output m by n array which contains the
    //     approximation to the Jacobian matrix evaluated at x.

    //   ldfjac is a positive integer input variable not less than m
    //     which specifies the leading dimension of the array fjac.

    //   epsfcn is an input variable used in determining a suitable
    //     step length for the forward-difference approximation. This
    //     approximation assumes that the relative errors in the
    //     functions are of the order of epsfcn. If epsfcn is less
    //     than the machine precision, it is assumed that the relative
    //     errors in the functions are of the order of the machine
    //     precision.

    //   wa is a work array of length m.

    // Subprograms called:

    //   user-supplied ...... fcn

    // Argonne National Laboratory. MINPACK project. March 1980.
    // Burton S. Garbow, Kenneth E. Hillstrom, Jorge J. More

    const double eps = sqrt(fmax(epsfcn, MACHINE_EPSILON));

    for (int j = 0; j < n; ++j)
    {
        const double temp = x[j];

        double h = eps * fabs(temp);

        if (h == 0.0)
        {
            h = eps;
        }

        x[j] = temp + h;

        (*fcn)(m, n, x, wa, LM_IFLAG_2, fcn_context);

        x[j] = temp;

        for (int i = 0; i < m; ++i)
        {
            fjac[i + j * ldfjac] = (wa[i] - fvec[i]) / h;
        }
    }
}
