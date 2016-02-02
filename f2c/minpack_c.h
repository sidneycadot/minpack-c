
// minpack_c.h

#ifndef minpack_c_h
#define minpack_c_h

#include <stdbool.h>

typedef int integer;
typedef double doublereal;

#include "dpmpar.h"

typedef void (*lmder_fcn)(const integer m, const integer n, doublereal *, doublereal *, doublereal *, const integer, integer *);
typedef void (*lmdif_fcn)(const integer m, const integer n, doublereal * x, doublereal * fvec, integer * iflag);

doublereal enorm(
        const integer n,
        const doublereal * x
    );

void qrsolv(
        const integer  n,
        doublereal    *r,
        const integer  ldr,
        const integer *ipvt,
        doublereal    *diag,
        doublereal    *qtb,
        doublereal    *x,
        doublereal    *sdiag,
        doublereal    *wa
    );

void fdjac2(
        lmdif_fcn         fcn,
        const integer     m,
        const integer     n,
        doublereal       *x,
        const doublereal *fvec,
        doublereal       *fjac,
        const integer     ldfjac,
        integer          *iflag,
        const doublereal  epsfcn,
        doublereal       *wa
    );

void lmpar(
        const integer     n,
        doublereal       *r,
        const integer     ldr,
        integer          *ipvt,
        doublereal       *diag,
        doublereal       *qtb,
        doublereal       *delta,
        doublereal       *par,
        doublereal       *x,
        doublereal       *sdiag,
        doublereal       *wa1,
        doublereal       *wa2
    );

void qrfac(
        const integer     m,
        const integer     n,
        doublereal       *a,
        const integer     lda,
        const bool        pivot,
        integer          *ipvt,
        const integer     lipvt,
        doublereal       *rdiag,
        doublereal       *acnorm,
        doublereal       *wa
    );

void lmder(
        lmder_fcn         fcn,
        const integer     m,
        const integer     n,
        doublereal       *x,
        doublereal       *fvec,
        doublereal       *fjac,
        const integer     ldfjac,
        const doublereal  ftol,
        const doublereal  xtol,
        const doublereal  gtol,
        const integer     maxfev,
        doublereal       *diag,
        const integer     mode,
        const doublereal  factor,
        const integer     nprint,
        integer          *info,
        integer          *nfev,
        integer          *njev,
        integer          *ipvt,
        doublereal       *qtf,
        doublereal       *wa1,
        doublereal       *wa2,
        doublereal       *wa3,
        doublereal       *wa4
    );

void lmder1(
        lmder_fcn         fcn,
        const integer     m,
        const integer     n,
        doublereal       *x,
        doublereal       *fvec,
        doublereal       *fjac,
        const integer     ldfjac,
        const doublereal  tol,
        integer          *info,
        integer          *ipvt,
        doublereal       *wa,
        const integer     lwa
    );

void lmdif(
        lmdif_fcn         fcn,
        const integer     m,
        const integer     n,
        doublereal       *x,
        doublereal       *fvec,
        const doublereal  ftol,
        const doublereal  xtol,
        const doublereal  gtol,
        const integer     maxfev,
        const doublereal  epsfcn,
        doublereal       *diag,
        const integer     mode,
        const doublereal  factor,
        const integer     nprint,
        integer          *info,
        integer          *nfev,
        doublereal       *fjac,
        const integer     ldfjac,
        integer          *ipvt,
        doublereal       *qtf,
        doublereal       *wa1,
        doublereal       *wa2,
        doublereal       *wa3,
        doublereal       *wa4
    );

void lmdif1(
        lmdif_fcn         fcn,
        const integer     m,
        const integer     n,
        doublereal       *x,
        doublereal       *fvec,
        const doublereal  tol,
        integer          *info,
        integer          *iwa,
        doublereal       *wa,
        const integer     lwa
    );

#endif // minpack_c_h
