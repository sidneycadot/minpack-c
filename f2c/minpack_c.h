
// minpack_c.h

#ifndef minpack_c_h
#define minpack_c_h

#include <stdbool.h>

typedef int integer;
typedef double doublereal;

typedef void (*S_fp)(integer *m, integer *n, doublereal *, doublereal *, doublereal *, integer *, integer *);
typedef int (*U_fp)(int);

doublereal dpmpar(integer *i__);

doublereal enorm(integer *n, doublereal *x);

double fmin(double a, double b);
double fmax(double a, double b);

int lmpar(
        integer    *n,
        doublereal *r__,
        integer    *ldr,
        integer    *ipvt,
        doublereal *diag,
        doublereal *qtb,
        doublereal *delta,
        doublereal *par,
        doublereal *x,
        doublereal *sdiag,
        doublereal *wa1,
        doublereal *wa2
    );

int qrfac(
        integer    *m,
        integer    *n,
        doublereal *a,
        integer    *lda,
        bool       *pivot,
        integer    *ipvt,
        integer    *lipvt,
        doublereal *rdiag,
        doublereal *acnorm,
        doublereal *wa
    );

int lmder(
        S_fp       fcn,
        integer    *m,
        integer    *n,
        doublereal *x,
        doublereal *fvec,
        doublereal *fjac,
        integer    *ldfjac,
        doublereal *ftol,
        doublereal *xtol,
        doublereal *gtol,
        integer    *maxfev,
        doublereal *diag,
        integer    *mode,
        doublereal *factor,
        integer    *nprint,
        integer    *info,
        integer    *nfev,
        integer    *njev,
        integer    *ipvt,
        doublereal *qtf,
        doublereal *wa1,
        doublereal *wa2,
        doublereal *wa3,
        doublereal *wa4
    );

int lmder1(
        S_fp       fcn,
        integer    *m,
        integer    *n,
        doublereal *x,
        doublereal *fvec,
        doublereal *fjac,
        integer    *ldfjac,
        doublereal *tol,
        integer    *info,
        integer    *ipvt,
        doublereal *wa,
        integer    *lwa
    );

#endif // minpack_c_h
