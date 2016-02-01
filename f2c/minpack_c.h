
// minpack_c.h

#ifndef minpack_c_h
#define minpack_c_h

#include <stdbool.h>

typedef int integer;
typedef double doublereal;

typedef void (*S_fp)(const integer m, const integer n, doublereal *, doublereal *, doublereal *, const integer, integer *);
//typedef int (*U_fp)(int);

typedef void (*lmdif_fcn)(const integer m, const integer n, doublereal * x, doublereal * fvec, integer * iflag);

doublereal dpmpar(integer *i__);

doublereal enorm(const integer n, doublereal *x);

double fmin(double a, double b);
double fmax(double a, double b);

void qrsolv(
        const integer n,
        doublereal *r__,
        const integer ldr,
        integer    *ipvt,
        doublereal *diag,
        doublereal *qtb,
        doublereal *x,
        doublereal *sdiag,
        doublereal *wa
    );

int fdjac2(
        lmdif_fcn  fcn,
        const integer m,
        const integer n,
        doublereal *x,
        doublereal *fvec,
        doublereal *fjac,
        const integer ldfjac,
        integer    *iflag, 
        doublereal *epsfcn,
        doublereal *wa
    );

int lmpar(
        const integer n,
        doublereal *r__,
        const integer ldr,
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
        const integer m,
        const integer n,
        doublereal *a,
        const integer lda,
        bool       *pivot,
        integer    *ipvt,
        const integer lipvt,
        doublereal *rdiag,
        doublereal *acnorm,
        doublereal *wa
    );

int lmder(
        S_fp       fcn,
        const integer m,
        const integer n,
        doublereal *x,
        doublereal *fvec,
        doublereal *fjac,
        const integer ldfjac,
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
        const integer m,
        const integer n,
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

int lmdif(
        lmdif_fcn   fcn,
        const integer m,
        const integer n,
        doublereal *x, 
        doublereal *fvec,
        doublereal *ftol,
        doublereal *xtol,
        doublereal *gtol,
        integer    *maxfev,
        doublereal *epsfcn,
        doublereal *diag,
        integer    * mode,
        doublereal *factor,
        integer    *nprint,
        integer    *info,
        integer    *nfev,
        doublereal *fjac,
        const integer ldfjac,
        integer    *ipvt,
        doublereal *qtf,
        doublereal *wa1,
        doublereal *wa2,
        doublereal *wa3,
        doublereal *wa4
    );

int lmdif1(
        lmdif_fcn   fcn,
        const integer m,
        const integer n,
        doublereal *x,
        doublereal *fvec,
        doublereal *tol,
        integer    *info,
        integer    *iwa,
        doublereal *wa,
        integer    *lwa
    );

#endif // minpack_c_h
