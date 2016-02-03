
// minpack_c.h

#ifndef minpack_c_h
#define minpack_c_h

#include <stdbool.h>

#define PIVOT_OFFSET 1

typedef void (*lmder_fcn)(const int m, const int n, double * x, double * fvec, double * fjac, const int ldfjac, int * iflag);
typedef void (*lmdif_fcn)(const int m, const int n, double * x, double * fvec,                                  int * iflag);

inline double square(const double x)
{
    return x * x;
}

double enorm(
        const int     n,
        const double *x
    );

void qrsolv(
        const int     n,
        double       *r,
        const int     ldr,
        const int    *ipvt,
        double       *diag,
        double       *qtb,
        double       *x,
        double       *sdiag,
        double       *wa
    );

void fdjac2(
        lmdif_fcn     fcn,
        const int     m,
        const int     n,
        double       *x,
        const double *fvec,
        double       *fjac,
        const int     ldfjac,
        int          *iflag,
        const double  epsfcn,
        double       *wa
    );

void lmpar(
        const int     n,
        double       *r,
        const int     ldr,
        int          *ipvt,
        double       *diag,
        double       *qtb,
        double       *delta,
        double       *par,
        double       *x,
        double       *sdiag,
        double       *wa1,
        double       *wa2
    );

void qrfac(
        const int     m,
        const int     n,
        double       *a,
        const int     lda,
        const bool    pivot,
        int          *ipvt,
        const int     lipvt,
        double       *rdiag,
        double       *acnorm,
        double       *wa
    );

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
    );

void lmder1(
        lmder_fcn     fcn,
        const int     m,
        const int     n,
        double       *x,
        double       *fvec,
        double       *fjac,
        const int     ldfjac,
        const double  tol,
        int          *info,
        int          *ipvt,
        double       *wa,
        const int     lwa
    );

void lmdif(
        lmdif_fcn     fcn,
        const int     m,
        const int     n,
        double       *x,
        double       *fvec,
        const double  ftol,
        const double  xtol,
        const double  gtol,
        const int     maxfev,
        const double  epsfcn,
        double       *diag,
        const int     mode,
        const double  factor,
        const int     nprint,
        int          *info,
        int          *nfev,
        double       *fjac,
        const int     ldfjac,
        int          *ipvt,
        double       *qtf,
        double       *wa1,
        double       *wa2,
        double       *wa3,
        double       *wa4
    );

void lmdif1(
        lmdif_fcn     fcn,
        const int     m,
        const int     n,
        double       *x,
        double       *fvec,
        const double  tol,
        int          *info,
        int          *iwa,
        double       *wa,
        const int     lwa
    );

#endif // minpack_c_h
