
/////////////////
// minpack_c.h //
/////////////////

#ifndef minpack_c_h
#define minpack_c_h

#include <math.h>
#include <stdbool.h>

typedef void (*lmdif_fcn_ptr)(const int m, const int n, double *x, double *fvec, int *iflag);

inline double square(const double x)
{
    return x * x;
}

static const int PIVOT_OFFSET = 1;

static const double MACHINE_EPSILON = 2.22044604926e-16;
static const double MACHINE_MINPOS  = 2.22507385852e-308;

double enorm(const int n, const double *x);
void fdjac2(lmdif_fcn_ptr fcn, const int m, const int n, double *x, const double *fvec, double *fjac, const int ldfjac, int *iflag, const double epsfcn, double *wa);
void qrfac(const int m, const int n, double *a, const int lda, const bool pivot, int *ipvt, const int lipvt, double *rdiag, double *acnorm, double *wa);
void qrsolv(const int n, double *r, const int ldr, const int *ipvt, const double *diag, const double *qtb, double *x, double *sdiag, double *wa);
void lmpar(const int n, double *r, const int ldr, const int *ipvt, const double *diag, double *qtb, const double delta, double *par, double *x, double *sdiag, double *wa1, double *wa2);

void lmdif(lmdif_fcn_ptr fcn, const int m, const int n, double *x, double *fvec, const double ftol, const double xtol, const double gtol, const int maxfev,
           const double epsfcn, double *diag, const int mode, const double factor, const int nprint, int *info, int *nfev,
           double *fjac, const int ldfjac, int *ipvt, double *qtf, double *wa1, double *wa2, double *wa3, double *wa4);

void lmdif1(lmdif_fcn_ptr fcn, const int m, const int n, double *x, double *fvec, const double tol, int *info, int *iwa, double *wa, const int lwa);

#endif // minpack_c_h
