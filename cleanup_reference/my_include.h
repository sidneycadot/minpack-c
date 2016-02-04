
#ifndef my_include_h
#define my_include_h

#include <math.h>
#include <stdbool.h>

typedef bool logical;
typedef double doublereal;
typedef int integer;

typedef void (*U_fp)(integer * m, integer *n, doublereal *x, doublereal *fvec, int *iflag);

inline double square(const double x)
{
    return x * x;
}

static const int PIVOT_OFFSET = 1;

static const double MACHINE_EPSILON =2.22044604926e-16;

double enorm(const int n, const double *x);
void fdjac2(U_fp fcn, int m, int n, double *x, double *fvec, double *fjac, const int ldfjac, int *iflag, const double epsfcn, double *wa);
void qrfac(const int m, const int n, double *a, const int lda, const bool pivot, int *ipvt, const int lipvt, double *rdiag, double *acnorm, double *wa);
void qrsolv(const int n, double *r, const int ldr, const int *ipvt, const double *diag, const double *qtb, double *x, double *sdiag, double *wa);

/* Subroutine */ int lmpar_(integer *n, doublereal *r__, integer *ldr, 
        integer *ipvt, doublereal *diag, doublereal *qtb, doublereal *delta, 
        doublereal *par, doublereal *x, doublereal *sdiag, doublereal *wa1, 
        doublereal *wa2);

int lmdif1_(U_fp fcn, integer *m, integer *n, doublereal *x,
        doublereal *fvec, doublereal *tol, integer *info, integer *iwa,
        doublereal *wa, integer *lwa);

#endif // my_include_h
