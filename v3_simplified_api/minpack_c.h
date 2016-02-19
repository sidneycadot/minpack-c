
/////////////////
// minpack_c.h //
/////////////////

#ifndef minpack_c_h
#define minpack_c_h

#include <math.h>
#include <stdbool.h>

enum lm_info_t
{
    LM_INFO_0 = 0, // Improper input parameters.
    LM_INFO_1 = 1, // Both actual and predicted relative reductions in the sum of squares are at most ftol.
    LM_INFO_2 = 2, // Relative error between two consecutive iterates is at most xtol.
    LM_INFO_3 = 3, // Conditions for info = 1 and info = 2 both hold
    LM_INFO_4 = 4, // The cosine of the angle between fvec and any column of the Jacobian is at most gtol in absolute value.
    LM_INFO_5 = 5, // Number of calls to fcn has reached or exceeded 200*(n+1).
    LM_INFO_6 = 6, // ftol is too small. No further reduction in the sum of squares is possible.
    LM_INFO_7 = 7, // xtol is too small. No further improvement in the approximate solution x is possible.
    LM_INFO_8 = 8  // gtol is too small. fvec is orthogonal to the columns of the Jacobian to machine precision.
};

enum lm_mode_t
{
    LM_MODE_1 = 1,
    LM_MODE_2 = 2
};

enum lm_iflag_t
{
    LM_IFLAG_1 = 1, // calculate function value
    LM_IFLAG_2 = 2  // calculate forward difference for Jacobian
};

typedef void (*lmdif_fcn_ptr)(const int m, const int n, double *x, double *fvec, const enum lm_iflag_t iflag, void * fcn_context);

inline double square(const double x)
{
    return x * x;
}

static const double MACHINE_EPSILON = 2.22044604926e-16;
static const double MACHINE_MINPOS  = 2.22507385852e-308;

double enorm(const int n, const double *x);

void fdjac2(
        lmdif_fcn_ptr fcn, void * fcn_context,
        const int m, const int n,
        double *x, const double *fvec, double *fjac, const int ldfjac, const double epsfcn, double *wa);

void qrfac(
        const int m, const int n,
        double *a, const int lda, int *ipvt, double *rdiag, double *acnorm, double *wa);

void qrsolv(
        const int n,
        double *r, const int ldr,
        const int *ipvt, const double *diag, const double *qtb, double *x, double *sdiag, double *wa);

void lmpar(
        const int n,
        double *r, const int ldr,
        const int *ipvt, const double *diag, const double *qtb,
        const double delta, double *par,
        double *x, double *sdiag, double *wa1, double *wa2);

enum lm_info_t lmdif(
        lmdif_fcn_ptr fcn, void *fcn_context,
        const int m, const int n,
        double *x, double *fvec,
        const double ftol, const double xtol, const double gtol,
        const int maxfev, const double epsfcn, double *diag, const enum lm_mode_t mode, const double factor, int *nfev,
        double *fjac, const int ldfjac, int *ipvt, double *qtf,
        double *wa1, double *wa2, double *wa3, double *wa4);

enum lm_info_t lmdif1(
        lmdif_fcn_ptr fcn, void *fcn_context,
        const int m, const int n,
        double *x, double *fvec,
        const double tol, int *iwa, double *wa, const int lwa);

#endif // minpack_c_h
