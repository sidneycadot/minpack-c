
#ifndef my_include_h
#define my_include_h

#include <math.h>
#include <stdbool.h>

typedef bool logical;
typedef double doublereal;
typedef int integer;

typedef void (*U_fp)(integer * m, integer *n, doublereal *x, doublereal *fvec, int *iflag);

int lmdif1_(U_fp fcn, integer *m, integer *n, doublereal *x,
	    doublereal *fvec, doublereal *tol, integer *info, integer *iwa,
	    doublereal *wa, integer *lwa);
  
#endif // my_include_h
