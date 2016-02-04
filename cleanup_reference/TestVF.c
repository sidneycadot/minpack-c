
//////////////
// TestVF.c //
//////////////

#include <math.h>
#include <assert.h>
#include <stdio.h>
#include "my_include.h"

double v_arr[3001];
double f_arr[3001];

double f_min;
double f_max;

void fcn(integer *m, integer *n, doublereal *x, doublereal *fvec, integer *iflag)
{
    const double PI = 3.141592653589793238462643383;

    printf("invocation fcn(m = %d, n = %d, x[], fvec[], iflag = %d)\n", *m, *n, *iflag);

    for (int i = 0; i < *n; ++i)
    {
        printf("x[%d] = %f\n", i, x[i]);
    }

    const double logFrequencySlackLeft  = x[0];
    const double logFrequencySlackRight = x[1];
    const double curviness              = x[2];
    const double vOffset                = x[3];
    const double vScale                 = x[4];

    const double tanh_curviness = tanh(curviness);

    for (int i = 0; i < *m; ++i)
    {
        const double f = f_arr[i];

        const double f1 = f_min - exp(logFrequencySlackLeft);
        const double f2 = f_max + exp(logFrequencySlackRight);

        const double x  = (f - f1) / (f2 - f1) * 2 - 1;
        const double xx = x + tanh_curviness * (1 - x) * (1 + x) / 2;
        const double y  = tan((PI / 2) * xx);
        const double v  = vOffset + vScale * y;

        fvec[i] = (v - v_arr[i]);
    }
}

doublereal wa[30000];
integer lwa = 30000;

int main(void)
{
    integer m = 3001; // number of functions (i.e., data points)
    integer n = 5; // number of parameters;

    FILE * f = fopen("vf.txt", "r");

    for (int i = 0; i < m; ++i)
    {
        int r = fscanf(f, "%lf%lf", &v_arr[i], &f_arr[i]);
        assert(r == 2);
        //printf("%f %f\n", v_arr[i], f_arr[i]);
    }
    fclose(f);

    f_min = +INFINITY;
    f_max = -INFINITY;

    for (int i = 0; i < m; ++i)
    {
        f_min = fmin(f_min, f_arr[i]);
        f_max = fmax(f_max, f_arr[i]);
    }

    printf("frequency range: %f .. %f\n", f_min, f_max);

    double x[5]; // parameters: (logFrequencySlackLeft, logFrequencySlackRight, curviness, vOffset, vScale, f_min, f_max);

    double fvec[3001];

    int iwa[5];

    x[0] = 20.3;
    x[1] = 20.0;
    x[2] = 1.0;
    x[3] = 1.5;
    x[4] = 0.40;

    double tol = 1e-12;
    int info;

    lmdif1_(fcn, &m, &n, x, fvec, &tol, &info, iwa, wa, &lwa);

    printf("info: %d\n", info);

    if (1)
    {
        printf("========================= x[] after lmdif1: {");
        for (int i = 0; i < n; ++i)
        {
            if (i != 0)
            {
                printf(", ");
            }
            printf("%.20g", x[i]);
        }
        printf("}\n");
    }

    printf("reference solution :                        {19.802463007970860787, 21.810082992538220026, 5.9223404727777602119, 1.5958684488737919782, 0.3859655857639207599}\n");

    return 0;
}
