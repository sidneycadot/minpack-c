
//////////////
// TestVF.c //
//////////////

#include <math.h>
#include <stdio.h>
#include "minpack_c.h"

double v_arr[3001];
double f_arr[3001];

void fcn(const int m, const int n, double *x, double *fvec, int *iflag)
{
    const double PI = 3.141592653589793238462643383;

    printf("invocation fcn(m = %d, n = %d, x[], fvec[], iflag = %d)\n", m, n, *iflag);

    for (int i = 0; i < 7; ++i)
    {
        printf("x[%d] = %f\n", i, x[i]);
    }

    const double logFrequencySlackLeft  = x[0];
    const double logFrequencySlackRight = x[1];
    const double curviness              = x[2];
    const double vOffset                = x[3];
    const double vScale                 = x[4];
    const double f_min                  = x[5];
    const double f_max                  = x[6];

    const double tanh_curviness = tanh(curviness);

    for (int i = 0; i < m; ++i)
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

double wa[30000];
int lwa = 30000;

int main(void)
{
    const int m = 3001; // number of functions (i.e., data points)
    const int n = 5; // number of parameters;

    FILE * f = fopen("vf.txt", "r");

    for (int i = 0; i < m; ++i)
    {
        fscanf(f, "%lf%lf", &v_arr[i], &f_arr[i]);
        //printf("%f %f\n", v_arr[i], f_arr[i]);
    }
    fclose(f);

    double f_min = +INFINITY;
    double f_max = -INFINITY;

    for (int i = 0; i < m; ++i)
    {
        f_min = fmin(f_min, f_arr[i]);
        f_max = fmax(f_max, f_arr[i]);

    }

    printf("frequency range: %f .. %f\n", f_min, f_max);

    double x[7]; // parameters: (logFrequencySlackLeft, logFrequencySlackRight, curviness, vOffset, vScale, f_min, f_max);

    double fvec[3001];

    int iwa[5];

    x[0] = 20.3;
    x[1] = 20.0;
    x[2] = 1.0;
    x[3] = 1.5;
    x[4] = 0.40;
    x[5] = f_min;
    x[6] = f_max;

    double tol = 1e-12;
    int info;

    lmdif1(fcn, m, n, x, fvec, tol, &info, iwa, wa, lwa);

    printf("info: %d\n", info);

    return 0;
}
