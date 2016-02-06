
//////////////////////////
// TestVF_Performance.c //
//////////////////////////

#include <math.h>
#include <assert.h>
#include <stdio.h>
#include <time.h>

#include "VF_Table.h"

#include "minpack_c.h"

struct fcn_context_t
{
    double f_min;
    double f_max;
    const struct vf_t * vf;
    int stride;
    int nfev;
};

void fcn(const int m, const int n, double *x, double *fvec, const enum lm_iflag_t iflag, void *fcn_context)
{
    (void)n; // 5 parameters
    (void)iflag;

    struct fcn_context_t * context = fcn_context;

    ++context->nfev;

    const double PI = 3.141592653589793238462643383;

    const double logFrequencySlackLeft  = x[0];
    const double logFrequencySlackRight = x[1];
    const double curviness              = x[2];
    const double vOffset                = x[3];
    const double vScale                 = x[4];

    const double tanh_curviness = tanh(curviness);

    for (int i = 0; i < m; ++i)
    {
        const double v = context->vf[i * context->stride].v;
        const double f = context->vf[i * context->stride].f;

        const double f1 = context->f_min - exp(logFrequencySlackLeft);
        const double f2 = context->f_max + exp(logFrequencySlackRight);

        const double x  = (f - f1) / (f2 - f1) * 2 - 1;
        const double xx = x + tanh_curviness * (1 - x) * (1 + x) / 2;
        const double y  = tan((PI / 2) * xx);

        const double vModel  = vOffset + vScale * y;

        fvec[i] = (vModel -v);
    }
}

double wa[30000];
const int lwa = 30000;

int main(void)
{
    const int n = 5; // number of parameters

    for (int m = 5; m <= VF_TABLE_SIZE; ++m)
    {
        if ((VF_TABLE_SIZE - 1) % (m - 1))
        {
            // No snug fit.
            continue;
        }

        // Snug fit! Proceed.

        int stride = (VF_TABLE_SIZE - 1) / (m - 1);

        double f_min = +INFINITY;
        double f_max = -INFINITY;

        for (int i = 0; i < m; ++i)
        {
            f_min = fmin(f_min, vf_table[i * stride].f);
            f_max = fmax(f_max, vf_table[i * stride].f);
        }

        printf("========== testcase: m = %d, stride = %d, frequency range = { %f .. %f GHz }\n\n", m, stride, f_min / 1e9, f_max / 1e9);

        struct fcn_context_t context;

        context.f_min  = f_min;
        context.f_max  = f_max;
        context.vf     = vf_table;
        context.stride = stride;

        const double tol = 1e-12;

        double x[5]; // parameters: (logFrequencySlackLeft, logFrequencySlackRight, curviness, vOffset, vScale)

        double fvec[m];

        int iwa[5];

        const int num_repeats = 5;

        double best_duration = INFINITY;

        for (int rep = 0; rep < num_repeats; ++rep)
        {
            x[0] = 20.3;
            x[1] = 20.0;
            x[2] = 1.0;
            x[3] = 1.5;
            x[4] = 0.40;

            const clock_t t1 = clock();

            context.nfev = 0;

            enum lm_info_t info = lmdif1(fcn, &context, m, n, x, fvec, tol, iwa, wa, lwa);

            const clock_t t2 = clock();

            const double duration = (double)(t2 - t1) / CLOCKS_PER_SEC;

            printf("    duration = %f seconds (nfev = %d)\n", duration, context.nfev);

            best_duration = fmin(best_duration, duration);

            if (rep == 0)
            {
                printf("    return value of solver: %d\n", info);

                printf("    solution = {");
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
        }
        printf("    best duration = %f seconds\n\n", best_duration);

    } // m loop

    return 0;
}
