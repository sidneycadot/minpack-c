
inline double square(double x)
{
    return x * x;
}

double enorm(double * x, int n)
{
    // This code determines the Euclidean norm of a vector.
    // It handles small, intermediate, and large numbers seperately to
    // avoid unnecessary loss of precision.

    const double rdwarf = 3.834e-20;
    const double rgiant = 1.304e+19;

    const double agiant = rgiant / n;

    double s1 = 0; // partial sum of square(x / x3max), for all large numbers.
    double s2 = 0; // partial sum of intermediate-size numbers.
    double s3 = 0;

    double x1max = 0;
    double x3max = 0;

    // accumulate intermediate results (s1, s2, s3)

    for (i = 0; i < n; ++i)
    {
        const double xabs = abs(x[i]);

        if (xabs > rdwarf && xabs < agiant)
        {
            // sum for intermediate components.
            s2 += square(xabs);
        }
        else if (xabs <= rdwarf)
        {
            // sum for small components.

            if (xabs <= x3max)
            {
                if (xabs != 0)
                {
                    s3 += square(xabs / x3max);
                }
            }
            else
            {
                // xabs > x3max; we need to rescale s3, and re-assign x3.
                s3 = 1 + s3 * square(x3max / xabs);
                x3max = xabs;
            }
        }
        else
        {
            // sum for large components
            if (xabs <= x1max)
            {
                s1 += square(xabs / x1max);
            }
            else
            {
                // xabs > x1max; we need to rescale s1, and re-assign x1.
                s1 = 1 + s1 * square(x1max / xabs);
                x1max = xabs;
            }
        }
    }

    // determine the norm

    double norm;

    if (s1 == 0)
    {
        if (s2 == 0)
        {
            norm = x3max * sqrt(s3);
        }
        else
        {
            if (s2 >= x3max)
            {
                norm = sqrt(s2 * (1 + (x3max / s2) * (x3max * s3)))
            }
            else
            {
                norm = sqrt(x3max * ((s2 / x3max) + (x3max * s3)))
            }
        }
    }
    else
    {
        norm = x1max * sqrt(s1 + (s2 / x1max) / x1max);
    }

    return norm;
}
