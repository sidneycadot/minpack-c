#! /usr/bin/env python3

import numpy as np

MACHINE_EPSILON = 2.22044604926e-16
PIVOT_OFFSET = 1

def qrfac(a, pivot_flag):

    a = a.copy()

    (m, n) = a.shape

    # Compute the initial column norms and initialize several arrays.

    acnorm = np.linalg.norm(a, axis = 0)

    rdiag  = acnorm.copy()
    wa     = rdiag.copy()

    if pivot_flag:
        ipvt = np.arange(n) + PIVOT_OFFSET
    else:
        ipvt = None

    # Reduce a to r with Householder transformations.

    minmn = min(m, n)

    for j in range(minmn):

        if pivot_flag:

            # Bring the column of largest norm into the pivot position.

            kmax = j

            for k in range(j, n):
                if rdiag[k] > rdiag[kmax]:
                    kmax = k

            if kmax != j:

                for i in range(0, m):
                    swap = a[i, j]
                    a[i, j] = a[i, kmax]
                    a[i, kmax] = swap

                rdiag[kmax] = rdiag[j]
                wa[kmax] = wa[j]

                swap = ipvt[j]
                ipvt[j] = ipvt[kmax]
                ipvt[kmax] = swap

        # Compute the Householder transformation to reduce the
        # j-th column of a to a multiple of the j-th unit vector.

        ajnorm = np.linalg.norm(a[j:, j])

        if ajnorm != 0.0:

            if a[j, j] < 0.0:
                ajnorm = -ajnorm

            for i in range(j, m):
                a[i, j] /= ajnorm

            a[j, j] += 1.0

            # Apply the transformation to the remaining columns
            # and update the norms.

            for k in range(j + 1, n):

                a[j:, k] -= np.dot(a[j:, j], a[j:, k]) / a[j, j] * a[j:, j]

                if pivot_flag and rdiag[k] != 0.0:

                    rdiag[k] *= np.sqrt(max(0.0, 1.0 - (a[j, k] / rdiag[k]) ** 2))

                    if 0.05 * (rdiag[k] / wa[k]) ** 2 <= MACHINE_EPSILON:
                        wa[k] = rdiag[k] = np.linalg.norm(a[j + 1:, k])

        rdiag[j] = -ajnorm

    return (a, ipvt, rdiag, acnorm)


def qrsolv(r, const int ldr, const int *ipvt, const diag, qtb, x, sdiag):

    assert r.shape[0] == r.shape[1]
    n = r.shape[0]

    for j in range(0, n):

        for i in range(j, n):
            r[i, j] = r[j, i]

        x[j] = r[j, j]
        wa[j] = qtb[j]

    # Eliminate the diagonal matrix d using a Givens rotation.

    for j in range(0, n):

        # Prepare the row of d to be eliminated, locating the
        # diagonal element using p from the QR factorization.

        l = ipvt[j] - PIVOT_OFFSET

        if diag[l] != 0.0:

            for k in range(j, n):
                sdiag[k] = 0.0

            sdiag[j] = diag[l]

            # The transformations to eliminate the row of d
            # modify only a single element of (q transpose)*b
            # beyond the first n, which is initially zero.

            qtbpj = 0.0

            for k in range(j, n):

                # Determine a givens rotation which eliminates the
                # appropriate element in the current row of d.

                if sdiag[k] != 0.0:

                    if abs(r[k, k]) < abs(sdiag[k]):
                        cot_ = r[k, k] / sdiag[k]
                        sin_ = 1.0 / sqrt(1.0 + cot_ ** 2)
                        cos_ = sin_ * cot_
                    else:
                        tan_ = sdiag[k] / r[k, k]
                        cos_ = 1.0 / sqrt(1.0 + tan_ ** 2)
                        sin_ = cos_ * tan_

                    # Compute the modified diagonal element of r and
                    # the modified element of ((q transpose)*b,0).

                    r[k, k] = cos_ * r[k, k] + sin_ * sdiag[k]

                    temp = wa[k]

                    wa[k] =  cos_ * temp + sin_ * qtbpj
                    qtbpj = -sin_ * temp + cos_ * qtbpj

                    # Accumulate the tranformation in the row of s.

                    for i in range(k + 1, n):
                        temp = r[i, k]
                        r[i, k]  =  cos_ * temp + sin_ * sdiag[i]
                        sdiag[i] = -sin_ * temp + cos_ * sdiag[i]

        # Store the diagonal element of s and restore
        # the corresponding diagonal element of r.

        sdiag[j] = r[j, j]
        r[j, j] = x[j]

    # Solve the triangular system for z. If the system is
    # singular, then obtain a least squares solution.

    nsing = n

    for j in range(0, n):

        if sdiag[j] == 0.0 and nsing == n:
            nsing = j

        if nsing < n:
            wa[j] = 0.0

    for k in range(1, nsing + 1):

        j = nsing - k

        summ = 0.0

        for i in range(j + 1, nsing):
            summ += r[i, j] * wa[i]

        wa[j] = (wa[j] - summ) / sdiag[j]

    # Permute the components of z back to components of x.

    for j in range(0, n):
        l = ipvt[j] - PIVOT_OFFSET
        x[l] = wa[j]


def lmpar(n, r, ipvt, const double *diag, double *qtb, delta, double par, double x, double sdiag):

    assert r.shape[0] == r.shape[1]
    n = r.shape[0]

    const double p1   = 0.1;
    const double p001 = 0.001;

    # Compute and store in x the Gauss-Newton direction. If the
    # Jacobian is rank-deficient, obtain a least squares solution.

    nsing = n;

    for j in range(0, n):

        wa1[j] = qtb[j]

        if r[j, j] == 0.0 and nsing == n:
            nsing = j

        if nsing < n:
            wa1[j] = 0

    for k in range(1, nsing):

        j = nsing - k

        wa1[j] /= r[j, j]

        temp = wa1[j]

        for i in range(0, j):
            wa1[i] -= r[i, j] * temp

        for j in range(0, n):
            l = ipvt[j] - PIVOT_OFFSET
            x[l] = wa1[j]

    # Initialize the iteration counter.
    # Evaluate the function at the origin, and test
    # for acceptance of the Gauss-Newton direction.

    iter_ = 0

    for j in range(0, n):
        wa2[j] = diag[j] * x[j]

    dxnorm = enorm(n, wa2)
    fp = dxnorm - delta

    if fp > p1 * delta:

        # If the Jacobian is not rank deficient, the Newton
        # step provides a lower bound, parl, for the zero of
        # the function. Otherwise set this bound to zero.

        parl = 0.0

        if nsing >= n:

            for j in range(0, n):
                l = ipvt[j] - PIVOT_OFFSET
                wa1[j] = diag[l] * (wa2[l] / dxnorm)

            for j in range(0, n):

                summ = 0.0
                for i in range(0, j):
                    summ += r[i, j] * wa1[i]

                wa1[j] = (wa1[j] - summ) / r[j, j]

            temp = enorm(n, wa1)

            parl = fp / delta / temp / temp

        # Calculate an upper bound, paru, for the zero of the function.

        for j in range(0, n):

            summ = 0.0

            for i in range(0, j + 1):
                summ += r[i, j] * qtb[i]

            l = ipvt[j] - PIVOT_OFFSET

            wa1[j] = summ / diag[l]

        gnorm = enorm(n, wa1)

        paru = gnorm / delta

        if paru == 0.0:
            paru = MACHINE_MINPOS / min(delta, p1)

        # If the input par lies outside of the interval (parl, paru),
        # set par to the closer endpoint.

        par = max(par, parl)
        par = min(par, paru)

        if par == 0.0:
            par = gnorm / dxnorm

        while True:

            # Beginning of an iteration.

            iter_ += 1

            # Evaluate the function at the current value of par.

            if par == 0.0:
                par = max(MACHINE_MINPOS, p001 * paru)

            temp = sqrt(par)

            for j in range(0, n):
                wa1[j] = temp * diag[j]

            qrsolv(n, r, ldr, ipvt, wa1, qtb, x, sdiag, wa2)

            for j in range(0, n):
                wa2[j] = diag[j] * x[j]

            dxnorm = enorm(n, wa2)
            temp = fp
            fp = dxnorm - delta

            # If the function is small enough, accept the current value
            # of par. Also test for the exceptional cases where parl
            # is zero or the number of iterations has reached 10.

            if abs(fp) <= p1 * delta or (parl == 0.0 and fp <= temp and temp < 0.0) or iter == 10:
                break

            # Compute the Newton correction.

            for j in range(0, n):
                l = ipvt[j] - PIVOT_OFFSET
                wa1[j] = diag[l] * (wa2[l] / dxnorm)

            for j in range(0, n):

                wa1[j] /= sdiag[j]

                temp = wa1[j]

                for i in range(j + 1, n):
                    wa1[i] -= r[i, j] * temp

            temp = enorm(n, wa1)

            parc = fp / delta / temp / temp

            # Depending on the sign of the function, update parl or paru.

            if fp > 0.0:
                parl = max(parl, par)
            else:
                paru = min(paru, par)

            # Compute an improved estimate for par.

            par = max(parl, par + parc)

            # End of an iteration.

    # Termination.

    if iter == 0:
        par = 0.0

    return par

def lmdif(lmdif_fcn_ptr fcn, const int m, const int n, double *x, double *fvec, const double ftol, const double xtol, const double gtol, const int maxfev,
           const double epsfcn, double *diag, const int mode, const double factor, const int nprint, int *info, int *nfev,
           double *fjac, const int ldfjac, int *ipvt, double *qtf, double *wa1, double *wa2, double *wa3, double *wa4:

    const double p1    = 0.1;
    const double p5    = 0.5;
    const double p25   = 0.25;
    const double p75   = 0.75;
    const double p0001 = 1e-4;

    # Local variables.

    double delta = 0.0; # Initialize to suppress warning.
    double xnorm = 0.0; # Initialize to suppress warning.

    # Parameter adjustments.

    *info = 0;
    *nfev = 0;

    iflag = 0

    # Check the input parameters for errors.

    if n <= 0 or m < n or ldfjac < m or ftol < 0.0 or xtol < 0.0 or gtol < 0.0 or maxfev <= 0 or factor <= 0.0:
        goto L_TERMINATE;

    if mode == 2:

        for (int j = 0; j < n; ++j)
        {
            if (diag[j] <= 0.0)
            {
                goto L_TERMINATE;
            }
        }

    # Evaluate the function at the starting point
    # and calculate its norm.

    iflag = 1;
    (*fcn)(m, n, x, fvec, &iflag);
    *nfev = 1;

    if iflag < 0:
        goto L_TERMINATE

    fnorm = enorm(m, fvec)

    # Initialize Levenberg-Marquardt parameter and iteration counter.

    par = 0.0
    iter_ = 1

    while True: # Outer loop.

        # Beginning of the outer loop.

        # Calculate the Jacobian matrix.

        iflag = 2;
        fdjac2(fcn, m, n, x, fvec, fjac, ldfjac, &iflag, epsfcn, wa4);
        *nfev += n;

        if iflag < 0:
            goto L_TERMINATE

        # If requested, call fcn to enable printing of iterates.

        if nprint > 0:
            iflag = 0
            if (iter_ - 1) % nprint == 0:
                (*fcn)(m, n, x, fvec, &iflag)

            if iflag < 0:
                goto L_TERMINATE

        # Compute the QR factorization of the Jacobian.

        qrfac(m, n, fjac, ldfjac, true, ipvt, n, wa1, wa2, wa3)

        # On the first iteration and if mode is 1, scale according
        # to the norms of the columns of the initial Jacobian.

        if iter_ == 1:

            if mode != 2:

                for (int j = 0; j < n; ++j)

                    diag[j] = wa2[j]

                    if wa2[j] == 0.0:
                        diag[j] = 1.0

            # On the first iteration, calculate the norm of the scaled x
            # and initialize the step bound delta.

            for (int j = 0; j < n; ++j)
            {
                wa3[j] = diag[j] * x[j];
            }

            xnorm = enorm(n, wa3)
            delta = factor * xnorm

            if delta == 0.0:
                delta = factor

        # Form (q transpose)*fvec and store the first n components in
        # qtf.

        for i in range(0, m):
            wa4[i] = fvec[i]

        for (int j = 0; j < n; ++j)
        {
            if fjac[j, j] != 0.0:
            {
                summ = 0.0

                for (int i = j; i < m; ++i)
                {
                    sum += fjac[i + j * ldfjac] * wa4[i];
                }

                temp = -summ / fjac[j, j];

                for (int i = j; i < m; ++i)
                {
                    wa4[i] += fjac[i, j] * temp
                }
            }

            fjac[j + j * ldfjac] = wa1[j];
            qtf[j] = wa4[j];
        }

        # Compute the norm of the scaled gradient.

        gnorm = 0.0

        if fnorm != 0.0:
            for j in range(0, n):
            {
                l = ipvt[j] - PIVOT_OFFSET

                if wa2[l] != 0.0:
                {
                    summ = 0.0

                    for (int i = 0; i < j + 1; ++i)
                    {
                        summ += fjac[i + j * ldfjac] * (qtf[i] / fnorm);
                    }
                    gnorm = max(gnorm, abs(summ / wa2[l]));
                }
            }

        # Test for convergence of the gradient norm.

        if gnorm <= gtol:
            info = 4

        if info != 0:
            goto L_TERMINATE

        # Rescale if necessary.

        if mode != 2:

            for j in range(0, n):
            {
                diag[j] = max(diag[j], wa2[j])
            }

        while True: # Inner loop.

            # Beginning of the inner loop.

            # Determine the Levenberg-Marquardt parameter.

            lmpar(n, fjac, ldfjac, ipvt, diag, qtf, delta, &par, wa1, wa2, wa3, wa4)

            # Store the direction p and x + p. calculate the norm of p.

            for j in range(0, n):
            {
                wa1[j] = -wa1[j];
                wa2[j] = x[j] + wa1[j];
                wa3[j] = diag[j] * wa1[j];
            }

            pnorm = enorm(n, wa3)

            # On the first iteration, adjust the initial step bound.

            if iter_ == 1:
                delta = min(delta, pnorm)

            # Evaluate the function at x + p and calculate its norm.

            iflag = 1;
            (*fcn)(m, n, wa2, wa4, &iflag);
            ++(*nfev);

            if iflag < 0:
            {
                goto L_TERMINATE;
            }

            fnorm1 = enorm(m, wa4)

            # Compute the scaled actual reduction.

            actred = -1.0

            if p1 * fnorm1 < fnorm:
                actred = 1.0 - square(fnorm1 / fnorm)

            # Compute the scaled predicted reduction and
            # the scaled directional derivative.

            for (int j = 0; j < n; ++j)
            {
                wa3[j] = 0.0;

                const int l = ipvt[j] - PIVOT_OFFSET;

                const double temp = wa1[l];

                for (int i = 0; i < j + 1; ++i)
                {
                    wa3[i] += fjac[i + j * ldfjac] * temp;
                }
            }

            temp1 = enorm(n, wa3) / fnorm
            temp2 = sqrt(par) * pnorm / fnorm

            prered = square(temp1) + 2 * square(temp2)
            dirder = -(square(temp1) + square(temp2))

            # Compute the ratio of the actual to the predicted
            # reduction.

            double ratio = 0.0

            if prered != 0.0:
                ratio = actred / prered

            # Update the step bound.

            if (ratio <= p25)
            {
                double temp;

                if (actred >= 0.0)
                {
                    temp = p5;
                }
                else:
                {
                    temp = p5 * dirder / (dirder + p5 * actred);
                }

                if (p1 * fnorm1 >= fnorm or temp < p1)
                {
                    temp = p1;
                }

                delta = temp * fmin(delta, pnorm / p1);
                par /= temp;
            }
            elif par == 0.0 or ratio >= p75:
            {
                delta = pnorm / p5
                par = p5 * par
            }

            # Test for successful iteration.

            if ratio >= p0001:

                # Successful iteration. Update x, fvec, and their norms.

                for (int j = 0; j < n; ++j)
                {
                    x[j] = wa2[j];
                    wa2[j] = diag[j] * x[j];
                }

                for (int i = 0; i < m; ++i)
                {
                    fvec[i] = wa4[i];
                }

                xnorm = enorm(n, wa2)
                fnorm = fnorm1

                iter_ += 1

            # Tests for convergence.

            if (abs(actred) <= ftol && prered <= ftol && p5 * ratio <= 1.0)
            {
                info = 1
            }

            if delta <= xtol * xnorm:
                info = 2

            if (abs(actred) <= ftol && prered <= ftol && p5 * ratio <= 1.0 && *info == 2)
            {
                info = 3
            }

            if (*info != 0)
            {
                goto L_TERMINATE;
            }

            # Tests for termination and stringent tolerances.

            if (*nfev >= maxfev)
            {
                *info = 5;
            }

            if abs(actred) <= MACHINE_EPSILON && prered <= MACHINE_EPSILON && p5 * ratio <= 1.0:
                info = 6

            if delta <= MACHINE_EPSILON * xnorm:
                info = 7

            if (gnorm <= MACHINE_EPSILON)
            {
                *info = 8;
            }

            if info != 0:
                goto L_TERMINATE;

            # End of the inner loop. Repeat if iteration unsuccessful.

            if (ratio >= p0001)
            {
                break;
            }


        # End of the outer loop.

L_TERMINATE:

    # Termination, either normal or user imposed.

    if iflag < 0:
    {
        info = iflag;
    }

    iflag = 0
    if nprint > 0:
    {
        f(m, n, x, fvec, iflag)
    }

    return info


np.set_printoptions(precision = 8, linewidth = np.inf, formatter = { 'float': lambda x : "{:12.6f}".format(x)} )

m = 3
n = 7

axi = np.expand_dims(1 + np.arange(m), 1).repeat(n, 1)
axj = np.expand_dims(1 + np.arange(n), 0).repeat(m, 0)

a = 1 / np.sin(10 * axi + axj)

print("a before qrfac:\n")
print(a)
print("\n")

(a, ipvt, rdiag, acnorm) = qrfac(a, True)

print("a after qrfac:\n")
print(a)
print("\n")

print("ipvt after qrfac:\n")
print(ipvt)
print("\n")

print("rdiag after qrfac:\n")
print(rdiag)
print("\n")

print("acnorm after qrfac:\n")
print(acnorm)
print("\n")
