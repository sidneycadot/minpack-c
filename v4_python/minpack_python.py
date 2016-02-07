#! /usr/bin/env python3

import sys
import numpy as np

MACHINE_EPSILON = 2.22044604926e-16
PIVOT_OFFSET = 0

def square(x):
    return x * x

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


def qrsolv(r, ipvt, diag, qtb):

    r = r.copy()

    (m, n) = r.shape

    x = np.zeros(n)
    sdiag = np.zeros(n)

    for j in range(0, n):

        for i in range(j, n):
            r[i, j] = r[j, i]

        x[j] = r[j, j]

    wa = qtb.copy()

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
                        sin_ = 1.0 / np.sqrt(1.0 + cot_ ** 2)
                        cos_ = sin_ * cot_
                    else:
                        tan_ = sdiag[k] / r[k, k]
                        cos_ = 1.0 / np.sqrt(1.0 + tan_ ** 2)
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

    return (r, x, sdiag)

def lmpar(r, ipvt, diag, qtb, delta, par):

    (m, n) = r.shape

    wa1 = np.full(n, np.nan)
    wa2 = np.full(n, np.nan)

    x     = np.full(n, np.nan)
    sdiag = np.full(n, np.nan)

    p1   = 0.1
    p001 = 0.001

    # Compute and store in x the Gauss-Newton direction. If the
    # Jacobian is rank-deficient, obtain a least squares solution.

    nsing = n

    for j in range(0, n):

        wa1[j] = qtb[j]

        if r[j, j] == 0.0 and nsing == n:
            nsing = j

        if nsing < n:
            wa1[j] = 0.0

    print("lmpar: nsing =", n)

    print("wa1 before", wa1)

    for k in range(1, nsing + 1):

        j = nsing - k

        wa1[j] /= r[j, j]

        temp = wa1[j]

        for i in range(0, j):
            wa1[i] -= r[i, j] * temp

    print("wa1 after", wa1)

    for j in range(0, n):
        l = ipvt[j] - PIVOT_OFFSET
        x[l] = wa1[j]

    # Initialize the iteration counter.
    # Evaluate the function at the origin, and test
    # for acceptance of the Gauss-Newton direction.

    iter_ = 0

    for j in range(0, n):
        wa2[j] = diag[j] * x[j]

    print("diag", diag)
    print("x", x)
    print("wa2", wa2)

    dxnorm = np.linalg.norm(wa2)
    fp = dxnorm - delta

    print("lmpar: dxnorm =", dxnorm)
    print("lmpar: fp =", fp)

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

            temp = np.linalg.norm(wa1)

            parl = fp / delta / temp / temp

        # Calculate an upper bound, paru, for the zero of the function.

        for j in range(0, n):

            summ = 0.0

            for i in range(0, j + 1):
                summ += r[i, j] * qtb[i]

            l = ipvt[j] - PIVOT_OFFSET

            wa1[j] = summ / diag[l]

        gnorm = np.linalg.norm(wa1)

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

            temp = np.sqrt(par)

            for j in range(0, n):
                wa1[j] = temp * diag[j]

            (r, x, sdiag) = qrsolv(r, ipvt, wa1, qtb)

            for j in range(0, n):
                wa2[j] = diag[j] * x[j]

            dxnorm = np.linalg.norm(wa2)

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

            temp = np.linalg.norm(wa1)

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

    if iter_ == 0:
        par = 0.0

    return (par, x, sdiag)

def fdjac2(fcn, x, fvec, epsfcn):

    eps = np.sqrt(max(epsfcn, MACHINE_EPSILON))

    m = len(fvec)
    n = len(x)

    fjac = np.empty((m, n))

    for j in range(0, n):

        temp = x[j]

        h = eps * abs(temp)

        if h == 0.0:
            h = eps

        x[j] = temp + h

        fjac[:, j] = (fcn(x) - fvec) / h

        x[j] = temp

    return fjac

def lmdif(fcn, m, n, x, ftol, xtol, gtol, maxfev, epsfcn, diag, mode, factor):

    p1    = 0.1
    p5    = 0.5
    p25   = 0.25
    p75   = 0.75
    p0001 = 1e-4

    info = 0
    nfev = 0

    iflag = 0

    # Check the input parameters for errors.

    if n <= 0 or m < n or ftol < 0.0 or xtol < 0.0 or gtol < 0.0 or maxfev <= 0 or factor <= 0.0:
        return

    if mode == 2:
        # only check values of 'diag' array if it was pre-supplied.
        for j in range(0, n):
            if diag[j] <= 0.0:
                return

    # Evaluate the function at the starting point
    # and calculate its norm.

    fvec = fcn(x)
    nfev = 1

    fnorm = np.linalg.norm(fvec)

    print("fnorm right after eval: {:.20f}".format(fnorm))

    # Initialize Levenberg-Marquardt parameter and iteration counter.

    par = 0.0
    iter_ = 1

    while True: # Outer loop.

        # Beginning of the outer loop.

        # Calculate the Jacobian matrix.

        fjac = fdjac2(fcn, x, fvec, epsfcn)

        if False:

            print("*** fjac ***")
            print()
            print(fjac)
            print()

        nfev += n

        # Compute the QR factorization of the Jacobian.

        (fjac, ipvt, rdiag, acnorm) = qrfac(fjac, True)

        if False:

            print("*** after qrfac() ***")

            print("---> fjac")
            print()
            print(fjac)
            print()

            print("---> ipvt")
            print()
            print(ipvt)
            print()

            print("---> rdiag")
            print()
            print(rdiag)
            print()

            print("---> acnorm")
            print()
            print(acnorm)
            print()

        # On the first iteration and if mode is 1, scale according
        # to the norms of the columns of the initial Jacobian.

        if iter_ == 1:

            if mode != 2:

                diag = np.empty(n)

                for j in range(0, n):
                    diag[j] = acnorm[j]
                    if acnorm[j] == 0.0:
                        diag[j] = 1.0

            # On the first iteration, calculate the norm of the scaled x
            # and initialize the step bound delta.

            xnorm = np.linalg.norm(diag * x)
            delta = factor * xnorm

            if delta == 0.0:
                delta = factor

        if True:
            print("diag", diag)
            print("xnorm", xnorm)
            print("delta", delta)

        # Form (q transpose)*fvec and store the first n components in
        # qtf.

        wa4 = fvec.copy()

        qtf = np.empty(n)

        for j in range(0, n):

            if fjac[j, j] != 0.0:

                summ = 0.0

                for i in range(j, m):

                    summ += fjac[i, j] * wa4[i]

                temp = -summ / fjac[j, j]

                for i in range(j, m):
                    wa4[i] += fjac[i, j] * temp

            fjac[j, j] = rdiag[j]

            qtf[j] = wa4[j]

        if False:
            print("HELLO")
            print("fjac")
            print()
            print(fjac)
            print()
            print("qtf")
            print(qtf)
            print()

        # Compute the norm of the scaled gradient.

        gnorm = 0.0

        if fnorm != 0.0:

            for j in range(0, n):

                l = ipvt[j] - PIVOT_OFFSET

                if acnorm[l] != 0.0:

                    summ = 0.0

                    for i in range(0, j + 1):
                        summ += fjac[i, j] * (qtf[i] / fnorm)

                    gnorm = max(gnorm, abs(summ / acnorm[l]))

        print ("gnorm:", gnorm)

        # Test for convergence of the gradient norm.

        if gnorm <= gtol:
            return 4

        # Rescale if necessary.

        if mode != 2:

            for j in range(0, n):
                diag[j] = max(diag[j], acnorm[j])

        while True: # Inner loop.

            # Beginning of the inner loop.

            # Determine the Levenberg-Marquardt parameter.

            (par, wa1, wa2) = lmpar(fjac, ipvt, diag, qtf, delta, par)

            # Store the direction p and x + p. calculate the norm of p.

            wa1 = -wa1
            wa2 = x + wa1

            wa3 = diag * wa1

            pnorm = np.linalg.norm(wa3)

            print(pnorm)

            # On the first iteration, adjust the initial step bound.

            if iter_ == 1:
                delta = min(delta, pnorm)

            # Evaluate the function at x + p and calculate its norm.

            wa4 = fcn(wa2)

            nfev += 1

            fnorm1 = np.linalg.norm(wa4)

            print("fnorm1", fnorm1)

            # Compute the scaled actual reduction.

            actred = -1.0

            if p1 * fnorm1 < fnorm:
                actred = 1.0 - square(fnorm1 / fnorm)

            # Compute the scaled predicted reduction and
            # the scaled directional derivative.

            for j in range(0, n):
                wa3[j] = 0.0

                l = ipvt[j] - PIVOT_OFFSET

                temp = wa1[l]

                for i in range(0, j + 1):
                    wa3[i] += fjac[i, j] * temp

            temp1 = np.linalg.norm(wa3) / fnorm
            temp2 = np.sqrt(par) * pnorm / fnorm

            prered = square(temp1) + 2 * square(temp2)
            dirder = -(square(temp1) + square(temp2))

            print("prered, dirder", prered, dirder)

            # Compute the ratio of the actual to the predicted
            # reduction.

            ratio = 0.0

            if prered != 0.0:
                ratio = actred / prered

            # Update the step bound.

            print("RATIO:", ratio)

            if ratio <= p25:

                if actred >= 0.0:
                    temp = p5
                else:
                    temp = p5 * dirder / (dirder + p5 * actred)

                if p1 * fnorm1 >= fnorm or temp < p1:
                    temp = p1

                delta = temp * min(delta, pnorm / p1)
                par /= temp

            elif par == 0.0 or ratio >= p75:
                delta = pnorm / p5
                par = p5 * par

            # Test for successful iteration.

            if ratio >= p0001:

                # Successful iteration. Update x, fvec, and their norms.

                for j in range(0, n):
                    x[j] = wa2[j]
                    wa2[j] = diag[j] * x[j]

                for i in range(0, m):
                    fvec[i] = wa4[i]

                xnorm = np.linalg.norm(wa2)
                fnorm = fnorm1

                iter_ += 1

            # Tests for convergence.

            if abs(actred) <= ftol and prered <= ftol and p5 * ratio <= 1.0 and info == 2:
                return (3, x, nfev)

            if abs(actred) <= ftol and prered <= ftol and p5 * ratio <= 1.0:
                return (1, x, nfev)

            if delta <= xtol * xnorm:
                return (2, x, nfev)

            # Tests for termination and stringent tolerances.

            if gnorm <= MACHINE_EPSILON:
                return (8, x, nfev)

            if delta <= MACHINE_EPSILON * xnorm:
                return (7, x, nfev)

            if abs(actred) <= MACHINE_EPSILON and prered <= MACHINE_EPSILON and p5 * ratio <= 1.0:
                return (6, x, nfev)

            if nfev >= maxfev:
                return (5, x, nfev)

            # End of the inner loop. Repeat if iteration unsuccessful.

            if ratio >= p0001:
                break

        # End of the outer loop.

    return (info, x, nfev)

def lmdif1(fcn, m, n, x, tol):

    if n <= 0 or m < n or tol < 0.0:
        return (info, x)

    # Call lmdif().

    maxfev = (n + 1) * 200
    ftol   = tol
    xtol   = tol
    gtol   = 0.0
    epsfcn = 0.0
    diag   = None
    mode   = 1
    factor = 100.0

    (info, x, nfev) = lmdif(fcn, m, n, x, ftol, xtol, gtol, maxfev, epsfcn, diag, mode, factor)

    if info == 8:
        info = 4

    return (info, x, nfev)

class FitFunc:
    def __init__(self, vf):
        self.v = vf["v"]
        self.f = vf["f"]
        self.f_min = min(self.f)
        self.f_max = max(self.f)

    def __call__(self, x):

        print("fcn: x = ", x)

        logFrequencySlackLeft  = x[0]
        logFrequencySlackRight = x[1]
        curviness              = x[2]
        vOffset                = x[3]
        vScale                 = x[4]

        f1 = self.f_min - np.exp(logFrequencySlackLeft)
        f2 = self.f_max + np.exp(logFrequencySlackRight)

        x = (self.f - f1) / (f2 - f1) * 2 - 1
        xx = x + np.tanh(curviness) * (1 - x) * (1 + x) / 2

        y  = np.tan((np.pi / 2) * xx)
        vModel = vOffset + vScale * y

        residuals = vModel - self.v

        #print("fcn residuals:", residuals)
        #print()

        return residuals

def test():

    vf_dtype = [("v", np.float64), ("f", np.float64)]

    vf = np.loadtxt("vf.txt", dtype = vf_dtype)

    selection = (np.abs(vf["v"] % 0.500) < 1e-10)

    vf = vf[selection]

    fcn = FitFunc(vf)

    x = np.array([20.0, 20.0, 0.0, 1.5, 0.40])

    tol = 1e-12

    (info, x, nfev) = lmdif1(fcn, len(vf), 5, x, tol)

    print()
    print()
    print("info:", info)
    print("x", x)
    print("nfev", nfev)

def main():
    np.set_printoptions(linewidth = np.inf, precision = 20)
    test()

if __name__ == "__main__":
    main()
