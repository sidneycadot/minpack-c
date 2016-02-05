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
