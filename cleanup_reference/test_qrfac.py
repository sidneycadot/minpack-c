#! /usr/bin/env python3

import numpy as np

m = 3
n = 7

axi = np.expand_dims(1 + np.arange(m), 1).repeat(n, 1)
axj = np.expand_dims(1 + np.arange(n), 0).repeat(m, 0)

a = 1 / np.sin(10 * axi + axj)

print(a)

