import numpy as np


def GaussElimination(A, b):

    A = np.array(A, dtype=float)
    b = np.array(b, dtype=float)

    Ab = np.hstack([A, b.reshape(-1, 1)])

    n = len(b)

    for i in range(n):
        Ab[i] = Ab[i] / Ab[i, i]

        for j in range(i + 1, n):
            factor = Ab[j, i]
            Ab[j] = Ab[j] - factor * Ab[i]

    x = np.zeros(n)
    for i in range(n - 1, -1, -1):
        x[i] = Ab[i, -1]
        for j in range(i + 1, n):
            x[i] -= Ab[i, j] * x[j]
        x[i] = x[i] / Ab[i, i]

    return x
