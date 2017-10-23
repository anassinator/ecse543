#!/usr/bin/env python

import math
from matrix import Matrix2D, SparseMatrix2D


def cholesky_solve(A, b, half_bandwidth=None):
    """Solves for x in Ax=b by Cholesky decomposition.

    Args:
        A: Real, symmetric, positive-definite matrix of order n.
        b: Column vector of size n.
        half_bandwidth: Half bandwidth for banded solution, optional.

    Returns:
        Column vector x of size n.
    """
    if not A.is_symmetric():
        raise ValueError("matrix is not symmetric")
    if not A.is_square():
        raise ValueError("matrix is not square: {}".format(A.shape))

    n = A.rows
    A = A.clone()
    L = SparseMatrix2D.zeros(n, n)
    y = b.clone()

    for j in range(n):
        if A[j, j] <= 0:
            raise ValueError("matrix is not positive-definite")

        L[j, j] = math.sqrt(A[j, j])
        y[j, 0] /= L[j, j]

        for i in range(j + 1, n):
            if half_bandwidth and i > j + half_bandwidth:
                break

            L[i, j] = A[i, j] / L[j, j]
            y[i, 0] -= L[i, j] * y[j, 0]

            for k in range(j + 1, i + 1):
                if half_bandwidth and k > j + half_bandwidth:
                    break
                A[i, k] -= L[i, j] * L[k, j]

    x = Matrix2D.zeros(n, 1)
    for i in range(n - 1, -1, -1):
        s = 0.
        for j in range(i + 1, n):
            s += L[j, i] * x[j, 0]

        x[i, 0] = (y[i, 0] - s) / L[i, i]

    return x


if __name__ == "__main__":
    def allclose(x1, x2, epsilon=1e-6):
        return all(map(lambda x: abs(x[0] - x[1]) <= epsilon, zip(x1, x2)))

    def test_solve(A, b, expected_x):
        print("A = {}".format(A.tolist()))
        print("b = {}".format(b.flatten()))
        x = cholesky_solve(A, b)
        print("x = {}".format(x.flatten()))

        if allclose(x.flatten(), expected_x.flatten()):
            print("correct")
        else:
            print("incorrect: expected {}".format(expected_x.flatten()))

    # Test with n = 3
    A = Matrix2D([[25, 15, -5], [15, 18, 0], [-5, 0, 11]])
    x = Matrix2D([[1], [2], [3]])
    b = A.dot(x)
    test_solve(A, b, x)

    # Test with n = 3
    A = Matrix2D([[25, 15, -5], [15, 18, 0], [-5, 0, 11]])
    x = Matrix2D([[0], [0], [0]])
    b = A.dot(x)
    test_solve(A, b, x)
