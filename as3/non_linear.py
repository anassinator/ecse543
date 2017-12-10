#!/usr/bin/env python3


def newton_raphson(f, f_prime, x, tol=1e-6, max_iter=50, on_iteration=None):
    for i in range(max_iter):
        x -= f(x) / f_prime(x)

        if on_iteration is not None:
            on_iteration(i + 1, x, f(x))

        if abs(f(x) / f(0)) <= tol:
            break

    return x


def vector_newton_raphson(f, jac, x, tol=1e-6, max_iter=50, on_iteration=None):
    x0 = Matrix2D.zeros(x.rows, x.cols)
    for i in range(max_iter):
        x -= jac(x).inverse().dot(f(x))

        if on_iteration is not None:
            on_iteration(i + 1, x, f(x))

        if sum(f(x).map(abs).flatten()) <= tol:
            break

    return x


def successive_substitution(f,
                            x,
                            tol=1e-6,
                            max_iter=50,
                            c=6.5e-9,
                            on_iteration=None):
    for i in range(max_iter):
        x -= c * f(x)

        if on_iteration is not None:
            on_iteration(i + 1, x, f(x))

        if abs(f(x) / f(0)) <= tol:
            break

    return x


if __name__ == '__main__':
    import math
    from matrix import Matrix2D
    from bh_curve import get_data
    import matplotlib.pyplot as plt
    from interpolation import piecewise_linear_interpolator

    S = 0.0001
    X, Y = get_data()
    interp = piecewise_linear_interpolator(X, Y)

    def H(psi):
        B = psi / S
        return interp(B)

    def Hprime(psi):
        B = psi / S
        for i in range(X.rows - 1):
            if B <= X[i + 1, 0]:
                break
        x0 = X[i, 0]
        x1 = X[i + 1, 0]

        y0 = Y[i, 0]
        y1 = Y[i + 1, 0]

        m = (y1 - y0) / (x1 - x0)
        return m / S

    def f(psi):
        return 3.9789e7 * psi + 0.30 * H(psi) - 8000

    def fprime(psi):
        return 3.9789e7 + 0.30 * Hprime(psi)

    def on_iteration(i, x, f_x):
        print("iteration: {}, x = {:E}, f(x) = {:+E}".format(i, x, f_x))

    print("Newton-Raphson method")
    x = newton_raphson(f, fprime, 0, on_iteration=on_iteration)
    print("Magnetic flux = {:E}".format(x), "Wb")
    print()

    print("Successive substitution method")
    x = successive_substitution(f, 0, on_iteration=on_iteration)
    print("Magnetic flux = {:E}".format(x), "Wb")
    print()

    vT = 0.025
    I_sA = 0.6e-6
    I_sB = 1.2e-6
    R = 500
    E = 0.220

    def f(v):
        v1 = v[0, 0]
        v2 = v[1, 0]
        return Matrix2D([
            [1 / R * v1 - E / R + I_sA * (math.exp((v1 - v2) / vT) - 1)],
            [
                I_sB * (math.exp(v2 / vT) - 1) - I_sA * (math.exp(
                    (v1 - v2) / vT) - 1)
            ],
        ])

    def fprime(v):
        v1 = v[0, 0]
        v2 = v[1, 0]
        return Matrix2D([
            [
                1 / R + I_sA / vT * math.exp((v1 - v2) / vT),
                -I_sA / vT * math.exp((v1 - v2) / vT)
            ],
            [
                -I_sA / vT * math.exp((v1 - v2) / vT), (I_sA * math.exp(
                    (v1 - v2) / vT) + I_sB * math.exp(v2 / vT)) / vT
            ],
        ])

    errors = []

    def on_iteration(i, v, f_v):
        print("iteration:", i)
        v. print(name="V")
        error = sum(f_v.map(abs).flatten())
        print("error: {:+E}".format(error))
        errors.append(error)

    def plot(errors):
        iterations = list(range(1, len(errors) + 1))
        plt.plot(iterations, errors, "k")
        plt.xlabel("Iteration")
        plt.ylabel("Error")
        plt.show()

    print("Question 3")
    v = Matrix2D.zeros(2, 1)
    v = vector_newton_raphson(f, fprime, v, on_iteration=on_iteration)
    v. print(name="V")
    plot(errors)
