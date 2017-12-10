#!/usr/bin/env python3

from matrix import Matrix2D


def lagrange_interpolator(X, Y):
    def multiplier(j, x):
        if isinstance(x, Matrix2D):
            y = Matrix2D.ones(x.rows, x.cols)
        else:
            y = 1.0

        for i in range(X.rows):
            if i == j:
                continue
            y *= (x - X[i, 0]) / (X[j, 0] - X[i, 0])
        return y

    def interpolate(x):
        if isinstance(x, Matrix2D):
            y = Matrix2D.zeros(x.rows, x.cols)
        else:
            y = 0.0

        for j in range(X.rows):
            y += multiplier(j, x) * Y[j, 0]
        return y

    return interpolate


def hermite_cubic_interpolator(X, Y):
    def get_coeffs(x0, x1, y0, y1, dy0_dx, dy1_dx):
        A = Matrix2D([
            [x0**3, x0**2, x0, 1.0],
            [x1**3, x1**2, x1, 1.0],
            [3 * x0**2, 2 * x0, 1.0, 0.0],
            [3 * x1**2, 2 * x1, 1.0, 0.0],
        ])
        b = Matrix2D([y0, y1, dy0_dx, dy1_dx]).T
        return A.solve(b).T

    def apply(coeffs, x):
        x = Matrix2D([x**3, x**2, x, 1.0]).T
        return coeffs.dot(x)

    coefficients = []
    dy1_dx = 0
    for i in range(X.rows - 1):
        x0 = X[i, 0]
        x1 = X[i + 1, 0]

        y0 = Y[i, 0]
        y1 = Y[i + 1, 0]

        dy0_dx = dy1_dx
        dy1_dx = (y1 - y0) / (x1 - x0)

        coefficients.append(get_coeffs(x0, x1, y0, y1, dy0_dx, dy1_dx))

    def interpolate(x):
        if isinstance(x, Matrix2D):
            return x.map(interpolate)

        for i in range(X.rows - 1):
            if x <= X[i + 1, 0]:
                break
        return apply(coefficients[i], x)

    return interpolate


def piecewise_linear_interpolator(X, Y):
    def interpolate(x):
        if isinstance(x, Matrix2D):
            return x.map(interpolate)

        for i in range(X.rows - 1):
            if x <= X[i + 1, 0]:
                break

        x0 = X[i, 0]
        x1 = X[i + 1, 0]

        y0 = Y[i, 0]
        y1 = Y[i + 1, 0]

        m = (y1 - y0) / (x1 - x0)
        y = m * (x - x0) + y0

        return y
    return interpolate


if __name__ == "__main__":
    from bh_curve import get_data
    import matplotlib.pyplot as plt

    def plot(b, h, B, H):
        plt.plot(H, B, "r.", label="data points")
        plt.plot(h, b, "k", label="interpolated curve")
        plt.xlabel("H (A/m)")
        plt.ylabel("B (T)")
        plt.legend()
        plt.show()

    x_test = Matrix2D.linspace(1.9)

    # a.
    print("a")
    indices = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]
    X, Y = get_data(indices)
    X.print(name="X")
    Y.print(name="Y")
    interp = lagrange_interpolator(X, Y)
    y_test = interp(x_test)
    plot(x_test, y_test, X, Y)

    # b.
    print("b")
    indices = [0.0, 1.3, 1.4, 1.7, 1.8, 1.9]
    X, Y = get_data(indices)
    X.print(name="X")
    Y.print(name="Y")
    interp = lagrange_interpolator(X, Y)
    y_test = interp(x_test)
    plot(x_test, y_test, X, Y)

    # c.
    print("c")
    indices = [0.0, 1.3, 1.4, 1.7, 1.8, 1.9]
    X, Y = get_data(indices)
    X.print(name="X")
    Y.print(name="Y")
    interp = hermite_cubic_interpolator(X, Y)
    y_test = interp(x_test)
    plot(x_test, y_test, X, Y)
