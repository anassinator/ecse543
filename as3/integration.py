#!/usr/bin/env python3

from matrix import linspace


def gauss_legendre(f, a, b, N, xs=None):
    if xs is None:
        xs = linspace(a, b, N)

    y = 0
    for prev_x, curr_x in zip(xs, xs[1:]):
        dx = curr_x - prev_x
        mid = (prev_x + curr_x) / 2
        y += f(mid) * dx

    return y


if __name__ == '__main__':
    import math
    import matplotlib.pyplot as plt

    def plot_log_errors(errors):
        ns, errs = zip(*errors)
        plt.plot(ns, errs, "ks")
        plt.plot(ns, errs, "k")
        plt.xlabel("Log segment count")
        plt.ylabel("Log error")
        plt.show()

    # sin(x)
    truth = 1 - math.cos(1)
    errors = []
    for N in range(1, 21):
        I = gauss_legendre(math.sin, 0, 1, N)
        error = abs(I - truth)
        errors.append((math.log10(N), math.log10(error)))
    plot_log_errors(errors)

    # ln(x)
    truth = -1
    errors = []
    for N in range(10, 210, 10):
        I = gauss_legendre(math.log, 0, 1, N)
        error = abs(I - truth)
        errors.append((math.log10(N), math.log10(error)))
    plot_log_errors(errors)

    # uneven segments ln(x)
    N = 10
    xs = [(x / N)**math.e for x in range(N + 1)]
    I_even = gauss_legendre(math.log, 0, 1, N)
    I_uneven = gauss_legendre(math.log, 0, 1, N, xs)
    print("evenly distributed error: {:E}".format(abs(I_even - truth)))
    print("unevenly distributed error: {:E}".format(abs(I_uneven - truth)))
