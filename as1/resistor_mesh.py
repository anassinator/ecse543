#!/usr/bin/env python3

from linear_resistive_network import LinearResistiveNetwork


def resistor_mesh(rows, cols, resistance, test_voltage=10.):
    n = rows * cols
    m = 2 * rows * cols - rows - cols

    resistance = float(resistance)

    branches = []
    for start in range(n):
        for end in (start + 1, start + cols):
            if (start + 1) % cols == 0 and end == start + 1:
                continue

            if end >= n:
                continue

            branch = start, end, 0., resistance, 0.
            branches.append(branch)

    if test_voltage:
        m += 1
        branches.append((n - 1, 0, 0., resistance, test_voltage))

    return LinearResistiveNetwork(n, m, branches)


def solve_resistance(N, network, resistance, test_voltage=10., banded=True):
    half_bandwidth = N + 1 if banded else None
    v = network.solve_for_voltages(half_bandwidth)
    v_measured = v[-1]
    z = v_measured * resistance / (test_voltage - v_measured)
    return z


if __name__ == "__main__":
    import time

    def print_solution_test(N, resistance, banded):
        print("N", N)
        print("Banded", banded)
        network = resistor_mesh(2 * N, N, resistance)
        start = time.time()
        total_resistance = solve_resistance(N, network, resistance,
                                            banded=banded)
        end = time.time()
        print("R_eq", total_resistance)
        print("time", end - start)

    resistance = 1000.0
    for N in range(2, 11):
        print_solution_test(N, resistance, True)
        print_solution_test(N, resistance, False)
        print("")
