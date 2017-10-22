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


def solve_resistance(network, resistance, test_voltage=10.):
    v = network.solve_for_voltages()
    v_measured = v[-1]
    z = v_measured * resistance / (test_voltage - v_measured)
    return z


if __name__ == "__main__":
    resistance = 1000.0
    for n in range(2, 11):
        print("N", n)
        network = resistor_mesh(n, 2 * n, resistance)
        total_resistance = solve_resistance(network, resistance)
        print("R_eq", total_resistance)
