#!/usr/bin/env python3

import sys
from matrix import Matrix2D
from cholesky import cholesky_solve


def incidence_matrix(n, m, branches):
    elem_index = 0
    incidence = Matrix2D.zeros(n, m, dtype=int)
    for branch in branches:
        start, end = branch[:2]
        incidence[end, elem_index] = 1
        incidence[start, elem_index] = -1
        elem_index += 1

    return incidence


def Y_matrix(n, m, branches):
    Y = Matrix2D.zeros(m, m)
    for i, branch in enumerate(branches):
        r_value = branch[3]
        if r_value:
            Y[i, i] = 1.0 / r_value
    return Y


def J_matrix(n, m, branches):
    J = Matrix2D.zeros(m, 1)
    for i, branch in enumerate(branches):
        j_value = branch[2]
        J[i, 0] = j_value
    return J


def E_matrix(n, m, branches):
    E = Matrix2D.zeros(m, 1)
    for i, branch in enumerate(branches):
        e_value = branch[4]
        E[i, 0] = e_value
    return E


def solve(n, m, branches):
    A = incidence_matrix(n, m, branches)
    Y = Y_matrix(n, m, branches)
    J = J_matrix(n, m, branches)
    E = E_matrix(n, m, branches)

    # Delete last row to get reduced incidence matrix.
    A = A[:-1]

    left_hand_side = (A * Y) * A.T
    right_hand_side = A * (J - Y * E)

    v = cholesky_solve(left_hand_side, right_hand_side)
    return v


def parse_count(line):
    return int(line.strip().split(" ")[1])


def parse_branch(line):
    start, end, j_value, r_value, e_value = line.strip().split(" ")
    return int(start), int(end), float(j_value), float(r_value), float(e_value)


def parse_file(f):
    # Parse node and branch count.
    n = parse_count(f.readline())
    m = parse_count(f.readline())

    # Parse branches.
    branches = [
        parse_branch(f.readline())
        for _ in range(m)
    ]

    return n, m, branches


def print_solve(n, m, branches):
    print("node count:", n)
    print("branch count:", m)
    v = solve(n, m, branches)
    print("voltages", v.flatten())


if __name__ == "__main__":
    if len(sys.argv) == 1:
        n, m, branches = parse_file(sys.stdin)
        print_solve(n, m, branches)
    else:
        for filename in sys.argv[1:]:
            print(filename)
            with open(filename) as f:
                n, m, branches = parse_file(f)
                print_solve(n, m, branches)
