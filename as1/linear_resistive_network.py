#!/usr/bin/env python3

import sys
from cholesky import cholesky_solve
from matrix import Matrix2D, SparseMatrix2D


class LinearResistiveNetwork(object):

    def __init__(self, n, m, branches, sparse=True):
        self.n = n
        self.m = m
        self.branches = branches

        self._matrix_class = SparseMatrix2D if sparse else Matrix2D

        self.A = self._A()
        self.Y = self._Y()
        self.J = self._J()
        self.E = self._E()

    def _A(self):
        elem_index = 0
        A = self._matrix_class.zeros(self.n, self.m, dtype=int)
        for branch in self.branches:
            start, end = branch[:2]
            A[end, elem_index] = 1
            A[start, elem_index] = -1
            elem_index += 1

        # Delete the first row to get reduced incidence matrix.
        A = A[1:]

        return A

    def _Y(self):
        Y = SparseMatrix2D.zeros(self.m, self.m)
        for i, branch in enumerate(self.branches):
            r_value = branch[3]
            if r_value:
                Y[i, i] = 1.0 / r_value
        return Y

    def _J(self):
        J = self._matrix_class.zeros(self.m, 1)
        for i, branch in enumerate(self.branches):
            j_value = branch[2]
            J[i, 0] = j_value
        return J

    def _E(self):
        E = self._matrix_class.zeros(self.m, 1)
        for i, branch in enumerate(self.branches):
            e_value = branch[4]
            E[i, 0] = e_value
        return E

    def solve_for_voltages(self, half_bandwidth=None):
        left_hand_side = (self.A * self.Y) * self.A.T
        right_hand_side = self.A * (self.J - self.Y * self.E)

        v = cholesky_solve(left_hand_side, right_hand_side, half_bandwidth)
        return v.flatten()

    @staticmethod
    def from_file(f):

        def parse_count(line):
            return int(line.strip().split(" ")[1])

        def parse_branch(line):
            values = line.strip().split(" ")
            start, end = map(int, values[:2])
            j_value, r_value, e_value = map(float, values[2:])
            return start, end, j_value, r_value, e_value

        # Parse node and branch count.
        n = parse_count(f.readline())
        m = parse_count(f.readline())

        # Parse branches.
        branches = [parse_branch(f.readline()) for _ in range(m)]

        return LinearResistiveNetwork(n, m, branches)


if __name__ == "__main__":

    def print_solve(network):
        print("node count:", network.n)
        print("branch count:", network.m)
        print("voltages", network.solve_for_voltages())

    if len(sys.argv) == 1:
        network = LinearResistiveNetwork.from_file(sys.stdin)
        print_solve(network)
    else:
        for filename in sys.argv[1:]:
            print(filename)
            with open(filename) as f:
                network = LinearResistiveNetwork.from_file(f)
                print_solve(network)
            print("")
