#!/usr/bin/env python3

import math
from matrix import Matrix2D
from rectangle import Rectangle


class ElectrostaticProblem(object):

    def __init__(self, outer_conductor, inner_conductor, outer_v, inner_v):
        self.outer_conductor = outer_conductor
        self.inner_conductor = inner_conductor

        self.outer_voltage = float(outer_v)
        self.inner_voltage = float(inner_v)

    def uniformly_discretize(self, h):
        discretized = UniformlyDiscretizedProblem(
            self.outer_conductor.width,
            self.outer_conductor.height, h)

        discretized[:] = self.outer_voltage

        bottom_left = self.inner_conductor.bottom_left
        left, bottom = discretized.to_grid_coordinate(*bottom_left)
        top_right = self.inner_conductor.top_right
        right, top = discretized.to_grid_coordinate(*top_right)

        discretized[bottom:top + 1, left:right + 1] = self.inner_voltage

        return discretized


class UniformlyDiscretizedProblem(object):

    def __init__(self, width, height, h):
        self.h = h
        self.rows = int(width / h) + 1
        self.cols = int(height / h) + 1
        self.grid = Matrix2D.zeros(self.rows, self.cols)

    def __getitem__(self, idx):
        return self.grid.__getitem__(idx)

    def __setitem__(self, idx, value):
        return self.grid.__setitem__(idx, value)

    def on_edge(self, x, y):
        return x in (0, self.rows - 1) or y in (0, self.cols - 1)

    def to_grid_coordinate(self, x, y):
        i = max(0, min(int(y / self.h), self.cols - 1))
        j = max(0, min(int(x / self.h), self.rows - 1))
        return i, j

    def to_original_coordinate(self, i, j):
        return float(j * self.h), float(i * self.h)


def sor(problem, w, h, tolerance=1e-5):
    grid = problem.uniformly_discretize(h)

    inner_y, inner_x = grid.to_grid_coordinate(
        *problem.inner_conductor.bottom_left)

    residual = float("inf")
    iterations = 0
    while residual > tolerance:
        iterations += 1
        for x in range(1, grid.rows):
            for y in range(1, grid.cols):
                if x >= inner_x and y >= inner_y:
                    continue

                if x == grid.rows - 1 and y == grid.cols - 1:
                    # Deal with top right corner.
                    # Shouldn't happen since the conductor is there.
                    s = 2 * grid[x - 1, y] + 2 * grid[x, y - 1]
                elif x == grid.rows - 1:
                    # Deal with right edge.
                    s = 2 * grid[x - 1, y] + grid[x, y + 1] + grid[x, y - 1]
                elif y == grid.cols - 1:
                    # Deal with top edge.
                    s = grid[x + 1, y] + grid[x - 1, y] + 2 * grid[x, y - 1]
                else:
                    s = grid[x + 1, y] + grid[x - 1, y] + \
                        grid[x, y + 1] + grid[x, y - 1]

                grid[x, y] = (1 - w) * grid[x, y] + (w / 4) * s

        residual = 0.0
        for x in range(1, grid.rows - 1):
            for y in range(1, grid.cols - 1):
                local_residual = -4 * grid[x, y]
                local_residual += grid[x + 1, y]
                local_residual += grid[x - 1, y]
                local_residual += grid[x, y + 1]
                local_residual += grid[x, y - 1]

                residual = max(residual, local_residual)

    return grid, iterations


if __name__ == "__main__":
    # Conductor sizes.
    OUTER_CONDUCTOR_SIZE = 0.2
    INNER_CONDUCTOR_WIDTH = 0.08
    INNER_CONDUCTOR_HEIGHT = 0.04

    # Conductor voltages.
    OUTER_CONDUCTOR_VOLTAGE = 0.0
    INNER_CONDUCTOR_VOLTAGE = 15.0

    # Since this is symmetrical we can just consider the bottom left quadrant
    # of the problem.
    outer_conductor = Rectangle(
        OUTER_CONDUCTOR_SIZE / 2,
        OUTER_CONDUCTOR_SIZE / 2,
        OUTER_CONDUCTOR_SIZE / 4,
        OUTER_CONDUCTOR_SIZE / 4)
    inner_conductor = Rectangle(
        INNER_CONDUCTOR_WIDTH / 2,
        INNER_CONDUCTOR_HEIGHT / 2,
        (OUTER_CONDUCTOR_SIZE - INNER_CONDUCTOR_WIDTH) / 2,
        (OUTER_CONDUCTOR_SIZE - INNER_CONDUCTOR_HEIGHT) / 2)
    problem = ElectrostaticProblem(
        outer_conductor,
        inner_conductor,
        OUTER_CONDUCTOR_VOLTAGE,
        INNER_CONDUCTOR_VOLTAGE)

    TARGET_COORDINATE = 0.06, 0.04

    def print_solution(discretized, iterations):
        coord = discretized.to_grid_coordinate(*TARGET_COORDINATE)
        v = discretized[coord]
        print("voltage at", TARGET_COORDINATE, "=", v, "V")
        print("iterations", iterations)

    print("SOR")
    print("varying w...")
    iterations_per_w = {}
    for m in range(10, 20):
        w = m / 10
        h = 0.02
        print("w =", w)
        print("h =", h)
        discretized, iterations = sor(problem, w, h)
        print_solution(discretized, iterations)
        iterations_per_w[w] = iterations
        discretized.grid.print()
        print("")

    # Minimum iterations from previous loop.
    min_w = min(iterations_per_w, key=lambda x: iterations_per_w[x])
    print("minimum w =", min_w, "@", iterations_per_w[min_w], "iterations")
    print("")

    print("SOR")
    print("varying h...")
    for m in range(5):
        w = min_w
        h = 0.02 / 2**m
        print("w =", w)
        print("h =", h)
        discretized, iterations = sor(problem, w, h)
        print_solution(discretized, iterations)
    print("")
