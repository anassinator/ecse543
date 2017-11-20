import math
from matrix import Matrix2D
from cholesky import cholesky_solve


def conjgrad(A, b, x=None, tol=1e-6, on_iteration=None):
    """Solves for x in Ax=b by the Conjugate Gradient method.

    Args:
        A: A matrix.
        b: b vector.
        x: Initial guess. Default: All zeros.
        tol: Tolerance. Default: 1e-10.
        on_iteration: Callback to call on each iteration.
            Signature: (iteration, x, residual) -> None.
            Default: None.

    Returns:
        x estimate.
    """
    if x is None:
        x = Matrix2D.zeros(*b.shape)

    r = b - A.dot(x);
    p = r;
    residual_old = r.T.dot(r)

    iteration = 0
    while True:
        iteration += 1
        Ap = A.dot(p);

        alpha = residual_old / (p.T.dot(Ap))
        x = x + p * alpha
        r = r - Ap * alpha

        if on_iteration is not None:
            on_iteration(iteration, x, r)

        residual_new = r.T.dot(r)
        if math.sqrt(residual_new) < tol:
            break

        p = r + p * (residual_new / residual_old)

        residual_old = residual_new

    return x


A = Matrix2D([
    # 0   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18
    [-4,  1,  2,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0],  # 0
    [ 1, -4,  0,  2,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0],  # 1
    [ 1,  0, -4,  1,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0],  # 2
    [ 0,  1,  1, -4,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0],  # 3
    [ 0,  0,  1,  0, -4,  1,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0],  # 4
    [ 0,  0,  0,  1,  1, -4,  1,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0],  # 5
    [ 0,  0,  0,  0,  0,  1, -4,  1,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0],  # 6
    [ 0,  0,  0,  0,  0,  0,  1, -4,  1,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0],  # 7
    [ 0,  0,  0,  0,  0,  0,  0,  2, -4,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0],  # 8
    [ 0,  0,  0,  0,  1,  0,  0,  0,  0, -4,  1,  0,  0,  0,  1,  0,  0,  0,  0],  # 9
    [ 0,  0,  0,  0,  0,  1,  0,  0,  0,  1, -4,  1,  0,  0,  0,  1,  0,  0,  0],  # 10
    [ 0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  1, -4,  1,  0,  0,  0,  1,  0,  0],  # 11
    [ 0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  1, -4,  1,  0,  0,  0,  1,  0],  # 12
    [ 0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  2, -4,  0,  0,  0,  0,  1],  # 13
    [ 0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0, -4,  1,  0,  0,  0],  # 14
    [ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  1, -4,  1,  0,  0],  # 15
    [ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  1, -4,  1,  0],  # 16
    [ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  1, -4,  1],  # 17
    [ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  2, -4],  # 18
])


b = Matrix2D([
    [  0],  # 0
    [-15],  # 1
    [  0],  # 2
    [-15],  # 3
    [  0],  # 4
    [  0],  # 5
    [-15],  # 6
    [-15],  # 7
    [-15],  # 8
    [  0],  # 9
    [  0],  # 10
    [  0],  # 11
    [  0],  # 12
    [  0],  # 13
    [  0],  # 14
    [  0],  # 15
    [  0],  # 16
    [  0],  # 17
    [  0],  # 18
])


if __name__ == "__main__":
    import sys
    import numpy as np
    import matplotlib.pyplot as plt

    # Store residuals.
    iterations = []
    residuals = []

    def on_iteration(iteration, x, r):
        iterations.append(iteration)
        residuals.append(r)

    print("Cholesky Decomposition solution")
    # Make it symmetric, positive definite by multiplying both sides by A.T.
    chol_sol = cholesky_solve(A.T * A, A.T * b)
    chol_sol.print()

    print("Conjugate Gradient solution")
    cg_sol = conjgrad(A, b, on_iteration=on_iteration)
    cg_sol.print()

    print("Error")
    (cg_sol - chol_sol).map(abs).print()

    # Print the value at the desired location (0.06, 0.04).
    # Here this is node 11.
    print("Value at (0.06, 0.04): {:7.4f} V".format(cg_sol[11, 0]))

    # Compute the capacitance.
    phi_sum = cg_sol[0, 0] / 2
    phi_sum += cg_sol[2, 0] + cg_sol[4, 0] + cg_sol[9, 0]
    phi_sum += cg_sol[14, 0] * 2
    phi_sum += cg_sol[15, 0] + cg_sol[16, 0] + cg_sol[17, 0]
    phi_sum += cg_sol[18, 0] / 2
    phi_sum *= 4
    Q = 8.854e-12 * phi_sum
    C = Q / 15.0
    print("Total capacitance: {:7.4e} F/m".format(C))

    # Print the residuals for plotting natively in LaTeX.
    if len(sys.argv) >= 2 and sys.argv[1] == "--plot":
        two_norm = np.linalg.norm(residuals, ord=2, axis=1)
        inf_norm = np.linalg.norm(residuals, ord=np.inf, axis=1)
        print("2-norm")
        print(list(zip(iterations, two_norm.flat)))
        print("Infinity norm")
        print(list(zip(iterations, inf_norm.flat)))
