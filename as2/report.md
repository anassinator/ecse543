---
header-includes:
  - \usepackage{listings}
  - \usepackage{tikz}
  - \usepackage{subcaption}
  - \usepackage{unicode-math}
  - \usepackage{pgfplots}
  - \usepackage[utf8]{inputenc}
---

\lstset{basicstyle=\small}

# ECSE 543 Numerical Methods - Assignment 2

Anass Al-Wohoush, 260575013

## Question 1
Let the area of a triangle be $A = \frac{0.02^2}{2} = 0.0002$.

Given

$$
\begin{bmatrix}
  U_1 \\
  U_2 \\
  U_3 \\
\end{bmatrix} = \begin{bmatrix}
  1 + x_1 + y_1 \\
  1 + x_2 + y_2 \\
  1 + x_3 + y_3 \\
\end{bmatrix} \begin{bmatrix}
  a \\
  b \\
  c \\
\end{bmatrix}
$$

The potential at each node of a triangle with nodes $1$, $2$, and $3$ can be
expressed as:

$$
U = \sum_{i = 1}^{3} U_i \alpha_i(x, y)
$$

where $\alpha_i$ is defined as:

$$
\begin{aligned}
\alpha_1(x, y) &= \frac{1}{2A}[(x_2 y_3 - x_3 y_2) + (y_2 - y_3)x + (x_3 - x_2) y] \\
\alpha_2(x, y) &= \frac{1}{2A}[(x_3 y_1 - x_1 y_3) + (y_3 - y_1)x + (x_1 - x_3) y] \\
\alpha_3(x, y) &= \frac{1}{2A}[(x_1 y_2 - x_2 y_1) + (y_1 - y_2)x + (x_2 - x_1) y] \\
\end{aligned}
$$ {#eq:alphas}

The total energy in the triangle is given by Equation @eq:energy where the
elements of the $\symbf{S}$ matrix is defined in Equation @eq:S_elements.

$$
W = \frac{1}{2} \symbf{U}^T \symbf{S} \symbf{U}
$$ {#eq:energy}

$$
S_{ij} = A \nabla \alpha_i \cdot \nabla \alpha_j
$$ {#eq:S_elements}

### Triangle $1$-$2$-$3$
The first triangle with nodes $1$, $2$, and $3$ have their coordinates as in
Table \ref{table:triangle_1}.

Node  $x$  $y$
---- ---- ----
1    0.00 0.02
2    0.00 0.00
3    0.02 0.00

Table: Triangle $1$-$2$-$3$'s nodes' coordinates.
\label{table:triangle_1}

Plugging this into Equation @eq:alphas yields the following:

$$
\begin{aligned}
\alpha_1(x, y) &= 50 y \\
\alpha_2(x, y) &= 1 - 50 x - 50 y] \\
\alpha_3(x, y) &= 50 x \\
\end{aligned}
$$ {#eq:alphas_triangle_1}

$$
\begin{aligned}
\nabla \alpha_1 &= \begin{bmatrix} 0 & 50 \end{bmatrix} \\
\nabla \alpha_2 &= \begin{bmatrix} -50 & -50 \end{bmatrix} \\
\nabla \alpha_3 &= \begin{bmatrix} 50 & 0 \end{bmatrix} \\
\end{aligned}
$$ {#eq:gradients_triangle_1}

This yields the following $\symbf{S}$ matrix:

$$
\symbf{S_{1,2,3}} = \begin{bmatrix}
  0.5 & -0.5 & 0 \\
  -0.5 & 1 & -0.5 \\
  0 & -0.5 & 0.5 \\
\end{bmatrix}
$$

### Triangle $4$-$5$-$6$
The first triangle with nodes $4$, $5$, and $6$ have their coordinates as in
Table \ref{table:triangle_2}.

Node  $x$  $y$
---- ---- ----
4    0.02 0.02
5    0.00 0.02
6    0.02 0.00

Table: Triangle $4$-$5$-$6$'s nodes' coordinates.
\label{table:triangle_1}

Plugging this into Equation @eq:alphas yields the following:

$$
\begin{aligned}
\alpha_4(x, y) &= -1 + 50 x + 50 y \\
\alpha_5(x, y) &= 1 - 50 x \\
\alpha_6(x, y) &= 1 - 50 y \\
\end{aligned}
$$ {#eq:alphas_triangle_1}

$$
\begin{aligned}
\nabla \alpha_4 &= \begin{bmatrix} 50 & 50 \end{bmatrix} \\
\nabla \alpha_5 &= \begin{bmatrix} -50 & 0 \end{bmatrix} \\
\nabla \alpha_6 &= \begin{bmatrix} 0 & -50 \end{bmatrix} \\
\end{aligned}
$$ {#eq:gradients_triangle_1}

This yields the following $\symbf{S}$ matrix:

$$
\symbf{S_{4,5,6}} = \begin{bmatrix}
  1 & -0.5 & -0.5 \\
  -0.5 & 0.5 & 0 \\
  -0.5 & 0 & 0.5 \\
\end{bmatrix}
$$

### Conjoining
Let $\symbf{C}$ be the matrix that maps the nodes from the conjoint triangles
to the disjoint triangles:

$$
\symbf{C} = \begin{bmatrix}
  1 & 0 & 0 & 0 \\
  0 & 1 & 0 & 0 \\
  0 & 0 & 1 & 0 \\
  0 & 0 & 0 & 1 \\
  1 & 0 & 0 & 0 \\
  0 & 0 & 1 & 0 \\
\end{bmatrix}
$$

such that:

$$
\symbf{U_{disjoint}} = \symbf{C} \symbf{U_{conjoint}}
$$

To find the energy of the system, we add the energy of each triangle as
follows:

$$
\begin{aligned}
W &= W_{1,2,3} + W_{4,5,6} \\
  &= \frac{1}{2} \symbf{U_{disjoint}}^T S_{disjoint} \symbf{U_{disjoint}} \\
  &= \frac{1}{2} (\symbf{C} \symbf{U_{conjoint}})^T \symbf{S_{disjoint}} (\symbf{C} \symbf{U_{conjoint}}) \\
  &= \frac{1}{2} \symbf{U_{conjoint}}^T (\symbf{C}^T \symbf{S_{disjoint}} \symbf{C}) \symbf{U_{conjoint}} \\
\end{aligned}
$$

where

$$
\symbf{S_{disjoint}} = \begin{bmatrix}
  \symbf{S_{1,2,3}} & \symbf{0} \\
  \symbf{0} & \symbf{S_{4,5,6}} \\
\end{bmatrix} = \begin{bmatrix}
  0.5 & -0.5 & 0 & 0 & 0 & 0 \\
  -0.5 & 1 & -0.5 & 0 & 0 & 0 \\
  0 & -0.5 & 0.5 & 0 & 0 & 0 \\
  0 & 0 & 0 & 1 & -0.5 & -0.5 \\
  0 & 0 & 0 & -0.5 & 0.5 & 0 \\
  0 & 0 & 0 & -0.5 & 0 & 0.5 \\
\end{bmatrix}
$$

This yields:

$$
\begin{aligned}
\symbf{S_{conjoint}} &= \symbf{C}^T \symbf{S_{disjoint}} \symbf{C} \\
 &= \begin{bmatrix}
    1 & -0.5 & 0 & -0.5 \\
    -0.5 & 1 & -0.5 & 0 \\
    0 & -0.5 & 1 & -0.5 \\
    -0.5 & 0 & -0.5 & 1 \\
 \end{bmatrix}
\end{aligned}
$$

## Question 2

### a

The problem was set up according to Figure \ref{fig:simple2d_mesh}. This yields
34 nodes, 46 mesh triangles, and 15 constraints.

`mesh.dat` in Appendix A shows how the input file for `SIMPLE2D` was setup.

\begin{figure}
\centering
\begin{subfigure}[t]{.5\textwidth}
  \centering
  \input{assets/simple2d_mesh_triangles.tikz}
  \caption{Mesh triangles configuration}
  \label{fig:simple2d_mesh_triangles}
\end{subfigure}%
\begin{subfigure}[t]{.5\textwidth}
  \centering
  \input{assets/simple2d_mesh_nodes.tikz}
  \caption{Mesh nodes configuration}
  \label{fig:simple2d_mesh_nodes}
\end{subfigure}
\caption{Mesh configuration of the bottom left quarter of the coaxial cable.
         The grey is at $15 V$ and the thick lines are grounded at $0 V$.}
\label{fig:simple2d_mesh}
\end{figure}

### b
Table \ref{table:simple2d_output} shows the formatted output of `SIMPLE2D`. The
raw output is available in Appendix A. This shows us that the potential at
$(x, y) = (0.06, 0.04)$ (i.e. Node 20) is $5.5263 V$. As expected, this matches
the estimate from Assignment 1.

Node    $x$    $y$ Potential (V)
---- ------ ------ -------------
1    0.0000 0.1000        0.0000
2    0.0200 0.1000        4.2525
3    0.0400 0.1000        9.0919
4    0.0600 0.1000       15.0000
5    0.0000 0.0800        0.0000
6    0.0200 0.0800        3.9590
7    0.0400 0.0800        8.5575
8    0.0600 0.0800       15.0000
9    0.0800 0.0800       15.0000
10   0.1000 0.0800       15.0000
11   0.0000 0.0600        0.0000
12   0.0200 0.0600        3.0262
13   0.0400 0.0600        6.1791
14   0.0600 0.0600        9.2492
15   0.0800 0.0600       10.2912
16   0.1000 0.0600       10.5490
17   0.0000 0.0400        0.0000
18   0.0200 0.0400        1.9667
19   0.0400 0.0400        3.8834
20   0.0600 0.0400        5.5263
21   0.0800 0.0400        6.3668
22   0.1000 0.0400        6.6135
23   0.0000 0.0200        0.0000
24   0.0200 0.0200        0.9571
25   0.0400 0.0200        1.8616
26   0.0600 0.0200        2.6060
27   0.0800 0.0200        3.0360
28   0.1000 0.0200        3.1714
29   0.0000 0.0000        0.0000
30   0.0200 0.0000        0.0000
31   0.0400 0.0000        0.0000
32   0.0600 0.0000        0.0000
33   0.0800 0.0000        0.0000
34   0.1000 0.0000        0.0000

Table: SIMPLE2D output.
\label{table:simple2d_output}

### c
The capacitance per unit length can be calculated from the total energy of the
system as in Equation @eq:capacitance.

$$
C = \frac{2 W \epsilon_0}{V_s^2}
$$ {#eq:capacitance}

where $V_s = 15 V$, $\epsilon_0 = 8.854 \times 10^{-12} \frac{F}{m}$ and 

$$
W = \frac{1}{2} \symbf{U}^T \symbf{S} \symbf{U}
$$

The $S$ matrix was extracted from a modified version of the `SIMPLE2D_M` code
as `Sglob`. The modified code can be found in Appendix A. Using this, $C$ was
computed in Matlab using the `compute_capacitance.m` script also available in
Appendix A. The resulting capacitance was $52.136 \frac{pF}{m}$. Note that
this is actually 4 times the value given by Equation @eq:capacitance, since we
only evaluated the bottom left quarter of the system.

## Question 3

All the code used to solve this question are available in Appendix B.

The problem was set up according to Figure \ref{fig:cg_mesh}. This yields a
total of 19 unconstrained nodes.

\begin{figure}
\centering
\begin{subfigure}[t]{.5\textwidth}
  \centering
  \input{assets/cg_mesh_axes.tikz}
  \caption{Mesh configuration}
  \label{fig:cg_mesh_axes}
\end{subfigure}%
\begin{subfigure}[t]{.5\textwidth}
  \centering
  \input{assets/cg_mesh_nodes.tikz}
  \caption{Mesh nodes configuration}
  \label{fig:cg_mesh_nodes}
\end{subfigure}
\caption{Mesh configuration of the bottom left quarter of the coaxial cable.
         The grey is at $15 V$ and the thick lines are grounded at $0 V$.}
\label{fig:cg_mesh}
\end{figure}

Through five point second order approximation, we yield $A$ and $b$ as in
Equations @eq:A_cg and @eq:b_cg. For the nodes along the edges of symmetry
(i.e. nodes $0$, $1$, $8$, $13$, $18$), some weights needed to be doubled.

$$
\setcounter{MaxMatrixCols}{20}
A = \mbox{\scriptsize%
$\begin{bmatrix}
  -4 &  1 &  2 &  0 &  0 &  0 &  0 &  0 &  0 &  0 &  0 &  0 &  0 &  0 &  0 &  0 &  0 &  0 &  0 \\
   1 & -4 &  0 &  2 &  0 &  0 &  0 &  0 &  0 &  0 &  0 &  0 &  0 &  0 &  0 &  0 &  0 &  0 &  0 \\
   1 &  0 & -4 &  1 &  1 &  0 &  0 &  0 &  0 &  0 &  0 &  0 &  0 &  0 &  0 &  0 &  0 &  0 &  0 \\
   0 &  1 &  1 & -4 &  0 &  1 &  0 &  0 &  0 &  0 &  0 &  0 &  0 &  0 &  0 &  0 &  0 &  0 &  0 \\
   0 &  0 &  1 &  0 & -4 &  1 &  0 &  0 &  0 &  1 &  0 &  0 &  0 &  0 &  0 &  0 &  0 &  0 &  0 \\
   0 &  0 &  0 &  1 &  1 & -4 &  1 &  0 &  0 &  0 &  1 &  0 &  0 &  0 &  0 &  0 &  0 &  0 &  0 \\
   0 &  0 &  0 &  0 &  0 &  1 & -4 &  1 &  0 &  0 &  0 &  1 &  0 &  0 &  0 &  0 &  0 &  0 &  0 \\
   0 &  0 &  0 &  0 &  0 &  0 &  1 & -4 &  1 &  0 &  0 &  0 &  1 &  0 &  0 &  0 &  0 &  0 &  0 \\
   0 &  0 &  0 &  0 &  0 &  0 &  0 &  2 & -4 &  0 &  0 &  0 &  0 &  1 &  0 &  0 &  0 &  0 &  0 \\
   0 &  0 &  0 &  0 &  1 &  0 &  0 &  0 &  0 & -4 &  1 &  0 &  0 &  0 &  1 &  0 &  0 &  0 &  0 \\
   0 &  0 &  0 &  0 &  0 &  1 &  0 &  0 &  0 &  1 & -4 &  1 &  0 &  0 &  0 &  1 &  0 &  0 &  0 \\
   0 &  0 &  0 &  0 &  0 &  0 &  1 &  0 &  0 &  0 &  1 & -4 &  1 &  0 &  0 &  0 &  1 &  0 &  0 \\
   0 &  0 &  0 &  0 &  0 &  0 &  0 &  1 &  0 &  0 &  0 &  1 & -4 &  1 &  0 &  0 &  0 &  1 &  0 \\
   0 &  0 &  0 &  0 &  0 &  0 &  0 &  0 &  1 &  0 &  0 &  0 &  2 & -4 &  0 &  0 &  0 &  0 &  1 \\
   0 &  0 &  0 &  0 &  0 &  0 &  0 &  0 &  0 &  1 &  0 &  0 &  0 &  0 & -4 &  1 &  0 &  0 &  0 \\
   0 &  0 &  0 &  0 &  0 &  0 &  0 &  0 &  0 &  0 &  1 &  0 &  0 &  0 &  1 & -4 &  1 &  0 &  0 \\
   0 &  0 &  0 &  0 &  0 &  0 &  0 &  0 &  0 &  0 &  0 &  1 &  0 &  0 &  0 &  1 & -4 &  1 &  0 \\
   0 &  0 &  0 &  0 &  0 &  0 &  0 &  0 &  0 &  0 &  0 &  0 &  1 &  0 &  0 &  0 &  1 & -4 &  1 \\
   0 &  0 &  0 &  0 &  0 &  0 &  0 &  0 &  0 &  0 &  0 &  0 &  0 &  1 &  0 &  0 &  0 &  2 & -4
\end{bmatrix}$}
$$ {#eq:A_cg}

$$
b = \begin{bmatrix}
    0 \\
  -15 \\
    0 \\
  -15 \\
    0 \\
    0 \\
  -15 \\
  -15 \\
  -15 \\
    0 \\
    0 \\
    0 \\
    0 \\
    0 \\
    0 \\
    0 \\
    0 \\
    0 \\
    0 \\
\end{bmatrix}
$$ {#eq:b_cg}

### a
Unfortunately, $A$ is neither symmetric, nor positive definite, so Cholesky
decomposition fails. Instead, we multiply both sides of $A x = b$ by the
transpose of $A$ as follows:

$$
\begin{aligned}
A x &= b \\
A^T A x &= A^T b \\
A_{new} x &= b_{new}
\end{aligned}
$$

and solve that system of equations instead. This works since $A^T A$ is always
positive definite.

### b
Solving the system above yields:

\centering
\begin{lstlisting}[language=bash,caption={Output of conjugate gradient},label={lst:cg_out}]
$ python3 conjugate_gradient.py
Cholesky Decomposition solution
|  4.2525 |
|  9.0919 |
|  3.9590 |
|  8.5575 |
|  3.0262 |
|  6.1791 |
|  9.2492 |
| 10.2912 |
| 10.5490 |
|  1.9667 |
|  3.8834 |
|  5.5263 |
|  6.3668 |
|  6.6135 |
|  0.9571 |
|  1.8616 |
|  2.6060 |
|  3.0360 |
|  3.1714 |
Conjugate Gradient solution
|  4.2525 |
|  9.0919 |
|  3.9590 |
|  8.5575 |
|  3.0262 |
|  6.1791 |
|  9.2492 |
| 10.2912 |
| 10.5490 |
|  1.9667 |
|  3.8834 |
|  5.5263 |
|  6.3668 |
|  6.6135 |
|  0.9571 |
|  1.8616 |
|  2.6060 |
|  3.0360 |
|  3.1714 |
Error
|  0.0000 |
|  0.0000 |
|  0.0000 |
|  0.0000 |
|  0.0000 |
|  0.0000 |
|  0.0000 |
|  0.0000 |
|  0.0000 |
|  0.0000 |
|  0.0000 |
|  0.0000 |
|  0.0000 |
|  0.0000 |
|  0.0000 |
|  0.0000 |
|  0.0000 |
|  0.0000 |
|  0.0000 |
Value at (0.06, 0.04):  5.5263 V
Total capacitance: 5.2136e-11 F/m
\end{lstlisting}

We see that both the Cholesky Decomposition and the Conjugate Gradient yield
the same result as expected. This also matches the results from `SIMPLE2D`.

### c

Figure \ref{fig:residual_plot} shows the $\infty$-norm and $2$-norm of the
residuals over time.

\begin{figure}
\centering
\caption{Evolution of residual norms by iteration}
\label{fig:residual_plot}
\begin{tikzpicture}
\begin{axis}[
    xlabel={Iteration},
    ylabel={Residual norm ($V$)},
    xmin=-1, xmax=33,
    ymin=-1, ymax=22,
    xtick={0,5,10,15,20,25,30},
    ytick={0,5,10,15,20},
    legend pos=north east,
    ymajorgrids=true,
    grid style=dashed,
    width=\textwidth
]
 
\addplot[
    color=blue,
    ]
    coordinates {
      (1, 19.960899278339141)
      (2, 14.753041714212179)
      (3, 8.9700878134619781)
      (4, 5.5345179686305626)
      (5, 2.4462407446730094)
      (6, 1.2433040280635403)
      (7, 0.60485398671941293)
      (8, 0.32108547774443086)
      (9, 0.19311996925441099)
      (10, 0.11912483627167843)
      (11, 0.047119403517898041)
      (12, 0.028500222562492093)
      (13, 0.020154149077830013)
      (14, 0.011182631562872653)
      (15, 0.0054315279811212747)
      (16, 0.0031093254034940164)
      (17, 0.0017689420129602073)
      (18, 0.0010700883320705872)
      (19, 0.00068384652102493175)
      (20, 0.00034149762955497522)
      (21, 0.00017279881017706451)
      (22, 9.066022464251618e-05)
      (23, 5.1641856546095237e-05)
      (24, 3.6039909679127454e-05)
      (25, 2.4587275802965769e-05)
      (26, 1.2896868646942411e-05)
      (27, 6.7314170779994477e-06)
      (28, 3.5126463230051873e-06)
      (29, 1.6110581270743485e-06)
      (30, 1.1389608637599635e-06)
      (31, 6.1173581060519434e-07)
    };
    \addlegendentry{Residual 2-norm}

\addplot[
    color=red,
    ]
    coordinates {
      (1, 12.5)
      (2, 6.2704918032786887)
      (3, 4.2097860955183304)
      (4, 2.5994846067364965)
      (5, 1.1674227607383778)
      (6, 0.70615806148334681)
      (7, 0.30947026785948578)
      (8, 0.20669205866585416)
      (9, 0.089494400648726172)
      (10, 0.050885043196393073)
      (11, 0.028682070610006782)
      (12, 0.013473948611356954)
      (13, 0.0095679974973694633)
      (14, 0.0058616647413585726)
      (15, 0.0036588788154160554)
      (16, 0.0018288900919320177)
      (17, 0.00083191485692417615)
      (18, 0.00046527331722984394)
      (19, 0.00044985649678882725)
      (20, 0.00019403565090673108)
      (21, 7.7131172875220612e-05)
      (22, 4.3057894892775447e-05)
      (23, 3.2130258470471804e-05)
      (24, 2.5072514707788462e-05)
      (25, 1.2816641148231336e-05)
      (26, 5.501848769877326e-06)
      (27, 2.9132076280168852e-06)
      (28, 1.9532439469126015e-06)
      (29, 1.023919469242361e-06)
      (30, 5.7597000692899292e-07)
      (31, 3.0568353949874851e-07)
    };
    \addlegendentry{Residual $\infty$-norm}

\end{axis}
\end{tikzpicture}
\end{figure}

### d
As seen in **b**, this system yields a potential of $5.5263 V$ at
$(x,y) = (0.06,0.04)$. This matches exactly with the estimate in both Question
1 and Assignment 1 as expected.

### e
The capacitance per unit length can be computed from the potential using
Equation @eq:cg_capacitance.

$$
C = Q / V
$$ {#eq:cg_capacitance}

where $Q$ can be calculated using Gauss’s Law:

$$
Q = \epsilon_0 \oint_L E dL
$$

Setting the contour $L$ to the line formed by the grounded conductor and
approximating the integral with a discrete sum, we get:

$$
Q = \epsilon_0 \sum_n E_n l
$$

where $l$ in this case is equal to $h = 0.02$. Given that $E = \phi / l$, this
can be simplified to:

$$
Q = \epsilon_0 \sum_n \phi_n
$$

On the contour $L$ and amongst the quarter we're looking at, the only nodes we
care about are: $n = {0, 2, 4, 9, 14, 15, 16, 17, 18}$. Considering the
symmetry, we can evaluate the entire contour by simply considering nodes $0$
and $18$ (which are at the edge) half as much. We must also account for node
$14$ twice. This yields:

$$
\begin{aligned}
Q &= 4 \epsilon_0 \sum_n \phi_n \\
  &= 4 \epsilon_0 (\frac{\phi_0}{2} + \phi_2 + \phi_4 + \phi_9 + 2 \phi_{14} + \phi_{15} + \phi_{16} + \phi_{17} \frac{\phi_{18}}{2}) \\
  &= 4 \epsilon_0 (\frac{4.2525}{2} + 3.9590 + 3.0262 + 1.9667 + 2 \times 0.9571 + 1.8616 + 2.6060 + 3.0360 + \frac{3.1714}{2}) \\
  &= 4 \epsilon_0 \times 22.08165 \\
  &= 7.8204 \times 10^-10 C
\end{aligned}
$$

Plugging back into Equation @eq:cg_capacitance, yields:

$$
C = Q / V = \frac{7.8204 \times 10^-10}{15} = 52.136 \frac{pF}{m}
$$

which matches the value computed in Question 2. The script
`conjugate_gradient.py` was altered to compute this capacitance as well, and
confirms this result. The output can be seen in Listing \ref{lst:cg_out}.

\newpage
## Appendix A

### `mesh.dat`

\lstinputlisting{q2/mesh.dat}

\newpage
### `compute_capacitance.m`

\lstinputlisting[language=Matlab]{q2/compute_capacitance.m}

\newpage
### `SIMPLE2D_M.m`

\lstinputlisting[language=Matlab]{q2/SIMPLE2D_M.m}

\newpage
### `SIMPLE2D` Raw Output

```matlab
>> SIMPLE2D_M('mesh.dat')

ans =

    1.0000         0    0.1000         0
    2.0000    0.0200    0.1000    4.2525
    3.0000    0.0400    0.1000    9.0919
    4.0000    0.0600    0.1000   15.0000
    5.0000         0    0.0800         0
    6.0000    0.0200    0.0800    3.9590
    7.0000    0.0400    0.0800    8.5575
    8.0000    0.0600    0.0800   15.0000
    9.0000    0.0800    0.0800   15.0000
   10.0000    0.1000    0.0800   15.0000
   11.0000         0    0.0600         0
   12.0000    0.0200    0.0600    3.0262
   13.0000    0.0400    0.0600    6.1791
   14.0000    0.0600    0.0600    9.2492
   15.0000    0.0800    0.0600   10.2912
   16.0000    0.1000    0.0600   10.5490
   17.0000         0    0.0400         0
   18.0000    0.0200    0.0400    1.9667
   19.0000    0.0400    0.0400    3.8834
   20.0000    0.0600    0.0400    5.5263
   21.0000    0.0800    0.0400    6.3668
   22.0000    0.1000    0.0400    6.6135
   23.0000         0    0.0200         0
   24.0000    0.0200    0.0200    0.9571
   25.0000    0.0400    0.0200    1.8616
   26.0000    0.0600    0.0200    2.6060
   27.0000    0.0800    0.0200    3.0360
   28.0000    0.1000    0.0200    3.1714
   29.0000         0         0         0
   30.0000    0.0200         0         0
   31.0000    0.0400         0         0
   32.0000    0.0600         0         0
   33.0000    0.0800         0         0
   34.0000    0.1000         0         0

```

\newpage
## Appendix B

### `conjugate_gradient.py`

\lstinputlisting[language=Python]{q3/conjugate_gradient.py}

\newpage
### `cholesky.py`

\lstinputlisting[language=Python]{q3/cholesky.py}

\newpage
### `matrix.py`

\lstinputlisting[language=Python]{q3/matrix.py}
