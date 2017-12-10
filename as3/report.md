---
header-includes:
  - \usepackage{listings}
cref: true
---

\lstset{basicstyle=\small}
\crefformat{equation}{(#2#1#3)}
\crefformat{figure}{Figure~#2#1#3}
\crefformat{table}{Table~#2#1#3}

# ECSE 543 Numerical Methods - Assignment 3

Anass Al-Wohoush, 260575013

## Question 1

### a
![Interpolated B-H curve for M19 steel using the full-domain Lagrange polynomials on the first 6 points.](assets/interp_a.png){#fig:interp-a width=75%}

@fig:interp-a seems like a plausible B-H curve over the given range. 

### b
![Interpolated B-H curve for M19 steel using the full-domain Lagrange polynomials at $B = 0, 1.3, 1.4, 1.7, 1.8, 1.9$.](assets/interp_b.png){#fig:interp-b width=75%}

@fig:interp-b does not seem like a plausible B-H curve over the given range. 

### c
![Interpolated B-H curve for M19 steel using the cubic Hermite polynomials at $B = 0, 1.3, 1.4, 1.7, 1.8, 1.9$.](assets/interp_c.png){#fig:interp-c width=75%}

@fig:interp-c is the generated interpolation plot by setting the slopes to:

$$
\frac{y_1 - y_0}{x_1 - x_0}.
$$ 

\newpage
## Appendix

### `interpolation.py`

\lstinputlisting[language=Python]{interpolation.py}

\newpage
### `matrix.py`

\lstinputlisting[language=Python]{matrix.py}
