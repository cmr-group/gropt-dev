### Problem

$$
\begin{aligned}
& \underset{x}{\text{minimize}}
& & \dfrac{1}{2} x^T P x + q^T x \\
& \text{subject to}
& & l \leq Ax \leq u \\
% & & & f_i(x) \leq b_i, \; i = 1, \ldots, m.
\end{aligned}
$$

### Residuals

$$
r_{prim} = Ax - z
$$

$$
r_{dual} = Px + q + A^Ty
$$

### Indirect Linear Solver

$$
(P + \sigma I + \rho A^T A)\tilde{x}^{k+1} = \sigma x^{k} - q + A^T(\rho z^k - y^k)
$$

$$
\tilde{z}^{k+1} = A \tilde{x}^{k+1}
$$

### rho Balancing

$$
\rho^{k+1} = \rho^{k} \sqrt{ 
    \frac{\left\Vert r_{prim}^k \right\Vert _{\infty} / \, \textrm{max} \{ \left\Vert Ax^k \right\Vert _{\infty}, \left\Vert z^k \right\Vert _{\infty}  \} }
    {\left\Vert r_{dual}^k \right\Vert _{\infty} / \,  \textrm{max} \{ \left\Vert Px^k \right\Vert _{\infty}, \left\Vert A^Ty^k \right\Vert _{\infty}, \left\Vert q \right\Vert _{\infty}  \} }
    }
$$