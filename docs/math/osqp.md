All copied from here: [Arxiv Link](https://arxiv.org/abs/1711.08013)

## Problem


Consider the following optimization problem

\begin{equation}\label{pb:main}
\begin{array}{ll}
	\mbox{minimize}   & (1/2) x^\tpose P x + q^\tpose x \\
	\mbox{subject to} & A x \in \mathcal{C},
\end{array}
\end{equation}

where $x\in \reals^{n}$ is the decision variable.

The objective function is defined by a positive semidefinite matrix $P \in \symm_{+}^{n}$ and a vector $q \in \reals^{n}$, and the constraints by a matrix $A\in \reals^{m \times n}$ and a nonempty, closed and convex set $\mathcal{C} \subseteq \reals^m$.

## Primal dual solution

We will find it convenient to rewrite problem~\eqref{pb:main} by introducing an additional decision variable $z\in\reals^m$, to obtain the equivalent problem

\begin{equation}\label{pb:equivalent}
\begin{array}{ll}
	\mbox{minimize}   & (1/2) x^\tpose P x + q^\tpose x \\
	\mbox{subject to} & A x = z \\
	& z \in \mathcal{C}.
\end{array}
\end{equation}

We can write the optimality conditions of problem $\eqref{pb:equivalent}$ as

\begin{align}
    & Ax = z, \label{eq:prim_feas} \\
    & Px + q + A^\tpose y = 0, \label{eq:dual_feas}\\
    & z \in \mathcal{C}, \quad y\in N_{\mathcal{C}}(z),\label{eq:separating_hyperplane}
\end{align}

where $y\in \reals^{m}$ is the Lagrange multiplier associated with the constraint $Ax=z$ and $N_{\mathcal{C}}(z)$ denotes the normal cone of $\mathcal{C}$ at $z$.

If there exist $x \in \reals^n$, $z \in \reals^{m}$ and $y \in \reals^{m}$ that satisfy the conditions above, then we say that $(x,z)$ is a *primal* and $y$ is a *dual* solution to problem $\eqref{pb:equivalent}$.

We define the primal and dual residuals of problem $\eqref{pb:main}$ as

\begin{align}
    r_{\rm prim} &\eqdef Ax - z, \label{eq:res_prim}\\
    r_{\rm dual} &\eqdef Px + q + A^\tpose y. \label{eq:res_dual}
\end{align}

## Algorithm

1.  **given** initial values $x^0$, $z^0$, $y^0$ and parameters $\rho>0$, $\sigma>0$, $\alpha \in (0,2)$
2.  **Repeat**
3.  &nbsp;&nbsp;&nbsp;&nbsp; $(\tilde{x}^{k+1}, \nu^{k+1}) \gets$ solve linear system $\begin{bmatrix} P + \sigma I & A^\tpose \\ A & -\rho^{-1}I \end{bmatrix} \begin{bmatrix} \tilde{x}^{k+1} \\ \nu^{k+1} \end{bmatrix}= \begin{bmatrix}\sigma x^{k} - q \\ z^{k} - \rho^{-1} y^k \end{bmatrix}$
4.  &nbsp;&nbsp;&nbsp;&nbsp; $\tilde{z}^{k+1} \gets z^k + \rho^{-1}(\nu^{k+1} - y^{k})$
5.  &nbsp;&nbsp;&nbsp;&nbsp; $x^{k+1} \gets \alpha \tilde{x}^{k+1} + (1-\alpha)x^{k}$
6.  &nbsp;&nbsp;&nbsp;&nbsp; $z^{k+1} \gets \Pi\left(\alpha \tilde{z}^{k+1} + (1-\alpha)z^{k} + \rho^{-1} y^{k} \right)$
7.  &nbsp;&nbsp;&nbsp;&nbsp; $y^{k+1} \gets y^{k} + \rho \left( \alpha \tilde{z}^{k+1} + (1-\alpha)z^{k} - z^{k+1} \right)$
8.  **Until** termination criterion is satisfied

where $\sigma>0$ and $\rho>0$ are the *step-size parameters*, $\alpha \in (0,2)$ is the *relaxation parameter*, and $\Pi$ denotes the Euclidean projection onto $\mathcal{C}$.

### Indirect Linear Solver

With regards to step 3 of the algorithm, in our applications we solve indirectly:

With large-scale QP, factoring linear system in step 3 might be prohibitive. In these cases it might be more convenient to use an indirect method by solving instead the linear system

\begin{equation*}
	\left(P + \sigma I + \rho A^\tpose A \right)\tilde{x}^{k+1} = \sigma x^{k} - q + A^\tpose (\rho z^{k} - y^{k})
\end{equation*}

obtained by eliminating $\nu^{k+1}$ from (step 3). We then compute $\tilde{z}^{k+1}$ as $\tilde{z}^{k+1} = A\tilde{x}^{k+1}$.

Note that the coefficient matrix in the above linear system is always positive definite.

The linear system can therefore be solved with an iterative scheme such as the conjugate gradient method.

When the linear system is solved up to some predefined accuracy, we terminate the method.
