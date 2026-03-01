<style>
table th, table td {
  font-size: 20px !important;
}

.md-typeset table:not([class]) {
  display: table;
  width: 100%;
}

</style>


| OSQP | PARSDMM | Current GroptSDMM |
| :---: | :---: | :---: |
| $\tilde{z}$ | $s$ | $s$ |
|  $\alpha\tilde{z} + (1-\alpha)z$ | $\bar{x} = \gamma s + (1-\gamma)y$ |  $\gamma s + (1-\gamma)z$ |
| $z$ | $y$ | $z$ |
| $y$ | $-v$ | $y$ |
| $\alpha$ | $\gamma$ | $\gamma,\gamma_x$ |

------

OSQP: Single $\alpha$ value used for overrelaxation of x and primal updates  
PARSDMM: Constraint specific $\gamma$ for primal update, and nothing on x   

------

Both $\rho$ are the same thing, and both are constriant specific (even though OSQP doesnt show it in algorithm)