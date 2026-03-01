All copied from here: [Arxiv Link](https://arxiv.org/abs/1902.09699)

## Algorithm

**Algorithm: PARSDMM**  
*(Projection Adaptive Relaxed Simultaneous Direction Method of Multipliers to compute the projection onto an intersection, including automatic selection of the penalty parameters and relaxation.)*  

---

**Inputs:**  
*   $m$ // point to project  
*   $A_1, A_2, \dots, A_{p}, A_{p+1}=I_N$   &nbsp;&nbsp;// linear operators  
*   $\operatorname{prox}_{f_i,\rho_i}(w) = \mathcal{P}_{\mathcal{C}_i}(w)$ for $i=1, 2, \dots, p$   &nbsp;&nbsp;// norm/bound/cardinality/... projectors  
*   $\operatorname{prox}_{f_{p+1},\rho_{p+1}}(w) = (m+\rho_{p+1} w)/(1+\rho_{p+1})$   &nbsp;&nbsp;// prox for the squared distance from $m$  
*   select $\rho_i^0$, $\gamma_i^0$, $\text{update-frequency}$ for $\gamma$ and $\rho$  
*   optional: initial guess for $x$, $y_i$ and $v_i$  

**Initialize:**  
*   $B_i = A_i^\top A_i$ for all $i$   &nbsp;&nbsp;// pre-compute  
*   $C = \sum_{i=1}^{p+1} (\rho_i B_i )$   &nbsp;&nbsp;// pre-compute  
*   $k=1$  
  
**Loop:**

**WHILE** not converged:  
&nbsp;&nbsp;&nbsp;&nbsp;$x^{k+1} = C^{-1} \sum_{i=1}^{p+1} \Big[A_i^\top(\rho_i^k y_i^{k} + v_i^k) \Big]$ // CG, stop when condition holds   
  
&nbsp;&nbsp;&nbsp;&nbsp;**FOR** $i=1, 2, \dots, p+1$   &nbsp;&nbsp;// compute in parallel  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;$s_i^{k+1} = A_i x^{k+1}$  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;$\bar{x}_i^{k+1} = \gamma_i^k s_i^{k+1} + ( 1-\gamma_i^k ) y_i^{k}$  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;$y_i^{k+1} = \operatorname{prox}_{f_i,\rho_i}(\bar{x}_i^{k+1} - \frac{v_i^k}{\rho_i^k} )$  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;$v_i^{k+1} = v_i^k + \rho_i^k (y_i^{k+1} - \bar{x}_i^{k+1})$  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;stop if stopping conditions hold  
  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;**IF** $\text{mod}(k,\text{update-frequency}) = 1$  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;$\{\rho_i^{k+1},\gamma_i^{k+1}\}=\text{adapt-rho-gamma}(v_i^k,v_i^{k+1},y_i^{k+1},s_i^{k+1},\rho_i^k)$  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;**END IF**  
&nbsp;&nbsp;&nbsp;&nbsp;**END FOR**  

&nbsp;&nbsp;&nbsp;&nbsp;**FOR** $i=1, 2, \dots, p+1$   &nbsp;&nbsp;// update C if necessary  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;**IF** $\rho_i^{k+1} \neq \rho_i^{k}$  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;$C \leftarrow C + B_i (\rho_i^{k+1} - \rho_i^{k})$  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;**END IF**  
&nbsp;&nbsp;&nbsp;&nbsp;**END FOR**  
  
&nbsp;&nbsp;&nbsp;&nbsp;$k \leftarrow k+1$  
**END WHILE**  
  
**Output:** $x$  