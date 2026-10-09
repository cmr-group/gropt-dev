## Convenience Function

::: gropt.solve
    options:
        heading_level: 3

## SolveResult

::: gropt.SolveResult
    options:
        filters: ["!^_"]
        heading_level: 3

## SolverGroptSDMM

The default solver. The options shared with `SolverOSQP` (iteration limits, the inner linear solver, debug histories and warm starts) are listed with it. See the [solver control example](../examples/solvers/solver_control.ipynb).

::: gropt.SolverGroptSDMM
    options:
        filters: ["!^_"]
        inherited_members: true
        heading_level: 3

## SolverOSQP

::: gropt.SolverOSQP
    options:
        filters: ["!^_"]
        heading_level: 3

## Convergence Plots

With `solver.extra_debug = True`, `solver.get_debug()` returns per-iteration histories, which these functions plot.

::: gropt.debug_plotting.plot_convergence
    options:
        heading_level: 3

::: gropt.debug_plotting.scrub_hist_X
    options:
        heading_level: 3

::: gropt.debug_plotting.scrub_hist_X_widget
    options:
        heading_level: 3

::: gropt.debug_plotting.scrub_hist_X_mpl
    options:
        heading_level: 3

## Equilibration

::: gropt.equilibrate
    options:
        heading_level: 3

::: gropt.estimate_row_col_norms
    options:
        heading_level: 3

::: gropt.rescale_eq_vecs
    options:
        heading_level: 3

::: gropt.estimate_spec_norm
    options:
        heading_level: 3

::: gropt.estimate_individual_spec_norm
    options:
        heading_level: 3
