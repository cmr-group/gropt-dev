## Solve and Result

`gropt.solve()` solves with a new `SolverGroptSDMM` and a few of its settings; it returns a `SolveResult`.

::: gropt.solve

::: gropt.SolveResult

## SolverGroptSDMM

The default solver. Every option is a property, listed in the table below and described after it; the options
shared with `SolverOSQP` (iteration limits, the inner linear solver, debug histories and warm starts) are included.
See the [solver control example](../examples/solvers/solver_control.ipynb).

::: gropt.SolverGroptSDMM
    options:
        show_root_heading: false
        inherited_members: true
        summary: true
        heading_level: 4

## SolverOSQP

::: gropt.SolverOSQP
    options:
        show_root_heading: false
        summary: true
        heading_level: 4

## Convergence Plots

With `solver.extra_debug = True`, `solver.get_debug()` returns per-iteration histories, which these functions plot.

::: gropt.debug_plotting.plot_convergence

::: gropt.debug_plotting.scrub_hist_X

::: gropt.debug_plotting.scrub_hist_X_widget

::: gropt.debug_plotting.scrub_hist_X_mpl

## Equilibration

::: gropt.equilibrate

::: gropt.estimate_row_col_norms

::: gropt.rescale_eq_vecs

::: gropt.estimate_spec_norm

::: gropt.estimate_individual_spec_norm
