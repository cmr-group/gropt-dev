## Defining the problem

The optimization problem is contained in a `GroptParams` object.  This contains descriptors of the waveform desired (size, raster time, number of axes, inversions, etc.).  The problem can be defined directly using these settings, or with helper functions such as `vec_init_simple` and `diff_init`.

For more information see [here](problem.md).

## Adding constraints and objectives

The `GroptParams` object has a range of functions that add constraints (gradient amplitude, slew rate, moments, b-value, PNS, eddy currents, etc.) and objectives.

Available constraints can be seen [here](constraints.md).

## Solving

Once the problem has been defined, and constraints and/or objective functions added, solve with the convenience `gropt.solve()` function, or by creating a solver object (`SolverGroptSDMM` or `SolverOSQP`) for finer control.

See [here](solvers.md) for details.