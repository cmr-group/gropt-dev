## Defining the problem

The optimization problem is contained in a `GroptParams` object.  This contains descriptors of the waveform desired (size, raster time, number of axes, inversions, etc.).  The problem can be defined directly using these settings, or with helper functions such as `vec_init_simple` and `diff_init`.

For more information see [here](problem.md).

## Adding constraints and objectives

The `GroptParams` object has a range of functions that add constraints (gradient amplitude, slew rate, moments, b-value, PNS, eddy currents, etc.) and objectives.

Available constraints can be seen [here](constraints.md).

## Solving

Once the problem has been defined, and constraints and/or objective functions added, solve with the convenience `gropt.solve()` function, or by creating a solver object (`SolverGroptSDMM` or `SolverOSQP`) for finer control.

See [here](solvers.md) for details.

GrOpt solves for a fixed number of samples. To find the shortest waveform that meets the constraints, search over `N`: if a duration is feasible, every longer one is too, so a bisection takes a few solves (see the [flow encoding example](../examples/moments/flow_encoding.ipynb)).

## Diffusion

`gropt.diffusion` builds spin-echo diffusion problems from a `DiffParams` description and solves them with tuned settings, with TE searches and sweeps on top. See [here](diffusion.md).

## Utilities

Plotting, waveform metrics (moments, b-value), resampling to the hardware raster and logging are described [here](utilities.md).
