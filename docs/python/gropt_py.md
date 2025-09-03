## Defining the problem

The optimization problem is contained in a `GroptParams` object.  This contains descriptors of the waveform desired (size, raster time, nummber of axis, inversions,...).  The problem can defined directly using these settings (see TODO: Add Example), or with helper functions, such as:


For more information see [here](problem.md)

## Adding constraints and objectives

The `GroptParams` object has a range of functions that add constraints.

Available constraints can be seen [here](constraints.md)

## Solving

Once the problem has been defined, and constraints and/or objective functions added, we solve for a solution with either the `GroptParams.solve()` function, which solves the problem with default options, or by specifying a solver and corresponding options as described here [here](solvers.md)