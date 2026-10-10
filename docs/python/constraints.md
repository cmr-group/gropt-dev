Every constraint is enforced by the ADMM solver, which counts it as met once it is within its tolerance.
Equality-type constraints (moments, eddy currents, concomitant balance) can instead take `project=True`, which
enforces them exactly on every iteration by a null-space projection; this is usually faster as well.

`weight_mod` sets a constraint's initial ADMM weight (default 1.0). The solver adapts the weights while it runs
(per-operator BB and global reweighting, see [Solving](solvers.md)), so the starting weights mostly change how
quickly a solve converges rather than its result.

All of these are methods or properties of `GroptParams`.

## Gradient Amplitude

::: gropt.GroptParams.add_gmax

::: gropt.GroptParams.add_gmax_vec

## Slew Rate

::: gropt.GroptParams.add_smax

::: gropt.GroptParams.add_smax_vec

## Moments

::: gropt.GroptParams.add_moment

## b-value

The b-value can be a constraint, or with `as_objective=True` an objective to maximize (see also [Diffusion](diffusion.md)).

::: gropt.GroptParams.add_bvalue

::: gropt.GroptParams.normalize_obj

::: gropt.GroptParams.add_diff_basin

## PNS (SAFE-model)

PNS parameters can be loaded from an .asc file using the `gropt.readasc.asc_to_safe` function, which will create a python dictionary that can be used as input to the SAFE constraints. They can also be built by hand with `gropt.make_safe_params`, or drawn at random (realistic, but not a real scanner) with `gropt.get_random_safe_params`. Alternatively, no parameters can be entered to use some default demonstration values.  Finally, to get a SAFE prediction, the `gropt.get_SAFE` function can be used.

With several gradient axes, `add_SAFE` limits the root-sum-square of the per-axis stimulation. The `safe_*` properties change how the model is discretized and solved; set them before `add_SAFE`. See the [PNS example](../examples/hardware/pns.ipynb).

::: gropt.GroptParams.add_SAFE

::: gropt.GroptParams.add_SAFE_vec

::: gropt.get_SAFE

::: gropt.readasc.asc_to_safe

::: gropt.readasc.make_safe_params

::: gropt.readasc.get_random_safe_params

::: gropt.GroptParams.safe_alpha_exact

::: gropt.GroptParams.safe_lifted

::: gropt.GroptParams.safe_eps

::: gropt.GroptParams.safe_signed13

## Eddy Currents

::: gropt.GroptParams.add_eddy

::: gropt.GroptParams.eq_proj_solver

::: gropt.GroptParams.eq_proj_rcond

## Concomitant Fields

::: gropt.GroptParams.add_concomitant

## Total Variation

::: gropt.GroptParams.add_TV

## Objectives

With only constraints, the solver returns the first feasible waveform it reaches. An objective chooses among the feasible waveforms: `add_obj_identity` (minimum energy) and the b-value objective above. See the [regularization example](../examples/objective_functions/regularizations.ipynb).

::: gropt.GroptParams.add_obj_identity
