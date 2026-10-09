Every constraint is enforced by the ADMM solver, which counts it as met once it is within its tolerance.
Equality-type constraints (moments, eddy currents, concomitant balance) can instead take `project=True`, which
enforces them exactly on every iteration by a null-space projection; this is usually faster as well.

`weight_mod` sets a constraint's initial ADMM weight (default 1.0). The solver adapts the weights while it runs
(per-operator BB and global reweighting, see [Solving](solvers.md)), so the starting weights mostly change how
quickly a solve converges rather than its result.

## Gradient Amplitude

::: gropt.GroptParams
    options:
        members:
            - add_gmax
            - add_gmax_vec
        heading_level: 3
        show_root_heading: false

## Slew Rate

::: gropt.GroptParams
    options:
        members:
            - add_smax
            - add_smax_vec
        heading_level: 3
        show_root_heading: false

## Moments

::: gropt.GroptParams
    options:
        members:
            - add_moment
        heading_level: 3
        show_root_heading: false

## b-value

The b-value can be a constraint, or with `as_objective=True` an objective to maximize (see also [Diffusion](diffusion.md)).

::: gropt.GroptParams
    options:
        members:
            - add_bvalue
            - normalize_obj
            - add_diff_basin
        heading_level: 3
        show_root_heading: false

## PNS (SAFE-model)

PNS parameters can be loaded from an .asc file using the `gropt.readasc.asc_to_safe` function, which will create a python dictionary that can be used as input to the SAFE constraints. They can also be built by hand with `gropt.make_safe_params`, or drawn at random (realistic, but not a real scanner) with `gropt.get_random_safe_params`. Alternatively, no parameters can be entered to use some default demonstration values.  Finally, to get a SAFE prediction, the `gropt.get_SAFE` function can be used.

With several gradient axes, `add_SAFE` limits the root-sum-square of the per-axis stimulation. The `safe_*` properties change how the model is discretized and solved; set them before `add_SAFE`. See the [PNS example](../examples/hardware/pns.ipynb).

::: gropt.readasc.asc_to_safe
    options:
        heading_level: 3

::: gropt.readasc.make_safe_params
    options:
        heading_level: 3

::: gropt.readasc.get_random_safe_params
    options:
        heading_level: 3

::: gropt.GroptParams
    options:
        members:
            - add_SAFE
            - add_SAFE_vec
            - safe_alpha_exact
            - safe_lifted
            - safe_eps
            - safe_signed13
        heading_level: 3
        show_root_heading: false

::: gropt.get_SAFE
    options:
        heading_level: 3

## Eddy Currents

::: gropt.GroptParams
    options:
        members:
            - add_eddy
            - eq_proj_solver
            - eq_proj_rcond
        heading_level: 3
        show_root_heading: false

## Concomitant Fields

::: gropt.GroptParams
    options:
        members:
            - add_concomitant
        heading_level: 3
        show_root_heading: false

## Total Variation

::: gropt.GroptParams
    options:
        members:
            - add_TV
        heading_level: 3
        show_root_heading: false

## Objectives

With only constraints, the solver returns the first feasible waveform it reaches. An objective chooses among the feasible waveforms: `add_obj_identity` (minimum energy) and the b-value objective above. See the [regularization example](../examples/objective_functions/regularizations.ipynb).

::: gropt.GroptParams
    options:
        members:
            - add_obj_identity
        heading_level: 3
        show_root_heading: false
